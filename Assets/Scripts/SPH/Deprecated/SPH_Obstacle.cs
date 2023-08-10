using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(MeshFilter))]
public class SPH_Obstacle : MonoBehaviour
{
    
    // This enum system ensure that the obstacle will or will not change in the duration of the simulation
    // Those who are considered `static` will not change their bounds or world-scale vertices/normals.
    // Only those who are `dynamic` will change.
    public enum ObstacleType { Static, Dynamic }

    // THis enum system ensures that the obstacle will cause a reaction to the particles in some way.
    public enum ObstacleEffect { Reflect, Stick }

    // Common structure for defining the transform of this obstacle, in a way that the GPU can understand
    public struct ObsTransform {
        public float3 position;         // The world position of the obstacle
        public float4 rotation;         // The world rotation of the obstacle
        public float3 scale;            // The world scale of the obstacle
    }

    // Common structure for the determination of a triangle in a given mesh.
    // Includes the indices of the vertexes (stored inside the `_vertices` array), ... 
    // ... the world-scale centroid, and the world-scale normal vector, and the distance of this plane from the origin
    [System.Serializable]
    public struct ParticleTriangle {
        public int3 vertexIndices;
        public float3 c;
        public float3 n;
        public float d;
    }

    // Common structure for the detenmrination of a projection towards an obstalce object
    [System.Serializable]
    public struct Projection {
        public int triangleIndex;
        public float3 position;
        public float d;
        public int count;
        public float3 check;
    }

    [Header("== BASE CLASS VARIABLES ==")]

    [SerializeField] protected internal ObstacleType _obstacleType = ObstacleType.Static;
    public ObstacleType obstacleType => _obstacleType;
    [SerializeField] protected internal ObstacleEffect _obstacleEffect = ObstacleEffect.Reflect;
    public ObstacleEffect obstacleEffect => _obstacleEffect;

    [SerializeField] protected internal float _particle_radius;

    // This might be for debugging purposes, but we'll store a list of Transforms as our debug "particles"
    [SerializeField] protected internal List<Transform> _debugParticles = new List<Transform>();

    public void CalculateTrianglesAndVertices(out float3[] verts, out ParticleTriangle[] tris) {
        // Grab a reference to the shared mesh from the MeshFilter component. Note that all vertices, normals, etc. are in local space
        Mesh mesh = GetComponent<MeshFilter>().sharedMesh;

        // initialize vertices and triangles
        var vs = mesh.vertices;
        var ts = mesh.triangles;
        Vector4 rot = new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);

        // Create a temp list to store all vertices
        verts = new float3[vs.Length];
        tris = new ParticleTriangle[ts.Length / 3];
        
        // Initialize vertices 1, 2, and 3, as well as new `particletriangle`, for the loop
        Vector3 v1, v2, v3;
        float3 v1f, v2f, v3f;
        Plane plane;
        ParticleTriangle triangle;
        
        // iterate through all triangles of mesh
        // The magic of this system is that we can derive centroid and normal purely based off of vertices alone
        // We don't even need to look at the normals at all!
        // Therefore, in concept, all we need to do is store the world space positions of each vertex.
        // We can derive normals from just those. We'll store the normals as well, to optimize on calculations
        // ... As well as centroids.
        for(int t = 0; t < ts.Length; t+=3) {
            // Initialize new triangle
            triangle = new ParticleTriangle();
            // Grab the vertices and convert them into world scale
            v1 = transform.TransformPoint(vs[ts[t]]);
            v2 = transform.TransformPoint(vs[ts[t+1]]);
            v3 = transform.TransformPoint(vs[ts[t+2]]);
            /*
            v1 =  PlaneObstacleManager.LocalPointToWorldPoint(
                transform.position, 
                rot, 
                transform.lossyScale, 
                vs[ts[t]]
            );
            v2 = PlaneObstacleManager.LocalPointToWorldPoint(
                transform.position, 
                rot, 
                transform.lossyScale, 
                vs[ts[t+1]]
            );
            v3 = PlaneObstacleManager.LocalPointToWorldPoint(
                transform.position, 
                rot, 
                transform.lossyScale, 
                vs[ts[t+2]]
            );
            */
            // Form float3 versions of v1, v2, and v3
            v1f = new(v1.x, v1.y, v1.z);
            v2f = new(v2.x, v2.y, v2.z);
            v3f = new (v3.x, v3.y, v3.z);
            // Add to `verts` if vertex not present.
            verts[ts[t]] = v1f;
            verts[ts[t+1]] = v2f;
            verts[ts[t+2]] = v3f;
            // Add to `vertexIndices` of current triangle
            triangle.vertexIndices = new(ts[t],ts[t+1],ts[t+2]);
            // Calculate centroid based on average of v1,v2,v3
            triangle.c = (v1f + v2f + v3f) / 3f;
            // Calculate normal based on normals of v1,v2,v3
            Vector3 normDir = Vector3.Cross(v2 - v1, v3 - v1);
            triangle.n = new(normDir.x, normDir.y, normDir.z);
            plane = new Plane(v1,v2,v3);
            triangle.d = plane.GetDistanceToPoint(Vector3.zero);
            // Add triangle to list of triangles we have
            tris[t/3] = triangle;
        }

        // since `verts` and `tris` are outgoing variables, nothing else needs to be done
    }

    public void UpdateBounds(out float[] _bounds) {
        // We can get world-scale bounds via renderer.boudns
        var bounds = GetComponent<MeshRenderer>().bounds;
        // To get 6-float bounds, we need to get bounds.min and bounds.max, which are Vector3s
        _bounds = new float[6] {
            bounds.min.x,
            bounds.min.y,
            bounds.min.z,
            bounds.max.x,
            bounds.max.y,
            bounds.max.z
        };
    }

    public static int GetRayProjectionOntoPlane(Vector3 rayOrigin, Vector3 rayDirection, float3 normal, Vector3 refPoint, out Vector3 projection) {
        projection = rayOrigin;
        Vector3 normalV3 = new Vector3(normal[0], normal[1], normal[2]);
        float ray_n_dot =  Vector3.Dot(normalV3, rayDirection.normalized);
        if (ray_n_dot == 0) return 0;

        float d = -Vector3.Dot(normalV3, refPoint);
        float dist_from_point_to_plane = -(Vector3.Dot(normalV3, rayOrigin) + d) / ray_n_dot;
            
        if (dist_from_point_to_plane < 0) return 0;
        
        projection = rayOrigin + dist_from_point_to_plane * rayDirection.normalized;
        return 1;
    }
    public static int CheckIfPointInTriangle(Vector3 point, Vector3 v1, Vector3 v2, Vector3 v3, Vector3 normal) {
        Vector3 edge0 = v2 - v1;
        Vector3 edge1 = v3 - v2;
        Vector3 edge2 = v1 - v3;
        Vector3 C0 = point - v1;
        Vector3 C1 = point - v2;
        Vector3 C2 = point - v3;
        if (
            Vector3.Dot(normal, Vector3.Cross(edge0, C0)) > 0 
            && Vector3.Dot(normal, Vector3.Cross(edge1, C1)) > 0 
            && Vector3.Dot(normal, Vector3.Cross(edge2, C2)) > 0
        ) return 1;
        return 0;
    }
}
