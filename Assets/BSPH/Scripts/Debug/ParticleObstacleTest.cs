using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

[RequireComponent(typeof(MeshRenderer))]
[RequireComponent(typeof(MeshFilter))]
public class ParticleObstacleTest : MonoBehaviour
{   
    // This enum system ensure that the obstacle will or will not change in the duration of the simulation
    // Those who are considered `static` will not change their bounds or world-scale vertices/normals.
    // Only those who are `dynamic` will change.
    public enum ObstacleType {
        Static,
        Dynamic
    }
    
    // Common structure for the determination of a triangle in a given mesh.
    // Includes the indices of the vertexes (stored inside the `_vertices` array), ... 
    // ... the world-scale centroid, and the world-scale normal vector
    public class ParticleTriangle {
        public int[] vertexIndices;
        public float3 c;
        public float3 n;
    }

    // The enum `ObstacleType`, in value form and accessible via the Inspector
    [SerializeField] private ObstacleType _obstacleType = ObstacleType.Static;

    // DEBUG ELEMENT: the particle that is used for debugging the system
    [SerializeField] private Transform _testParticle = null;

    // privately reference the shared mesh used by the mesh filter for this obstacle
    private Mesh _mesh;
    // The boundaries of the obstacle in world-scale. Extracted from the mesh renderer, not from any vertex calculation
    [SerializeField] private float[] _bounds;
    // The world-scale vertices of the mesh.
    [SerializeField] private float3[] _vertices;
    // List to store the triangles of the simulation.
    private List<ParticleTriangle> _triangles = new List<ParticleTriangle>();

    // Is the debug particle within the boundaries of this obstacle?
    [SerializeField, ReadOnly] private bool _isIntersecting = false;
    // Counter for the intersection calculation. If the counter > 0, then the partice is inside the mesh
    [SerializeField, ReadOnly] private int _counter = 0;
    // DEBUG ELEMENT: stores all projections of the current debug particle so that we can render them via Gizmos
    private List<Vector3> projections = new List<Vector3>();
    // DEBUG ELEMENT: stores the closest projection and is used for the raycast direction.
    private Vector3 closestPoint;

    [SerializeField]
    private ComputeShader _obstacleShader = null;

    void OnDrawGizmos() {
        if (_bounds.Length != 6) return;
        Vector3 b = new Vector3(_bounds[3]-_bounds[0], _bounds[4]-_bounds[1], _bounds[5]-_bounds[2]);
        Vector3 c = new Vector3(_bounds[3]+_bounds[0], _bounds[4]+_bounds[1], _bounds[5]+_bounds[2]) / 2f;
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(c,b);
        
        Gizmos.color = Color.blue;
        Gizmos.DrawRay(_testParticle.position, closestPoint - _testParticle.position);

        if (_isIntersecting) {
            if (projections.Count > 0) {
                Gizmos.color = Color.yellow;
                foreach(Vector3 p in projections) {
                    Gizmos.DrawSphere(p,0.05f);

                }
            }
        }
        
    }

    // == This is called when the scene first starts up. ==
    void Start() {
        // Get a reference to the mesh that we will be using throughout the simulation
        _mesh = GetComponent<MeshFilter>().sharedMesh;

        // We want to initialize the kernels that refer to the various functions inside the obstacle shader, if it exists.
        // This is necessary as a first step because we'll be updating the buffers and variables at a later point and will need access to those kernels.
        if (_obstacleShader != null) {   
            InitializeKernels();
            UpdateParticleVariables();
        }

        // We need to calculate the world-space triangles and vertices 
        // We don't insert this into an `if (_obstacleShader != null)` check - these functions individaully will do the task for us
        CalculateTrianglesAndVertices();
        UpdateBounds();

        closestPoint = (_testParticle != null) ? _testParticle.position : transform.position;
        transform.hasChanged = false;
    }

    private int _CHECK_IN_BOUNDS_KERNEL;
    private int _FIND_CLOSEST_POINT_KERNEL;
    private void InitializeKernels() {
        _CHECK_IN_BOUNDS_KERNEL = _obstacleShader.FindKernel("CheckInBounds");
        _FIND_CLOSEST_POINT_KERNEL = _obstacleShader.FindKernel("FindClosestPoint");
    }

    

    // THIS IS PURELY FOR DEBUGGING PURPOSES
    private void UpdateParticleVariables() {
        if (DEBUG_PARTICLES_BUFFER == null) {
            DEBUG_PARTICLES_BUFFER = new ComputeBuffer(1, sizeof(float)*3);
            _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "particles", DEBUG_PARTICLES_BUFFER);
            _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "particles", DEBUG_PARTICLES_BUFFER);
        }

        _obstacleShader.SetInt("numParticles", 1);
        float3[] debug_particle_array = new float3[1] {
            new(_testParticle.position.x, _testParticle.position.y, _testParticle.position.z)
        };
        DEBUG_PARTICLES_BUFFER.SetData(debug_particle_array);

    }

    // Update is called once per frame
    void Update() {
        if (_obstacleType == ObstacleType.Dynamic && transform.hasChanged) {
            UpdateBounds();
            CalculateTrianglesAndVertices();
            transform.hasChanged = false;
        }
        if (_testParticle ==  null) {
            closestPoint = transform.position;
            return;
        }
        
        if (_obstacleShader != null) {
            // We'll be optimizing using the GPU
            GPU_UPDATE();
        } else {
            // We don't have a GPU shader for calculations. We'll fall back on CPU
            CPU_Update();
        }
    }

    private void GPU_UPDATE() {
        UpdateParticleVariables();
        _isIntersecting = CheckIfInBounds();

    }

    private void CPU_Update() {
        _isIntersecting = CheckIfInBounds();
        if (_isIntersecting) {
            FindClosestPoint();
            if (CheckIfIntersecting()) Debug.Log("Intersecting!");
        } else {
            closestPoint = _testParticle.position;
        }
    }

    private void CalculateTrianglesAndVertices() {
        // Create a temp list to store all vertices
        List<float3> fixedVerts = new List<float3>();
        List<ParticleTriangle> fixedTriangles = new List<ParticleTriangle>();

        // initialize vertices and triangles
        var vs = _mesh.vertices;
        var ts = _mesh.triangles;
        Vector4 rot = new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);
        // Initialize vertices 1, 2, and 3, as well as new `particletriangle`, for the loop
        Vector3 v1, v2, v3;
        float3 v1f, v2f, v3f;
        ParticleTriangle triangle;
        // iterate through all triangles of mesh
        for(int t = 0; t < ts.Length; t+=3) {
            // Initialize new triangle
            triangle = new ParticleTriangle();
            triangle.vertexIndices = new int[3];
            // Grab the vertices
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
            v1f = new(v1.x, v1.y, v1.z);
            v2f = new(v2.x, v2.y, v2.z);
            v3f = new (v3.x, v3.y, v3.z);
            // Add to `fixedVerts` if vertex not present.
            if (!fixedVerts.Contains(v1f)) fixedVerts.Add(v1f);
            if (!fixedVerts.Contains(v2f)) fixedVerts.Add(v2f);
            if (!fixedVerts.Contains(v3f)) fixedVerts.Add(v3f);
            // Add to `vertexIndices` of current triangle
            triangle.vertexIndices[0] = fixedVerts.IndexOf(v1f);
            triangle.vertexIndices[1] = fixedVerts.IndexOf(v2f);
            triangle.vertexIndices[2] = fixedVerts.IndexOf(v3f);
            // Calculate centroid based on average of v1,v2,v3
            triangle.c = (v1f + v2f + v3f) / 3f;
            // Calcualte normal based on normals of v1,v2,v3
            Vector3 normDir = Vector3.Cross(v2 - v1, v3 - v1);
            triangle.n = new(normDir.x, normDir.y, normDir.z);
            // Add triangle to list of triangles we have
            fixedTriangles.Add(triangle);
        }

        // Last step: Set `_triangles` as `fixedTriangles`, and convert `fixedfVerts` into array
        if (_obstacleShader != null) {
            UpdateMeshVariables(fixedTriangles, fixedVerts);
        }
        else {
            _triangles = fixedTriangles;
            _vertices = fixedVerts.ToArray();
        }
    }

    public ComputeBuffer TRIANGLE_BUFFER = null;
    public ComputeBuffer VERTICES_BUFFER = null;
    public ComputeBuffer DEBUG_PARTICLES_BUFFER = null;
    private void UpdateMeshVariables(List<ParticleTriangle> tris, List<float3> verts) {
        // Set variables
        _obstacleShader.SetInt("numTriangles", tris.Count);
        _obstacleShader.SetInt("numVertices", verts.Count);

        if (TRIANGLE_BUFFER != null) TRIANGLE_BUFFER.Release();
        TRIANGLE_BUFFER = new ComputeBuffer(tris.Count, sizeof(int)*3 + sizeof(float)*6);
        TRIANGLE_BUFFER.SetData(tris.ToArray());

        if (VERTICES_BUFFER != null) VERTICES_BUFFER.Release();
        VERTICES_BUFFER = new ComputeBuffer(verts.Count, sizeof(float)*3);
        VERTICES_BUFFER.SetData(verts.ToArray());

        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "triangles", TRIANGLE_BUFFER);
        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "vertices", VERTICES_BUFFER);
    }

    private void UpdateBounds() {
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
        if (_obstacleShader != null) UpdateBoundsVariables(_bounds);
    }

    public ComputeBuffer BOUNDS_BUFFER = null;
    private int _NUM_BLOCKS = 1024;
    private int _NUM_BLOCKS_PARTICLES;
    private void UpdateBoundsVariables(float[] bs) {
        if (BOUNDS_BUFFER == null) {
            BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
            _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "bounds", BOUNDS_BUFFER);
        }
        BOUNDS_BUFFER.SetData(bs);
        _NUM_BLOCKS_PARTICLES = Mathf.CeilToInt(1f / (float)_NUM_BLOCKS);
    }

    private bool CheckIfInBounds() {
            return (
                _bounds[0] <= _testParticle.position.x 
                && _bounds[1] <= _testParticle.position.y
                && _bounds[2] <= _testParticle.position.z
                && _bounds[3] >= _testParticle.position.x
                && _bounds[4] >= _testParticle.position.y
                && _bounds[5] >= _testParticle.position.z
            );
    }

    private void FindClosestPoint() {
        Vector3 v1, v2, v3, c, n;
        Vector3 projection;

        closestPoint = _testParticle.position;
        float dist, closestDistance = Mathf.Infinity;

        for(int ti = 0; ti < _triangles.Count; ti++) {
            v1 = _vertices[_triangles[ti].vertexIndices[0]];
            v2 = _vertices[_triangles[ti].vertexIndices[1]];
            v3 = _vertices[_triangles[ti].vertexIndices[2]];
            c = new Vector3(_triangles[ti].c[0], _triangles[ti].c[1], _triangles[ti].c[2]);
            n = new Vector3(_triangles[ti].n[0], _triangles[ti].n[1], _triangles[ti].n[2]).normalized;

            GetRayProjectionOntoPlane(
                _testParticle.position, 
                Mathf.Sign(Vector3.Dot(n, c - _testParticle.position)) * n, 
                n, 
                _vertices[_triangles[ti].vertexIndices[0]], 
                out projection
            );
            if (CheckIfPointInTriangle(projection, v1, v2, v3, n) == 1) {
                dist = Vector3.Distance(_testParticle.position, projection);
                if (dist < closestDistance) {
                    closestDistance = dist;
                    closestPoint = projection;
                }
            }
        }
    }

    private bool CheckIfIntersecting() {        
        int isValid;
        Vector3 v1, v2, v3, c, n;
        Vector3 projection;

        _counter = 0;
        projections = new List<Vector3>();
        Vector4 rot = new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);

        Vector3 raycast = closestPoint - _testParticle.position;
        
        foreach(ParticleTriangle t in _triangles) {
            v1 = _vertices[t.vertexIndices[0]];
            v2 = _vertices[t.vertexIndices[1]];
            v3 = _vertices[t.vertexIndices[2]];
            c = new Vector3(t.c[0], t.c[1], t.c[2]);
            n = new Vector3(t.n[0], t.n[1], t.n[2]).normalized;

            isValid = GetRayProjectionOntoPlane(_testParticle.position, raycast, n, v1, out projection);
            if (isValid == 0) continue;
        
            if (CheckIfPointInTriangle(projection, v1, v2, v3, n) == 1) {
                _counter += (Vector3.Dot(raycast,n) <= 0) ? -1 : 1;
                projections.Add(projection);
            }
        }
        return _counter > 0;
    }

    private int GetRayProjectionOntoPlane(Vector3 rayOrigin, Vector3 rayDirection, float3 normal, Vector3 refPoint, out Vector3 projection) {
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

    private static int CheckIfPointInTriangle(Vector3 point, Vector3 v1, Vector3 v2, Vector3 v3, Vector3 normal) {
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

    private void OnDestroy() {
        if (BOUNDS_BUFFER != null) BOUNDS_BUFFER.Release();
        if (TRIANGLE_BUFFER != null) TRIANGLE_BUFFER.Release();
        if (VERTICES_BUFFER != null) TRIANGLE_BUFFER.Release();
    }
}
