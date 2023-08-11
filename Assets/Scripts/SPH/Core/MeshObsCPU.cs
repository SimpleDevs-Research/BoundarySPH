using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using Unity.Mathematics;
using UnityEditor;
using OP = ObstaclePrimitives.Classes;

public class MeshObsCPU : MonoBehaviour
{
    // This one stays in the CPU and simply stores basic info about the obstacle. 
    // This one is NOT sent over to the 
    [System.Serializable]
    public class TestObstacle {
        public Transform obstacle;
        public bool show_first_triangle_only = false;

        public bool show_face_normals = true;
        public bool show_vertex_normals = true;
        public bool show_3D_edge_normals = true;
        public bool show_2D_edge_normals = true;
        public bool show_centroids = true;
        public bool show_bounds = true;
        public bool show_triangle_bounds = true;
        public bool show_gizmos => show_face_normals 
            || show_vertex_normals 
            || show_3D_edge_normals 
            || show_2D_edge_normals
            || show_centroids
            || show_bounds
            || show_triangle_bounds
            || show_triangle_bounds;
    }
    
    public List<TestObstacle> obstacles;

    [SerializeField] private List<OP.ObstacleStatic> obstacles_static = new List<OP.ObstacleStatic>();
    [SerializeField] private List<OP.ObstacleDynamic> obstacles_dynamic = new List<OP.ObstacleDynamic>();
    [SerializeField] private List<OP.VertexStatic> vertices_static = new List<OP.VertexStatic>();
    [SerializeField] private List<OP.VertexDynamic> vertices_dynamic = new List<OP.VertexDynamic>();
    [SerializeField] private List<OP.TriangleStatic> triangles_static = new List<OP.TriangleStatic>();
    [SerializeField] private List<OP.TriangleDynamic> triangles_dynamic = new List<OP.TriangleDynamic>();
    [SerializeField] private List<OP.EdgeStatic> edges_static = new List<OP.EdgeStatic>();
    [SerializeField] private List<float3> edges_dynamic = new List<float3>();

    [SerializeField] private bool drawGizmos = false;

    void OnDrawGizmos() {
        if (!drawGizmos) return;
        OP.ObstacleStatic o_static;
        OP.ObstacleDynamic o_dynamic;
        OP.VertexDynamic v_dynamic;
        for(int i = 0; i < obstacles_static.Count; i++) {
            if (!obstacles[i].show_gizmos) continue;
            o_static = obstacles_static[i];
            o_dynamic = obstacles_dynamic[i];
            // Render vertices
            for(int vi = o_static.vs[0]; vi < o_static.vs[0] + o_static.vs[1]; vi++) {
                v_dynamic = vertices_dynamic[vi];
                Vector3 vp = new Vector3(v_dynamic.position[0], v_dynamic.position[1], v_dynamic.position[2]);
                Vector3 vn = new Vector3(v_dynamic.normal[0],v_dynamic.normal[1],v_dynamic.normal[2]).normalized;
                Gizmos.color = Color.red;
                Gizmos.DrawSphere(vp,1f);
                Handles.color = Color.blue;
                Handles.DrawLine(vp, vp + vn*5f, 3);
            }
            // Render triangles
            for(int ti = o_static.ts[0]; ti < o_static.ts[0] + o_static.ts[1]; ti++) {
                // Draw the center
                Gizmos.color = Color.red;
                Gizmos.DrawSphere(triangles_dynamic[ti].center,0.75f);
                // Draw the 3D face normal
                Handles.color = Color.blue;
                Handles.DrawLine(
                    triangles_dynamic[ti].center,
                    triangles_dynamic[ti].center + triangles_dynamic[ti].normal*10f, 
                    3
                );
                // Draw the 2D edge normals
                Handles.color = Color.yellow;
                Handles.DrawLine(
                    triangles_dynamic[ti].center,
                    triangles_dynamic[ti].center + triangles_dynamic[ti].v1v2n*3f,
                    3
                );
                Handles.DrawLine(
                    triangles_dynamic[ti].center,
                    triangles_dynamic[ti].center + triangles_dynamic[ti].v2v3n*3f,
                    3
                );
                Handles.DrawLine(
                    triangles_dynamic[ti].center,
                    triangles_dynamic[ti].center + triangles_dynamic[ti].v1v3n*3f,
                    3
                );
                
                Gizmos.color = Color.grey;
                Gizmos.DrawWireCube(
                    (triangles_dynamic[ti].upperBound + triangles_dynamic[ti].lowerBound)/2f,
                    triangles_dynamic[ti].upperBound - triangles_dynamic[ti].lowerBound
                );
                
            }
            // Render the edge normals
            Handles.color = Color.blue;
            for(int ei = o_static.es[0]; ei < o_static.es[0] + o_static.es[1]; ei++) {
                o_static = obstacles_static[edges_static[ei].obstacleIndex];
                float3 v1 = vertices_dynamic[o_static.vs[0] + edges_static[ei].vertices[0]].position;
                float3 v2 = vertices_dynamic[o_static.vs[0] + edges_static[ei].vertices[1]].position;
                float3 m = (v1+v2)/2f;
                Vector3 midpoint = new Vector3(m[0],m[1],m[2]);
                float3 n = edges_dynamic[ei];
                Vector3 normal = new Vector3(n[0],n[1],n[2]);
                Handles.DrawLine(midpoint, midpoint + normal, 3);
            }
            // Render bounds
            Gizmos.color = Color.black;
            Gizmos.DrawWireCube(
                (o_dynamic.upperBound + o_dynamic.lowerBound)/2f, 
                o_dynamic.upperBound - o_dynamic.lowerBound
            );
        }
    }

    // Start is called before the first frame update
    void Start() {
        PreprocessObstacles();
    }

    public void PreprocessObstacles() {
        obstacles_static = new List<OP.ObstacleStatic>();
        obstacles_dynamic = new List<OP.ObstacleDynamic>();
        vertices_static = new List<OP.VertexStatic>();
        vertices_dynamic = new List<OP.VertexDynamic>();
        triangles_static = new List<OP.TriangleStatic>();
        triangles_dynamic = new List<OP.TriangleDynamic>();
        edges_static = new List<OP.EdgeStatic>();
        edges_dynamic = new List<float3>();
        for(int i = 0; i < obstacles.Count; i++) {
            TestObstacle obstacle = obstacles[i];
            PreprocessObstacle(obstacle,i);
            obstacle.obstacle.hasChanged = true;
        }
    }

    public void UpdateObstacles() {
        Transform t;
        bool needsUpdating = false;
        for(int i = 0; i < obstacles.Count; i++) {
            // Check if the transform has been changed in any way
            t = obstacles[i].obstacle; 
            if (t.hasChanged) {
                //  We need to update the associated TestObjectDynamic
                obstacles_dynamic[i].position = new(t.position.x, t.position.y, t.position.z);
                obstacles_dynamic[i].rotation = new(t.rotation.x, t.rotation.y, t.rotation.z, t.rotation.w);
                obstacles_dynamic[i].scale = new(t.lossyScale.x, t.lossyScale.y, t.lossyScale.z);
                obstacles_dynamic[i].lowerBound = obstacles_dynamic[i].position;
                obstacles_dynamic[i].upperBound = obstacles_dynamic[i].position;
                obstacles_dynamic[i].hasChanged = 1;
                needsUpdating = true;
                t.hasChanged = false;
            }
        }

        if (needsUpdating) {
            // First, let's update the vertices
            UpdateVertices(); 
            // Second, let's reset the edge normals
            UpdateEdges();
            // Thirdly, let's update the triangles
            UpdateTriangles();

            // Finally reset
            for(int i = 0; i < obstacles.Count; i++) {
                obstacles[i].obstacle.hasChanged = false;
                obstacles_dynamic[i].hasChanged = 0;
                obstacles_dynamic[i].lowerBound = obstacles_dynamic[i].lowerBound - (float3)new(1f,1f,1f);
                obstacles_dynamic[i].upperBound = obstacles_dynamic[i].upperBound + (float3)new(1f,1f,1f);
            }
        }
    }

    private void UpdateVertices() {
        // To update the vertices, we take the localPosition, convert it into world position, and reset the normal to 0,0,0
        OP.ObstacleDynamic obstacle;
        float3 localPosition;
        Vector3 worldPosition,worldNormalEnd,worldNormal;
        float[] xs;
        float[] ys;
        float[] zs;
        for(int i = 0; i < vertices_dynamic.Count; i++) {
            // Get the associated OP.ObstacleDynamic
            obstacle = obstacles_dynamic[vertices_dynamic[i].obstacleIndex];
            // Only update if the obstacle needs to be updated
            if (obstacles_dynamic[vertices_dynamic[i].obstacleIndex].hasChanged == 0) continue;
            // Get the local position
            localPosition = vertices_static[i].localPosition;
            // Convert into world position
            worldPosition = LocalPointToWorldPoint(obstacle.position,obstacle.rotation,obstacle.scale,localPosition);
            // Update the world position of the current vertex
            vertices_dynamic[i].position = worldPosition;
            // Update the world normal of the current vertex
            worldNormalEnd = LocalPointToWorldPoint(obstacle.position,obstacle.rotation,obstacle.scale,vertices_static[i].localNormalEnd);
            worldNormal = (worldNormalEnd-worldPosition).normalized;
            vertices_dynamic[i].normal = new(worldNormal.x,worldNormal.y,worldNormal.z);
            // update the bounds of the dynamic obstacle
            xs = new float[3] {obstacle.lowerBound[0], obstacle.upperBound[0], worldPosition[0]};
            ys = new float[3] {obstacle.lowerBound[1], obstacle.upperBound[1], worldPosition[1]};
            zs = new float[3] {obstacle.lowerBound[2], obstacle.upperBound[2], worldPosition[2]};
            obstacle.lowerBound = new(xs.Min(), ys.Min(), zs.Min());
            obstacle.upperBound = new(xs.Max(), ys.Max(), zs.Max());
        }
    }

    private void UpdateEdges() {
        OP.EdgeStatic e_static;
        OP.ObstacleDynamic o_dynamic;
        OP.ObstacleStatic o_static;
        Vector3 midpoint, normal;
        for(int i = 0; i < edges_static.Count; i++) {
            // References
            e_static = edges_static[i];
            o_static = obstacles_static[e_static.obstacleIndex];
            o_dynamic = obstacles_dynamic[e_static.obstacleIndex];
            // Update only if the current object has changed
            if (o_dynamic.hasChanged == 1) {
                midpoint = LocalPointToWorldPoint(o_dynamic.position,o_dynamic.rotation,o_dynamic.scale,e_static.midpoint); 
                normal = (LocalPointToWorldPoint(o_dynamic.position, o_dynamic.rotation, o_dynamic.scale,e_static.localNormalEnd) - midpoint).normalized;
                edges_dynamic[i] = new(normal.x,normal.y,normal.z);
            }
        }
    }

    private void UpdateTriangles() {
        OP.TriangleStatic t_static;
        OP.TriangleDynamic t_dynamic;
        OP.ObstacleStatic o_static;
        OP.ObstacleDynamic o_dynamic;
        OP.VertexDynamic v1f,v2f,v3f;
        Vector3 v1,v2,v3,c,nEnd;
        Vector3 n,v1v2,v1c,v1v2n,v2v3,v2c,v2v3n,v1v3,v1v3n;
        float[] xs;
        float[] ys;
        float[] zs;
        for(int i = 0; i < triangles_dynamic.Count; i++) {
            // References
            t_static = triangles_static[i];
            t_dynamic = triangles_dynamic[i];
            o_dynamic = obstacles_dynamic[t_dynamic.obstacleIndex];
            o_static = obstacles_static[t_dynamic.obstacleIndex];
            // Skip if no need to update
            if (o_dynamic.hasChanged == 0) continue;
            // Getting updated vertex positions
            v1f = vertices_dynamic[o_static.vs[0] + t_static.vertices[0]];
            v2f = vertices_dynamic[o_static.vs[0] + t_static.vertices[1]];
            v3f = vertices_dynamic[o_static.vs[0] + t_static.vertices[2]];
            v1 = new Vector3(v1f.position[0], v1f.position[1], v1f.position[2]);
            v2 = new Vector3(v2f.position[0], v2f.position[1], v2f.position[2]);
            v3 = new Vector3(v3f.position[0], v3f.position[1], v3f.position[2]);
            // Updating center
            //c = (v1+v2+v3)/3f;
            c = LocalPointToWorldPoint(o_dynamic.position, o_dynamic.rotation, o_dynamic.scale, t_static.localCenter);
            //t_dynamic.center = (v1f.position + v2f.position + v3f.position) / 3f;
            //t_dynamic.center = new(c.x,c.y,c.z);
            t_dynamic.center = new(c.x,c.y,c.z);
            // Updating normal
            //Vector3 n = Vector3.Cross(v2 - v1, v3 - v1).normalized;
            nEnd = LocalPointToWorldPoint(o_dynamic.position, o_dynamic.rotation, o_dynamic.scale, t_static.localNormalEnd);
            n = (nEnd - c).normalized;
            t_dynamic.normal = new(n.x,n.y,n.z);
            // Updating d
            t_dynamic.d = Vector3.Dot(-1f * c, n);
            // Update lowerBound and upperBound;
            xs = new float[3] {v1.x,v2.x,v3.x};
            ys = new float[3] {v1.y,v2.y,v3.y};
            zs = new float[3] {v1.z,v2.z,v3.z};
            t_dynamic.lowerBound = new(xs.Min(), ys.Min(), zs.Min());
            t_dynamic.upperBound = new(xs.Max(), ys.Max(), zs.Max());
            // updating v1v2n, v1v3n, v2v3n
            v1v2 = v2 - v1;
            v1c = c - v1;
            v1v2n = -Vector3.Normalize(v1c - (Vector3.Dot(v1c,v1v2)/Vector3.Dot(v1v2,v1v2))*v1v2);
            v2v3 = v3 - v2;
            v2c = c - v2;
            v2v3n = -Vector3.Normalize(v2c - (Vector3.Dot(v2c,v2v3)/Vector3.Dot(v2v3,v2v3))*v2v3);
            v1v3 = v3-v1;
            v1v3n = -Vector3.Normalize(v1c - (Vector3.Dot(v1c,v1v3)/Vector3.Dot(v1v3,v1v3))*v1v3);
            t_dynamic.v1v2n = new(v1v2n.x, v1v2n.y, v1v2n.z);
            t_dynamic.v1v3n = new(v1v3n.x, v1v3n.y, v1v3n.z);
            t_dynamic.v2v3n = new(v2v3n.x, v2v3n.y, v2v3n.z);
        }
    }

    void PreprocessObstacle(TestObstacle obstacle, int index) {
        
        // Before anything, we need to extract the mesh data! 
        // We also need to extract the vertex anbd triangle data from that mesh
        Mesh mesh =  obstacle.obstacle.GetComponent<MeshFilter>().sharedMesh;
        var vs = mesh.vertices;
        var ts = mesh.triangles;

        // We initialize the details for the obstacles themselves
        OP.ObstacleStatic o_static = new OP.ObstacleStatic();
        OP.ObstacleDynamic o_dynamic = new OP.ObstacleDynamic();
        o_static.index = index;
        o_dynamic.index = index;
        o_dynamic.hasChanged = 1;

        // ===== GENERATING VERTICES DATA ===== //

        // We need to condense the vertices into a smaller list
        // We do this by creating:
        //  1) `fixed_vs` : List<float3> = A list of vertices simialr to `vs` but without duplicates
        //  2) `vs_map`: int[] = A int-int mapping between `vs` and `fixed_vs`
        //  3) `vs_static` : List<OP.VertexStatic> = the static details of each vertex in `fixed_vs`
        //  4) `vs_dynamic` : List<OP.VertexDynamic> = the dynamic details of each vertex in `fixed_vs`
        List<float3> fixed_vs = new List<float3>();
        int[] vs_map = new int[vs.Length];
        List<OP.VertexStatic> vs_static = new List<OP.VertexStatic>();
        List<OP.VertexDynamic> vs_dynamic = new List<OP.VertexDynamic>();

        // We can now iterate through each verex in `vs` and filter out duplicates
        int mi;
        for(int i = 0; i < vs.Length; i++) {
            // We save the current vertex's local position
            mi = fixed_vs.IndexOf(vs[i]);
            // If the mapped index is -1, then there's no entry in `fixed_vs`
            if (mi == -1) {
                // Grab the current index, and add the local position to `fixed_vs`
                mi = fixed_vs.Count;
                fixed_vs.Add(vs[i]);
                // Initialize the corresponding entriers in `vs_static` and `vs_dynamic`
                OP.VertexStatic v_static = new OP.VertexStatic();
                OP.VertexDynamic v_dynamic = new OP.VertexDynamic();
                // We set the obstacle index to the current index
                v_static.obstacleIndex = index;
                v_dynamic.obstacleIndex = index;
                // We set the local position of the vertex
                v_static.localPosition = vs[i];
                // We set the ending point to the local position, currently. When determining the world-space normal in the update loop, we'll use this to determine the normal direction
                v_static.localNormalEnd = vs[i];
                //v_dynamic.normal = new(0f,0f,0f);
                // Add `v_static` into `vs_static` and `v_dynamic` into `vs_dynamic`
                vs_static.Add(v_static);
                vs_dynamic.Add(v_dynamic);
            }
            // Now, we merely have to update `vs_map`
            vs_map[i] = mi;
        }
        // We update `o_static` with number of filtered vertices
        // we know [0] because it's merely the current count of `vertices_static`
        o_static.vs = new(vertices_static.Count,fixed_vs.Count);
        vertices_static.AddRange(vs_static);
        vertices_dynamic.AddRange(vs_dynamic);


        // ===== GENERATING TRIANGLES AND EDGES DATA ==== //
        
        // We need to generate the data for this obstacle's triangles, and their associated edges
        // Note that an edge is associated with only two triangles
        //  1) `ts_static` : List<OP.TriangleStatic> = list of static details for each mesh triangle
        //  2) `ts_dynamic` : List<OP.TriangleDynamic> = list of dynamic details for each triangle
        //  3) `es_static` : List<OP.EdgeStatic> = list of static details for each edge
        //  4) `es_dynamic` : List<float3> = list of dynamic details for each edge (for now, the 3D normal)
        //  5) `es_static_map` : Dictionary<int2,int> = temporary dictionary linking int2 to index in `es_static`
        OP.TriangleStatic[] ts_static = new OP.TriangleStatic[ts.Length/3];
        OP.TriangleDynamic[] ts_dynamic = new OP.TriangleDynamic[ts.Length/3];
        List<OP.EdgeStatic> es_static = new List<OP.EdgeStatic>();
        List<float3> es_dynamic = new List<float3>();
        Dictionary<int2,int> es_static_map = new Dictionary<int2,int>();

        // Need to iterate through each triangle
        int ti;
        for(int i = 0; i < ts.Length; i+=3) {
            // Index for current triangle
            ti = i/3;
            // Generate a new triangle at the approriate index for both `ts_static` and `ts_dynamic`
            ts_static[ti] = new OP.TriangleStatic();
            ts_dynamic[ti] = new OP.TriangleDynamic();
            ts_static[ti].obstacleIndex = index;
            ts_dynamic[ti].obstacleIndex = index;

            // We get the corresponding indices of `vs_fixed` for this triangle. We need to store them in `ts_static[ti]`
            // In case it's confusing, know that `ts[i]` is an integer index for the original entry in `vs`.
            //  We have to get the corresponding index in `fixed_vs`. Hence `vs_map[ts[i]]`
            ts_static[ti].vertices = new(vs_map[ts[i]], vs_map[ts[i+1]], vs_map[ts[i+2]]);
            
            // We need to calculate the angles for each vertex. To do so, we need to reference the local positions of each vertex.
            float3 v1f = vs_static[ts_static[ti].vertices[0]].localPosition;
            float3 v2f = vs_static[ts_static[ti].vertices[1]].localPosition;
            float3 v3f = vs_static[ts_static[ti].vertices[2]].localPosition;
            Vector3 v1 = new Vector3(v1f[0],v1f[1],v1f[2]);
            Vector3 v2 = new Vector3(v2f[0],v2f[1],v2f[2]);
            Vector3 v3 = new Vector3(v3f[0],v3f[1],v3f[2]);
            // We can grab angles using `AngleFromVectors()`
            ts_static[ti].angles = new(
                AngleFromVectors(v3f-v1f,v2f-v1f),
                AngleFromVectors(v1f-v2f,v3f-v2f),
                AngleFromVectors(v2f-v3f,v1f-v3f)
            );

            // We can calculate the local centroid of this triangle using the average of all three vertex local positions
            ts_static[ti].localCenter = (v1f+v2f+v3f)/3f;
            // Now, we can calculate the LOCAL normal's end point for the triangle face. We do this by calculating the local centroid + normalized cross product.
            Vector3 localNormal = Vector3.Cross(v2-v1, v3-v1).normalized;
            float3 localNormalF = new(localNormal.x,localNormal.y,localNormal.z);
            ts_static[ti].localNormalEnd = ts_static[ti].localCenter + localNormalF;

            // We need to generate the normal vectors for our vertices too. We do this by, for each vertex, adding the normal weighed by the influence angle
            vs_static[ts_static[ti].vertices[0]].localNormalEnd += localNormalF * ts_static[ti].angles[0];
            vs_static[ts_static[ti].vertices[1]].localNormalEnd += localNormalF * ts_static[ti].angles[1];
            vs_static[ts_static[ti].vertices[2]].localNormalEnd += localNormalF * ts_static[ti].angles[2];
            
            // Now we have to update `es_static` and `es_dynamic`
            // The honestly crummy thing is that int2 can't be used as a key... if it's a list
            // So we used a dictionary `es_static_map` to temporarily link int2 pairs to edges in `es_static`

            int es_map_index_v1v2, es_map_index_v1v3, es_map_index_v2v3;
            // v1v2
            int2 edge1_index = (ts_static[ti].vertices[0] < ts_static[ti].vertices[1])
                ? new(ts_static[ti].vertices[0], ts_static[ti].vertices[1])
                : new(ts_static[ti].vertices[1], ts_static[ti].vertices[0]);
            // If `v1v2` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge1_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v1v2 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = index;
                e_static.vertices = edge1_index;
                e_static.triangles = new(ti,-1);
                es_static.Add(e_static);
                es_dynamic.Add(new(0f,0f,0f));
                es_static_map.Add(edge1_index,es_map_index_v1v2);
            } else {
                // The entry already exists for this edge. This is the 2nd triangle associated with this edge
                es_map_index_v1v2 = es_static_map[edge1_index];
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = es_static[es_map_index_v1v2].obstacleIndex;
                e_static.vertices = edge1_index;
                e_static.triangles = new(es_static[es_map_index_v1v2].triangles[0], ti);
                es_static[es_map_index_v1v2] = e_static;
            }
            // v1v3
            int2 edge2_index = (ts_static[ti].vertices[0] < ts_static[ti].vertices[2])
                ? new(ts_static[ti].vertices[0], ts_static[ti].vertices[2])
                : new(ts_static[ti].vertices[2], ts_static[ti].vertices[0]);
            // If `v1v3` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge2_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v1v3 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = index;
                e_static.vertices = edge2_index;
                e_static.triangles = new(ti,-1);
                es_static.Add(e_static);
                es_dynamic.Add(new(0f,0f,0f));
                es_static_map.Add(edge2_index,es_map_index_v1v3);
            } else {
                // The entry already exists for this edge. This is the 2nd triangle associated with this edge
                es_map_index_v1v3 = es_static_map[edge2_index];
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = es_static[es_map_index_v1v3].obstacleIndex;
                e_static.vertices = edge2_index;
                e_static.triangles = new(es_static[es_map_index_v1v3].triangles[0], ti);
                es_static[es_map_index_v1v3] = e_static;
            }
            // v2v3
            int2 edge3_index = (ts_static[ti].vertices[1] < ts_static[ti].vertices[2])
                ? new(ts_static[ti].vertices[1], ts_static[ti].vertices[2])
                : new(ts_static[ti].vertices[2], ts_static[ti].vertices[1]);
            // If `v2v3` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge3_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v2v3 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = index;
                e_static.vertices = edge3_index;
                e_static.triangles = new(ti,-1);
                es_static.Add(e_static);
                es_dynamic.Add(new(0f,0f,0f));
                es_static_map.Add(edge3_index,es_map_index_v2v3);
            } else {
                // The entry already exists for this edge. This is the 2nd triangle associated with this edge
                es_map_index_v2v3 = es_static_map[edge3_index];
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = es_static[es_map_index_v2v3].obstacleIndex;
                e_static.vertices = edge3_index;
                e_static.triangles = new(es_static[es_map_index_v2v3].triangles[0], ti);
                es_static[es_map_index_v2v3] = e_static;
            }

            // Add `es_map_index_v1v2` , `es_map_index_v1v3`, and `es_map_index_v2v3` to `ts_static[ti].edges
            ts_static[ti].edges = new(es_map_index_v1v2,es_map_index_v1v3,es_map_index_v2v3);
        }

        // After all that, we need to update `o_static` on the updates to triangles and edges
        // [0] is the current count of `triangles_static` and `edges_static`, respectively
        o_static.ts = new(triangles_static.Count, ts_static.Length);
        o_static.es = new(edges_static.Count, es_static.Count);
        
        // We update the global triangles array
        triangles_static.AddRange(ts_static);
        triangles_dynamic.AddRange(ts_dynamic);
        
        // We have to iterate through all edges, find their midpoints, and their localNormalEnd`s
        float3 midpoint,normal;
        for(int ei = 0; ei < es_static.Count; ei++) {
            int2 vertices = es_static[ei].vertices;
            midpoint = (vs_static[vertices[0]].localPosition + vs_static[vertices[1]].localPosition)/2f;
            // We calculate the normal by approximating the normals of each triangle assoicated with this edge. We have to be careful about edges with just one triangle (open edges).
            int2 triangles = es_static[ei].triangles;
            normal = ts_static[triangles[0]].localNormalEnd - ts_static[triangles[0]].localCenter;
            if (triangles[1] != -1) normal += (ts_static[triangles[1]].localNormalEnd - ts_static[triangles[1]].localCenter);
            // We can now update the midpoint and localNormalEnd of this edge
            es_static[ei].midpoint = midpoint;
            es_static[ei].localNormalEnd = midpoint + Unity.Mathematics.math.normalize(normal);
        }
        // And we add the edges to our global list of edges
        edges_static.AddRange(es_static);
        edges_dynamic.AddRange(es_dynamic);

        // Finally, we append `o_static` and `o_dynamic` to `obstacles_static` and `obstacles_dynamic`, respectively
        obstacles_static.Add(o_static);
        obstacles_dynamic.Add(o_dynamic);
    }

    // Update is called once per frame
    void Update() {
        UpdateObstacles();
    }

    // https://forum.unity.com/threads/whats-the-math-behind-transform-transformpoint.107401/
    public static Vector3 LocalPointToWorldPoint(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 localPoint) {
        Vector3 s = new Vector3(localPoint.x * scale.x, localPoint.y * scale.y, localPoint.z * scale.z);
        return RotMultVec3(rot, s) + pos;
    }
    public static Vector3 LocalPointToWorldPoint(Vector3 pos, Vector4 rot, Vector3 scale, float3 localPoint) {
        Vector3 s = new Vector3(localPoint[0] * scale.x, localPoint[1] * scale.y, localPoint[2] * scale.z);
        return RotMultVec3(rot, s) + pos;
    }
    public static Vector3 LocalPointToWorldPoint(float3 pf, float4 rf, float3 scale, float3 localPoint) {
        Vector3 pos = new Vector3(pf[0],pf[1],pf[2]);
        Vector4 rot = new Vector4(rf[0],rf[1],rf[2],rf[3]);
        Vector3 s = new Vector3(localPoint[0] * scale[0], localPoint[1] * scale[1], localPoint[2] * scale[2]);
        return RotMultVec3(rot, s) + pos;
    }    
    // https://answers.unity.com/questions/372371/multiply-quaternion-by-vector3-how-is-done.html
    public static Vector3 RotMultVec3(Vector4 quat, Vector3 vec){
        float num = quat.x * 2f;
        float num2 = quat.y * 2f;
        float num3 = quat.z * 2f;
        float num4 = quat.x * num;
        float num5 = quat.y * num2;
        float num6 = quat.z * num3;
        float num7 = quat.x * num2;
        float num8 = quat.x * num3;
        float num9 = quat.y * num3;
        float num10 = quat.w * num;
        float num11 = quat.w * num2;
        float num12 = quat.w * num3;
        Vector3 result;
        result.x = (1f - (num5 + num6)) * vec.x + (num7 - num12) * vec.y + (num8 + num11) * vec.z;
        result.y = (num7 + num12) * vec.x + (1f - (num4 + num6)) * vec.y + (num9 - num10) * vec.z;
        result.z = (num8 - num11) * vec.x + (num9 + num10) * vec.y + (1f - (num4 + num5)) * vec.z;
        return result;
    }
    public static Vector3 RotMultVec3(float4 quat, Vector3 vec){
        float num = quat[0] * 2f;
        float num2 = quat[1] * 2f;
        float num3 = quat[2] * 2f;
        float num4 = quat[0] * num;
        float num5 = quat[1] * num2;
        float num6 = quat[2] * num3;
        float num7 = quat[0] * num2;
        float num8 = quat[0] * num3;
        float num9 = quat[1] * num3;
        float num10 = quat[3] * num;
        float num11 = quat[3] * num2;
        float num12 = quat[3] * num3;
        Vector3 result;
        result.x = (1f - (num5 + num6)) * vec.x + (num7 - num12) * vec.y + (num8 + num11) * vec.z;
        result.y = (num7 + num12) * vec.x + (1f - (num4 + num6)) * vec.y + (num9 - num10) * vec.z;
        result.z = (num8 - num11) * vec.x + (num9 + num10) * vec.y + (1f - (num4 + num5)) * vec.z;
        return result;
    }
    public static Vector3 LocalVectorToWorldVector(Vector3 pos, Vector4 rot, Vector3 scale, Vector3 start, Vector3 end) {
        Vector3 startWorld = LocalPointToWorldPoint(pos, rot, scale, start);
        Vector3 endWorld = LocalPointToWorldPoint(pos, rot, scale, end);
        return endWorld - startWorld;
    }
    // https://github.com/HardlyDifficult/Tutorials/blob/master/Quaternions.md#361-quaternioninverse
    public static Quaternion InverseQuaternion(Quaternion rotation) {
        // Split the Quaternion component
        Vector3 vector = new Vector3(rotation.x, rotation.y, rotation.z);
        float scalar = rotation.w;
        // Calculate inverse
        vector = -vector;
        // Store results
        Quaternion inverseRotation = new Quaternion(vector.x, vector.y, vector.z, scalar);
        return inverseRotation;
    }
    public static float AngleFromVectors(Vector3 from, Vector3 to) {
        // sqrt(a) * sqrt(b) = sqrt(a * b) -- valid for real numbers
        float kEpsilonNormalSqrt = 1e-15F;
        float denominator = (float)Mathf.Sqrt(from.sqrMagnitude * to.sqrMagnitude);
        if (denominator < kEpsilonNormalSqrt) return 0f;
        float dot = Mathf.Clamp(Vector3.Dot(from, to) / denominator, -1f, 1f);
        return ((float)Mathf.Acos(dot)) * Mathf.Rad2Deg;
    }
    public static float AngleFromVectors(float3 from, float3 to) {
        // sqrt(a) * sqrt(b) = sqrt(a * b) -- valid for real numbers
        float kEpsilonNormalSqrt = 1e-15F;
        float fromSqrMg = from[0]*from[0] + from[1]*from[1] + from[2]*from[2];
        float toSqrMg = to[0]*to[0] + to[1]*to[1] + to[2]*to[2];

        Vector3 fromV3 = new Vector3(from[0],from[1],from[2]);
        Vector3 toV3 = new Vector3(to[0],to[1],to[2]);
        float denominator = (float)Mathf.Sqrt(fromSqrMg * toSqrMg);
        if (denominator < kEpsilonNormalSqrt) return 0f;
        float dot = Mathf.Clamp(Vector3.Dot(fromV3, toV3) / denominator, -1f, 1f);
        return ((float)Mathf.Acos(dot)) * Mathf.Rad2Deg;
    }
}
