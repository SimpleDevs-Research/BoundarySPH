using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class SimulationObjManager : MonoBehaviour
{
    [System.Serializable]
    public class SimulationObstacle {
        public MeshFilter obstacle;
        [Range(0f,1f)]
        public float friction_coefficient = 0f;
        [Range(0f,1f)]
        public float restitution_coefficient = 0.25f;
    }

    public struct Vertex {
        public uint obstacleIndex;
        public float3 position;
        public float3 normal;
        public Vertex(uint obstacleIndex, float3 pos) {
            this.obstacleIndex = obstacleIndex;
            this.position = pos;
            this.normal = new float3(0f,0f,0f);
        }
    }
    public struct Edge {
        public uint obstacleIndex;
        public uint2 vertices;
        public uint2 triangles;
        public float3 midpoint;
        public float3 normal;
        public Edge(uint obstacleIndex, uint2 vertices, uint2 triangles) {
            this.obstacleIndex = obstacleIndex;
            this.vertices = vertices;
            this.triangles = triangles;
            this.midpoint = new float3(0f,0f,0f);
            this.normal = new float3(0f,0f,0f);
        }
    }
    public struct Triangle {
        public uint obstacleIndex;
        public uint3 vertices;
        public float3 angles;
        public uint3 edges;
        public float3 center;
        public float3 normal;
        public float distance_to_center;
        public float3 v1v2n;
        public float3 v2v3n;
        public float3 v1v3n;
        public Triangle(
                uint obstacleIndex,
                uint3 vertices,
                float3 angles,
                uint3 edges,
                float3 center,
                float3 normal,
                float distance_to_center,
                float3 v1v2n,
                float3 v2v3n,
                float3 v1v3n
        ) {
            this.obstacleIndex = obstacleIndex;
            this.vertices = vertices;
            this.angles = angles;
            this.edges = edges;
            this.center = center;
            this.normal = normal;
            this.distance_to_center = distance_to_center;
            this.v1v2n = v1v2n;
            this.v2v3n = v2v3n;
            this.v1v3n = v1v3n;
        }

    }
    
    [System.Serializable]
    public struct Obstacle {
        public uint index;
        public float friction_coefficient;
        public float restitution_coefficient;
        public uint2 vertex_range;      // [0] = starting index in global vertices, [1] = # of vertices associated with this obstacle. _+
        public uint2 triangles_range;   // [0] = starting index in global triangles, [1] = # of vertices associated with this obstacle.
        public uint2 edges_range;       // [0] = starting index in global edges, [1] = # of edges associated with this obstacle
        
        public Obstacle(
                uint index, 
                float friction, float restitution,
                uint2 vrange,
                uint2 erange,
                uint2 trange
        ) {
            this.index = index;
            this.friction_coefficient = friction;
            this.restitution_coefficient = restitution;
            this.vertex_range = vrange;
            this.edges_range = erange;
            this.triangles_range = trange;
        }
    }

    public BufferManager _BM;
    public Simulation3D _SIM;
    public SimulationObstacle[] meshObjects;
    
    //Lists that will be pushed into the GPU
    public List<Obstacle> global_obstacles;
    private List<Vertex> global_vertices;
    private List<Edge> global_edges;
    private List<Triangle> global_triangles;

    void OnDrawGizmos() {
        if (!Application.isPlaying) return;
        Gizmos.color = Color.blue;
        /*
        foreach(Vertex v in global_vertices) {
            int obstacle_id = (int)v.obstacleIndex;
            Debug.Log($"Obstacle ID: {obstacle_id}");
            Vector3 vpos = meshObjects[obstacle_id].obstacle.transform.TransformPoint(v.position);
            Gizmos.DrawSphere(vpos, 0.25f);
            Vector3 normP = meshObjects[obstacle_id].obstacle.transform.TransformPoint(v.normal);
            Gizmos.DrawLine(vpos, vpos+normP);
        }
        foreach(Edge e in global_edges) {
            int obstacle_id = (int)e.obstacleIndex;
            int min_vertex_index = (int)global_obstacles[obstacle_id].vertex_range[0];
            Vector3 v1 = meshObjects[obstacle_id].obstacle.transform.TransformPoint(global_vertices[min_vertex_index + (int)e.vertices[0]].position);
            Vector3 v2 = meshObjects[obstacle_id].obstacle.transform.TransformPoint(global_vertices[min_vertex_index + (int)e.vertices[1]].position);
            Gizmos.DrawLine(v1,v2);
            Vector3 midpoint = meshObjects[obstacle_id].obstacle.transform.TransformPoint(e.midpoint);
            Gizmos.DrawSphere(midpoint, 0.25f);
            Vector3 norm = meshObjects[obstacle_id].obstacle.transform.TransformDirection(e.normal);
            Gizmos.DrawLine(midpoint, midpoint+norm);
        }
        */
        //foreach(Triangle t in global_triangles) {
            Triangle t = global_triangles[0];
            int obstacle_id = (int)t.obstacleIndex;
            int min_vertex_index = (int)global_obstacles[obstacle_id].vertex_range[0];
            int min_edge_index = (int)global_obstacles[obstacle_id].edges_range[0];
            
            Vector3 v1 = meshObjects[obstacle_id].obstacle.transform.TransformPoint(global_vertices[min_vertex_index + (int)t.vertices[0]].position);
            Vector3 v2 = meshObjects[obstacle_id].obstacle.transform.TransformPoint(global_vertices[min_vertex_index + (int)t.vertices[1]].position);
            Vector3 v3 = meshObjects[obstacle_id].obstacle.transform.TransformPoint(global_vertices[min_vertex_index + (int)t.vertices[2]].position);
            Gizmos.color = Color.blue;
            Gizmos.DrawSphere(v1, 0.25f);
            Gizmos.color = Color.red;
            Gizmos.DrawSphere(v2, 0.25f);
            Gizmos.color = Color.black;
            Gizmos.DrawSphere(v3, 0.25f);
            
            Vector3 midpoint = meshObjects[obstacle_id].obstacle.transform.TransformPoint(t.center);
            Gizmos.DrawSphere(midpoint, 0.2f);
            Vector3 norm = meshObjects[obstacle_id].obstacle.transform.TransformDirection(t.normal);
            Gizmos.DrawLine(midpoint, midpoint+norm);

            Edge e1 = global_edges[min_edge_index + (int)t.edges[0]];
            Edge e2 = global_edges[min_edge_index + (int)t.edges[1]];
            Edge e3 = global_edges[min_edge_index + (int)t.edges[2]];
            Vector3 e1m = meshObjects[obstacle_id].obstacle.transform.TransformPoint(e1.midpoint);
            Vector3 e2m = meshObjects[obstacle_id].obstacle.transform.TransformPoint(e2.midpoint);
            Vector3 e3m = meshObjects[obstacle_id].obstacle.transform.TransformPoint(e3.midpoint);
            Vector3 v1v2n = meshObjects[obstacle_id].obstacle.transform.TransformDirection(t.v1v2n);
            Vector3 v2v3n = meshObjects[obstacle_id].obstacle.transform.TransformDirection(t.v2v3n);
            Vector3 v1v3n = meshObjects[obstacle_id].obstacle.transform.TransformDirection(t.v1v3n);
            Gizmos.DrawLine(e1m, e1m+v1v2n);
            Gizmos.DrawLine(e2m, e2m+v1v3n);
            Gizmos.DrawLine(e3m, e3m+v2v3n);
            
            Gizmos.DrawLine(v1, v3);
            Gizmos.DrawLine(v2, v3);
            /*
            Edge e1 = global_edges[min_edge_index + (int)t.edges[0]];
            int evmin = (int)global_obstacles[(int)e1.obstacleIndex].vertex_range[0];
            Vector3 ev1 = meshObjects[(int)e1.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e1.vertices[0]].position);
            Vector3 ev2 = meshObjects[(int)e1.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e1.vertices[1]].position);
            Gizmos.DrawLine(ev1, ev2);
            Edge e2 = global_edges[min_edge_index + (int)t.edges[1]];
            evmin = (int)global_obstacles[(int)e2.obstacleIndex].vertex_range[0];
            ev1 = meshObjects[(int)e2.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e2.vertices[0]].position);
            ev2 = meshObjects[(int)e2.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e2.vertices[1]].position);
            Gizmos.DrawLine(ev1, ev2);
            Edge e3 = global_edges[min_edge_index + (int)t.edges[2]];
            evmin = (int)global_obstacles[(int)e3.obstacleIndex].vertex_range[0];
            ev1 = meshObjects[(int)e3.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e3.vertices[0]].position);
            ev2 = meshObjects[(int)e3.obstacleIndex].obstacle.transform.TransformPoint(global_vertices[evmin+(int)e3.vertices[1]].position);
            Gizmos.DrawLine(ev1, ev2);
            */
        //}
    }

    private void Start() {
        PreprocessObstacles();
    }

    private void PreprocessObstacles() {
        // Initialize the global lists
        global_obstacles = new List<Obstacle>();
        global_vertices = new List<Vertex>();
        global_edges = new List<Edge>();
        global_triangles = new List<Triangle>();

        // Iterate through all obstacles, initialzie their data, save global results to the arrays above
        for(int i = 0; i < meshObjects.Length; i++) {
            PreprocessObj(
                i, meshObjects[i],
                ref global_obstacles, ref global_vertices, ref global_edges, ref global_triangles
            );
        }

        // With all our data properly encapsulated, let's put them into the necessary buffers... assuming we have a buffer manager linked
        if (_BM == null) return;
    }

    private void PreprocessObj(
            int index, SimulationObstacle simObstacle,
            ref List<Obstacle> _obstacles, ref List<Vertex> _vertices, ref List<Edge> _edges, ref List<Triangle> _triangles
    ) {

        // Get mesh data
        Mesh m = simObstacle.obstacle.sharedMesh;
        var vs = m.vertices;
        var ts = m.triangles;

        // Generate the result arrays for the vertices, edges, and triangles.
        // These lists will extend `_obstacles`, `_vertices`, `_edges`, and `_triangles`. So the ultimate goal is to populate these arrays and details.
        uint obstacleIndex = (uint)index;
        Triangle[] triangles = new Triangle[ts.Length/3];
        List<Vertex> vertices = new List<Vertex>();
        List<Edge> edges = new List<Edge>();
        
        // Let's iterate through our triangles. We'll look through all our triangles, looking at our vertices and edges, and iterating the vertices and edges in the meanwhile.
        // SET UIP VERTEX-RELATED MAPPERS
        List<Vector3> condensed_vertex_positions = new List<Vector3>();
        int[] vs_to_condensed_map = new int[vs.Length];
        // SET UP EDGE-RELATED MAPPERS
        Dictionary<uint2, int> es_map = new Dictionary<uint2, int>();
        for(int i = 0; i < ts.Length; i+=3) {
            int ti = i/3;                                           // Get the index for the current triangle
            
            // PROCESS TRIANGLE VERTICES
            int[] tvs = new int[] { ts[i], ts[i+1], ts[i+2] };      // get index of first, second, and third vertices.
            for(int vi = 0; vi < 3; vi++) {                         // Loop through each of the triangle's vertex indices
                Vector3 v = vs[tvs[vi]];
                int mi = condensed_vertex_positions.IndexOf(v);           // Double-check if we have processed this vertex yet.
                if (mi == -1) {                                         // If not, we need to process it.
                    mi = condensed_vertex_positions.Count;            // Use the length of the condensed vertex positions as this vertex's key, then add this vertex
                    condensed_vertex_positions.Add(v);
                    float3 vpos =new float3(v.x, v.y, v.z);             // Create float3 version of vertex position, correct for small float point issues, and add to `vertices` too
                    if (Mathf.Abs(vpos[0]) < 0.000000000001f) vpos[0] = 0f;
                    if (Mathf.Abs(vpos[1]) < 0.000000000001f) vpos[1] = 0f;
                    if (Mathf.Abs(vpos[2]) < 0.000000000001f) vpos[2] = 0f;
                    vertices.Add(new Vertex(obstacleIndex, vpos));
                }
                vs_to_condensed_map[tvs[vi]] = mi;                  // Also add the vertex index (the original one, from the triangle) to `vs_condensed_map`.
            }
            print($"Triangle {ti} # vertices: {vertices.Count}");

            // TRIANGLE PROPERTIES
            // With the vertices, we can now calculate properties like the triangle's center and normal
            int3 tvsm = new int3(vs_to_condensed_map[tvs[0]], vs_to_condensed_map[tvs[1]], vs_to_condensed_map[tvs[2]]);
            print($"Triangle {ti} vertex indices: {tvsm}");
            Vector3 v1p = condensed_vertex_positions[tvsm[0]];       // Calculate the vertex positions of each triangle
            Vector3 v2p = condensed_vertex_positions[tvsm[1]];
            Vector3 v3p = condensed_vertex_positions[tvsm[2]];
            float3 v1f = new float3(v1p.x, v1p.y, v1p.z);
            float3 v2f = new float3(v2p.x, v2p.y, v2p.z);
            float3 v3f = new float3(v3p.x, v3p.y, v3p.z);
            float3 v1fn = Unity.Mathematics.math.normalize(v1f);
            float3 v2fn = Unity.Mathematics.math.normalize(v2f);
            float3 v3fn = Unity.Mathematics.math.normalize(v3f);
            // Calculate the inner angles of each vertex of this triangle
            float3 tangles = new float3(
                AngleFromVectors(v3fn-v1fn,v2fn-v1fn),
                AngleFromVectors(v1fn-v2fn,v3fn-v2fn),
                AngleFromVectors(v2fn-v3fn,v1fn-v3fn)
            );
            print($"Triangle {ti} vertex angles: {tangles}");
            float3 centerF = (v1f+v2f+v3f)/3f;                       // Calculate the center of this triangle
            Vector3 n = Vector3.Cross(v2p-v1p, v3p-v1p).normalized;     // Calculate the normal vector, following the RH rule of Unity
            float3 normalF = new(n.x, n.y, n.z);
            // Update the normal vectors of each vertex, using normalF and tangles
            Vertex v1 = vertices[vs_to_condensed_map[tvs[0]]];
            v1.normal += normalF * tangles[0];
            vertices[vs_to_condensed_map[tvs[0]]] = v1;
            Vertex v2 = vertices[vs_to_condensed_map[tvs[1]]];
            v2.normal += normalF * tangles[1];
            vertices[vs_to_condensed_map[tvs[1]]] = v2;
            Vertex v3 = vertices[vs_to_condensed_map[tvs[2]]];
            v3.normal += normalF * tangles[2];
            vertices[vs_to_condensed_map[tvs[2]]] = v3;

            // PROCESS TRIANGLE EDGES
            // The honestly crummy thing is that int2 can't be used as a key... if it's a list
            // So we used a dictionary `es_map` to temporarily link int2 pairs to edges in `edges`
            // v1v2
            // Step 1: sort the indices to guarantee increasing order 
            int es_map_index_v1v2, es_map_index_v1v3, es_map_index_v2v3;
            uint2 edge1_index = (tvsm[0] < tvsm[1]) ? new uint2((uint)tvsm[0], (uint)tvsm[1]) : new uint2((uint)tvsm[1], (uint)tvsm[0]);      
            if (!es_map.ContainsKey(edge1_index)) {                                             // Step 2: If `edge1_index` is not inside `es_map`, then it's also not inside `edges`.
                es_map_index_v1v2 = edges.Count;                                                // We create new corresponding entries in both `edges` and `es_map`
                Edge e = new Edge(obstacleIndex, edge1_index, new uint2((uint)ti,(uint)ts.Length));   // We pass a uint2 for triangles that is actually only the first. We'll update this once we see this edge again later.
                edges.Add(e);
                es_map.Add(edge1_index, es_map_index_v1v2);                                     // TODO: remember what `es_dynamic.Add(new(0f,0f,0f));` was for... it's legacy code
            } else {
                es_map_index_v1v2 = es_map[edge1_index];                                    // The entry already exists for this edge. this is the 2nd triangle associated with this edge
                Edge e = edges[es_map_index_v1v2]; 
                e.triangles = new uint2(e.triangles[0], (uint)ti);                              // Since this is the second triangle of this edge, we need to update its entry for the 2nd triangle index.
                //Edge e = new Edge(obstacleIndex, edge1_index, new uint2(edges[es_map_index_v1v2].triangles[0], ti));
                edges[es_map_index_v1v2] = e;
            }
            // v1v3
            // Step 1: sort the indices to guarantee increasing order
            uint2 edge2_index = (tvsm[0] < tvsm[2]) ? new uint2((uint)tvsm[0], (uint)tvsm[2]) : new uint2((uint)tvsm[2], (uint)tvsm[0]);
            if (!es_map.ContainsKey(edge2_index)) {                                             // Step 2: If `edge2_index` is not inside `es_map`, then it's also not inside `edges`.
                es_map_index_v1v3 = edges.Count;                                            // We create new corresponding entries in both `edges` and `es_map`
                Edge e = new Edge(obstacleIndex, edge2_index, new uint2((uint)ti,(uint)ts.Length));
                edges.Add(e);                                                                   // TODO: remember what `es_dynamic.Add(new(0f,0f,0f));` was for... it's legacy code
                es_map.Add(edge2_index, es_map_index_v1v3);
            } else {
                es_map_index_v1v3 = es_map[edge2_index];                                    // The entry already exists for this edge. This is the 2nd triangle associated with this edge
                Edge e = edges[es_map_index_v1v3];
                e.triangles = new uint2(e.triangles[0], (uint)ti);
                //Edge e = new Edge(obstacleIndex, edge2_index, new uint2(edges[es_map_index_v1v3].triangles[0], ti));
                edges[es_map_index_v1v3] = e;
            }
            // v2v3
            // Step 1: sort the indices to guarantee increasing order
            uint2 edge3_index = (tvsm[1] < tvsm[2]) ? new uint2((uint)tvsm[1], (uint)tvsm[2]) : new uint2((uint)tvsm[2], (uint)tvsm[1]);
            if (!es_map.ContainsKey(edge3_index)) {                                             // Step 2: If `edge3_index` is not inside `es_map`, then it's also not inside `edges`.
                es_map_index_v2v3 = edges.Count;                                            // We create new corresponding entries in both `edges` and `es_map`
                Edge e = new Edge(obstacleIndex, edge3_index, new((uint)ti,(uint)ts.Length));
                edges.Add(e);                                                                   // TODO: remember what `es_dynamic.Add(new(0f,0f,0f));` was for... it's legacy code       
                es_map.Add(edge3_index, es_map_index_v2v3);
            } else {
                es_map_index_v2v3 = es_map[edge3_index];                                    // The entry already exists for this edge. This is the 2nd triangle associated with this edge
                Edge e = edges[es_map_index_v2v3];
                e.triangles = new uint2(e.triangles[0], (uint)ti);
                //Edge e = new Edge(obstacleIndex, edge3_index, new uint2(edges[es_map_index_v2v3].triangles[0], ti));
                edges[es_map_index_v2v3] = e;
            }

            // CALCULATE 2D NORMALS OF EDGES
            // A pre-processing step we can conduct is the pre-computatin of v1v2n, v2v3n, v1v3n
            float3 v1v2 = v2.position - v1.position;
            float3 v1c = centerF - v1.position;
            float3 v1v2n = -Unity.Mathematics.math.normalize(v1c - (Unity.Mathematics.math.dot(v1c,v1v2)/Unity.Mathematics.math.dot(v1v2,v1v2))*v1v2);
            float3 v1v3 = v3.position - v1.position;
            float3 v1v3n = -Unity.Mathematics.math.normalize(v1c - (Unity.Mathematics.math.dot(v1c,v1v3)/Unity.Mathematics.math.dot(v1v3,v1v3))*v1v3);
            float3 v2v3 = v3.position - v2.position;
            float3 v2c = centerF - v2.position;
            float3 v2v3n = -Unity.Mathematics.math.normalize(v2c - (Unity.Mathematics.math.dot(v2c,v2v3)/Unity.Mathematics.math.dot(v2v3,v2v3))*v2v3);

            // CREATE TRIANGLES
            // We have everything to create our trianges
            triangles[ti] = new Triangle(
                obstacleIndex,                                              // obstacle index
                new uint3((uint)tvsm[0], (uint)tvsm[1], (uint)tvsm[2]),     // vertices
                tangles,                                                    // angles
                new uint3((uint)es_map_index_v1v2,(uint)es_map_index_v1v3,(uint)es_map_index_v2v3),     // edges
                centerF,                                                    // center
                normalF,                                                     // normal vector,
                Unity.Mathematics.math.dot(-1f * centerF, normalF),       // signed distance from object origin to center
                v1v2n, 
                v2v3n,                
                v1v3n
            );
        }

        // FURTHER CORRECT EDGES
        // The vertices and triangles are fine. Now we just need to calculate the normal vectors for each edge
        for(int i = 0; i < edges.Count; i++) {
            uint2 evs = edges[i].vertices;
            float3 midpoint = (vertices[(int)evs[0]].position + vertices[(int)evs[1]].position)/2f;
            // We calculate the normal by approximating the normals of each triangle assoicated with this edge. We have to be careful about edges with just one triangle (open edges).
            uint2 ets = edges[i].triangles;
            float3 norm = triangles[ets[0]].normal;
            if (ets[1] != (uint)ts.Length) norm += (triangles[ets[1]].normal);
            // We can now update the midpoint and localNormalEnd of this edge
            //OP.EdgeStatic e_static = new OP.EdgeStatic();
            Edge e = edges[i];
            e.midpoint = midpoint;
            e.normal = Unity.Mathematics.math.normalize(norm);
            edges[i] = e;
            print($"Edge {i}: {evs}");
        }

        // UPDATE OBSTACLES
        // Before we do anything else, we need to update our obstacle. We do this now because before we add our vertices, edges, and triangles to the global lists
        // we want to use the length of the global lists as the starting index of each of this object's ranges
        Obstacle obs = new Obstacle(
            obstacleIndex, 
            simObstacle.friction_coefficient, 
            simObstacle.restitution_coefficient,
            new uint2((uint)_vertices.Count, (uint)vertices.Count),        // vertex range
            new uint2((uint)_edges.Count, (uint)edges.Count),              // edge range
            new uint2((uint)_triangles.Count, (uint)triangles.Length)      // triangle range
        );
        
        // Finally, we update our output lists
        _vertices.AddRange(vertices);
        _edges.AddRange(edges);
        _triangles.AddRange(triangles);
        _obstacles.Add(obs);
    }

    public void UpdateBuffers() {

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
        return Mathf.Abs(((float)Mathf.Acos(dot)) * Mathf.Rad2Deg);
    }
}
