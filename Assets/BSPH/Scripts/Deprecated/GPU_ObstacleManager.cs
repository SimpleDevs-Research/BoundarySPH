using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using Random = UnityEngine.Random;
using System.Linq;
using OP = ObstaclePrimitives.Structs;

public class GPU_ObstacleManager : MonoBehaviour
{
    public BufferManager _BM;

    public enum ShowHideEnum { Show, Hide }
    public ShowHideEnum showEditorControls = ShowHideEnum.Hide;

    public enum UpdateOrCoroutine { Update, Coroutine }
    public UpdateOrCoroutine updateMethod = UpdateOrCoroutine.Coroutine;

    [System.Serializable]
    public struct ObstacleBounds {
        public float3 lowerBounds;
        public float3 upperBounds;
        public int isActive;
        public float frictionCoefficient;
    }

    [System.Serializable]
    public struct ObstacleTriangle {
        public int obstacleIndex;
        public int3 vertexIndices;
        public float3 c;
        public float3 n;
        public float d;
        public float3 v1v2, v1v3, v2v3;
        public float3 v1v2n, v1v3n, v2v3n;
        public float3 lowerBound, upperBound;
    }

    [System.Serializable]
    public struct ObstacleVertex {
        public int obstacleIndex;
        public float3 position;
        public float3 normal;
    }
    [System.Serializable]
    public class VertexNormal {
        public List<int> vertexIndices = new List<int>();
        public List<Vector3> norms = new List<Vector3>();
        public VertexNormal() {
            vertexIndices = new List<int>();
            norms = new List<Vector3>();
        }
        public VertexNormal(int index, Vector3 norm) {
            vertexIndices = new List<int>();
            norms = new List<Vector3>();
            vertexIndices.Add(index);
            norms.Add(norm);
        }
        public void Add(int index, Vector3 norm) {
            if (!vertexIndices.Contains(index)) vertexIndices.Add(index);
            //if (!norms.Contains(norm)) norms.Add(norm);
            norms.Add(norm);
        }
        public Vector3 GetAverageNorm(bool normalize = true) {
            Vector3 sum = Vector3.zero;
            foreach(Vector3 n in norms) sum += n;
            //Vector3 avg = new Vector3(sum.x / norms.Count, sum.y / norms.Count, sum.z / norms.Count);
            if (normalize) return sum.normalized;
            return sum;
        }
    }

    [System.Serializable] 
    public class Obstacle {
        public Transform obstacle;
        public bool isDynamic;
        [Range(0f,1f)] public float frictionCoefficient = 0f;
        [HideInInspector] public int2 triangleIndices;
        [HideInInspector] public int2 vertexIndices;
        [HideInInspector] public ObstacleVertex[] vertices;
        [HideInInspector] public List<ObstacleTriangle> triangles;
        [HideInInspector] public ObstacleBounds bounds;
        [HideInInspector] public int index;

        [HideInInspector] public float prevFriction;
    }

    [Header("== REFERENCES ==")]
    [Tooltip("Reference to a ParticleGrid component that acts as this controller's grid system")]
    public Grid _GRID = null;
    public ComputeShader _SHADER = null;
    public ParticleController _PARTICLE_CONTROLLER = null;

    [SerializeField] private List<Obstacle> _OBSTACLES = new List<Obstacle>();
    private int _numVertices, _numTriangles;
    [SerializeField, ReadOnlyInsp] private bool _initialized = false;
    [SerializeField] private bool _toggle_gizmos = true;
    [SerializeField] private bool _toggle_debug = false;

    [SerializeField] private OP.Projection[] debug_projections;
    [SerializeField] private Color gizmos_bounds_color = Color.yellow;
    private bool show_bounds => gizmos_bounds_color.a > 0f;
    [SerializeField] private Color gizmos_vertices_color = Color.white;
    private bool show_vertex_details => gizmos_vertices_color.a > 0f;
    [SerializeField] private Color gizmos_triangles_color = Color.red;
    private bool show_triangle_details => gizmos_triangles_color.a > 0f;
    [SerializeField] private Color gizmos_projection_positions_color = Color.green;
    private bool show_projection_positions => gizmos_projection_positions_color.a > 0f;
    [SerializeField] private Color gizmos_projection_normals_color = Color.green;
    private bool show_projection_normals => gizmos_projection_normals_color.a > 0f;
    private bool show_projections => show_projection_positions || show_projection_normals;
    private bool show_gizmos => _toggle_gizmos && (show_bounds || show_vertex_details || show_triangle_details || show_projections);


    private void OnDrawGizmos() {
        if (!Application.isPlaying) return;
        if (!_initialized) return;
        if (_OBSTACLES.Count == 0) return;        
        if (!show_gizmos) return;

        if (show_bounds) {
            ObstacleBounds[] bounds = new ObstacleBounds[_OBSTACLES.Count];
            _BOUNDS_BUFFER.GetData(bounds);
            Gizmos.color = gizmos_bounds_color;
            foreach(ObstacleBounds b in bounds) {
                float3 center = (b.upperBounds + b.lowerBounds)/2f;
                Gizmos.DrawWireCube(center, b.upperBounds - b.lowerBounds);
            }
        }

        if (show_vertex_details) {
            ObstacleVertex[] vertices = new ObstacleVertex[_numVertices];
            _VERTICES_BUFFER.GetData(vertices);
            Gizmos.color = gizmos_vertices_color;
            foreach(ObstacleVertex v in vertices) {
                Gizmos.DrawSphere(v.position,0.01f);
                Gizmos.DrawRay(v.position,v.normal);
            }
        }

        if (show_triangle_details) {
             ObstacleVertex[] vertices = new ObstacleVertex[_numVertices];
            _VERTICES_BUFFER.GetData(vertices);
            ObstacleTriangle[] triangles = new ObstacleTriangle[_numTriangles];
            _TRIANGLES_BUFFER.GetData(triangles);
            foreach(ObstacleTriangle triangle in triangles) {
                Gizmos.color = gizmos_triangles_color;
                Gizmos.DrawSphere(triangle.c,0.02f);
                float3 v1v2mid = (vertices[triangle.vertexIndices[0]].position + vertices[triangle.vertexIndices[1]].position) / 2f;
                Gizmos.DrawRay(v1v2mid, triangle.v1v2n * 2f);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[0]].position, triangle.v1v2n);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[1]].position, triangle.v1v2n);
                float3 v1v3mid = (vertices[triangle.vertexIndices[0]].position + vertices[triangle.vertexIndices[2]].position) / 2f;
                Gizmos.DrawRay(v1v3mid, triangle.v1v3n * 2f);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[0]].position, triangle.v1v3n);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[2]].position, triangle.v1v3n);
                float3 v2v3mid = (vertices[triangle.vertexIndices[1]].position + vertices[triangle.vertexIndices[2]].position) / 2f;
                Gizmos.DrawRay(v2v3mid, triangle.v2v3n * 2f);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[1]].position, triangle.v2v3n);
                //Gizmos.DrawRay(vertices[triangle.vertexIndices[2]].position, triangle.v2v3n);
                Gizmos.color = Color.blue;
                float3 v1v2_3Dn = vertices[triangle.vertexIndices[0]].normal + vertices[triangle.vertexIndices[1]].normal;
                Gizmos.DrawRay(v1v2mid, v1v2_3Dn * 2f);
                float3 v1v3_3Dn = vertices[triangle.vertexIndices[0]].normal + vertices[triangle.vertexIndices[2]].normal;
                Gizmos.DrawRay(v1v3mid, v1v3_3Dn * 2f);
                float3 v2v3_3Dn = vertices[triangle.vertexIndices[1]].normal + vertices[triangle.vertexIndices[2]].normal;
                Gizmos.DrawRay(v2v3mid, v2v3_3Dn * 2f);
                Gizmos.DrawRay(vertices[triangle.vertexIndices[0]].position, vertices[triangle.vertexIndices[0]].normal * 2f);
                Gizmos.DrawRay(vertices[triangle.vertexIndices[1]].position, vertices[triangle.vertexIndices[1]].normal * 2f);
                Gizmos.DrawRay(vertices[triangle.vertexIndices[2]].position, vertices[triangle.vertexIndices[2]].normal * 2f);

                Gizmos.DrawRay(triangle.c, triangle.n * 2f);
            }
        }

        if (_PARTICLE_CONTROLLER != null && show_projections) {
            OP.Particle[] particles = new OP.Particle[_PARTICLE_CONTROLLER.numParticles];
            _BM.PARTICLES_BUFFER.GetData(particles);
            debug_projections = new OP.Projection[_PARTICLE_CONTROLLER.numParticles];
            _BM.PARTICLES_EXTERNAL_FORCES_BUFFER.GetData(debug_projections);
            debug_projections = debug_projections.Where(p=>p.counter > 0).ToArray();

            OP.Projection proj;
            OP.Particle part;
            for(int i = 0; i < debug_projections.Length; i++) {
                part = particles[i];
                proj = debug_projections[i];
                if (proj.counter <= 0) continue;
                if (show_projection_positions) {
                    Gizmos.color = gizmos_projection_positions_color;
                    Gizmos.DrawSphere(proj.position,0.1f);
                }
                if (show_projection_normals) {
                    Vector3 n = new Vector3((float)proj.normal[0] / 1024f, (float)proj.normal[1] / 1024f, (float)proj.normal[2] / 1024f);
                    Gizmos.color = gizmos_projection_normals_color;
                    Gizmos.DrawRay(proj.position, n);
                }
            }
        }
    }

    private void Awake() {
        // Since this is just an onEventListener for when the simulation starts, we just call `Initialize()`.
        Initialize();
    }

    // This function can be called from either the inspector via the "Manual Initialization" button, or when the scene first runs
    // This function performs some key functionalities, per se:
    // 1. Updates all obstacles such that their entries in the _OBSTACLES list are filled with the necessary information for the shader
    public void Initialize() {
        // There are two checks we need to make:
        // 1. We double-check to make sure that the shader reference is set.
        // 2. We double-check the that _OBSTACLES list is actually populated with at least one transform element
        // Without thise two checks passing, we cannot initiate the coroutine for the update loop.
        if (_SHADER == null) {
            Debug.LogError("ERROR: No Compute Shader referenced. Please add a reference to a compatible Compute Shader");
            return;
        }
        if (_OBSTACLES.Count == 0) {
            Debug.LogError("ERROR: No obstacle transforms have been added in the inspector.");
            return;
        }

        // Assuming we get past the first two checks, we then proceed with the next: calling `UpdateObstacles` to check that we have the necessary obstacles
        // Keep in mind that this one will also remove any obstacles that do not have any MeshRenderer and MeshFilter.
        // We'll perform a 2nd check to make sure that the updated _OBSTACLES list still has any remaining obstacles.
        UpdateObstacles();
        if (_OBSTACLES.Count == 0) {
            Debug.LogError("ERROR: After filtering, no obstacles fit the necessary criteria. Please add obstacles that meet the criteria.");
            return;
        }

        // At this point, it's good to start initializing our shader.
        // Firstly, we initialize kernels. We CANNOT do initialization of variables yet - we'll do that in `InitializeBuffers()`.
        InitializeKernels();

        // Secondly, we initialize the buffers. This function is humongous because it haandles a lot at once. Namely, it:
        // 1. Initializes the ComputeBuffers needed for the simulation
        // 2. Combines the data from _OBSTACLES into a format that fits for the buffers
        // 3. populates the buffers with the data from step #2
        InitializeShader();

        // Finally, we declare that we've been initialized. 
        // This makes it easier for when, during the update loop, whenever an obstacle gets updated (ex. moved, scaled, rotated), we can let the system know to replace the necessary values in their respective buffers
        _initialized = true;
    }

    // Always note that the # of vertices and triangles WILL NOT CHANGE even when updating the obstacle here.
    // This means that this function will always be able to update _VERTICES, _TRIANGLES, and _BOUNDS 
    public void UpdateObstacles() {
        // Firstly, filter out any obstacles that don't fit the requirements.
        _OBSTACLES = _OBSTACLES.Where(obs => 
            obs.obstacle != null 
            && ((obs.obstacle.GetComponent<MeshFilter>() != null && obs.obstacle.GetComponent<MeshRenderer>() != null)||obs.obstacle.GetComponent<SkinnedMeshRenderer>() != null)
        ).ToList();
        // If the resulting list is empty, we end early
        if (_OBSTACLES.Count == 0) return;
        // For each obstacle, we call `UpdateObstacle` on them.
        for(int i = 0; i < _OBSTACLES.Count; i++) {
            UpdateObstacle(_OBSTACLES[i], i);
        }
    }

    public void PreprocessObstacle(Obstacle obs, int obsIndex) {
        // Grab the mesh data
        Mesh mesh = (obs.obstacle.GetComponent<MeshFilter>() != null) 
            ? obs.obstacle.GetComponent<MeshFilter>().sharedMesh
            : obs.obstacle.GetComponent<SkinnedMeshRenderer>().sharedMesh;

        // Extract the mesh vertices and triangles
        var vs = mesh.vertices;
        var ts = mesh.triangles;

        // Prepare CommonVertices list, which stores common vertices
        List<Vector3> commonVertices = new List<Vector3>();
        List<List<int>> neighborVertices = new List<List<int>>();
        List<List<int>> connectedTriangles = new List<List<int>>();

        int i, neighborI;
        int commonI, neighborCommonI;
        Vector3 localPosition, neighborPosition;

        // Iterate through each triangle:
        for(int ti = 0; ti < ts.Length; ti+=3) {
            for(int vi = 0; vi < 3; vi++) {
                // Combine ti and vi to get i, the current vertex index in ts
                i = ts[ti + vi];
                // Get the local position of the vertex associated with index `i`
                localPosition = vs[i];
                // Check if this local position is located inside of our list of commonVertices
                commonI = commonVertices.IndexOf(localPosition);
                if (commonI == -1) {
                    // Add the current vertex!
                    commonI = commonVertices.Count;
                    commonVertices.Add(localPosition);
                    neighborVertices.Add(new List<int>());
                    connectedTriangles.Add(new List<int>());
                }
                // Add the current triangle if it doesn't exist already
                if (!connectedTriangles[commonI].Contains(ti)) connectedTriangles[commonI].Add(ti);
                // We need to identify which are neighbor vertices. We know based on current triangle.
                for(int vi2 = 0; vi2 < 3; vi2++) {
                    // We ignore ourselves
                    if (vi2 == vi) continue;
                    // Grab neighbor's index in vs
                    neighborI = ts[ti + vi2];
                    // Grab the neighbor's local position
                    neighborPosition = vs[neighborI];
                    // Check if this local position is inside commonVertices
                    neighborCommonI = commonVertices.IndexOf(neighborPosition);
                    // if that neighbor doesn't exist, we add it
                    if (neighborCommonI == -1) {
                        neighborCommonI = commonVertices.Count;
                        commonVertices.Add(neighborPosition);
                        neighborVertices.Add(new List<int>());
                        connectedTriangles.Add(new List<int>());
                    }
                    if (!connectedTriangles[neighborCommonI].Contains(ti)) connectedTriangles[commonI].Add(ti);
                    // We can now add the neighbor to our neighbors list

                }
            }
            // Extract local space position of 
        }
    }

    public void UpdateObstacle(Obstacle obs, int obsIndex) {
        Mesh mesh = (obs.obstacle.GetComponent<MeshFilter>() != null) 
            ? obs.obstacle.GetComponent<MeshFilter>().sharedMesh
            : obs.obstacle.GetComponent<SkinnedMeshRenderer>().sharedMesh;

        var vs = mesh.vertices;
        var ts = mesh.triangles;
        /*
        Vector4 rot = new Vector4(
            obs.obstacle.rotation.x, 
            obs.obstacle.rotation.y, 
            obs.obstacle.rotation.z, 
            obs.obstacle.rotation.w
        );
        */

        // Initialize some variables
        ObstacleTriangle triangle;
        float[] valX,valY,valZ;
        Vector3 v1, v2, v3;
        ObstacleVertex v1f, v2f, v3f;

        ObstacleVertex[] vertices = new ObstacleVertex[vs.Length];
        Dictionary<Vector3,VertexNormal> vertexNormalDict = new Dictionary<Vector3,VertexNormal>();
        List<ObstacleTriangle> triangles = new List<ObstacleTriangle>();

        // Loop through triangles. We iterate through `ts.Length` this time.
        for(int t = 0; t < ts.Length; t+=3) {
            // Create a new triangle
            triangle = new ObstacleTriangle();
            // Grab the vertices and convert them into world scale
            v1 = obs.obstacle.TransformPoint(vs[ts[t]]);
            v2 = obs.obstacle.TransformPoint(vs[ts[t+1]]);
            v3 = obs.obstacle.TransformPoint(vs[ts[t+2]]);
            // Get the bounds of the triangle
            valX = new [] {v1.x,v2.x,v3.x};
            valY = new [] {v1.y,v2.y,v3.y};
            valZ = new [] {v1.z,v2.z,v3.z};
            triangle.lowerBound = new(
                valX.Min() - _PARTICLE_CONTROLLER.particleRenderRadius, 
                valY.Min() - _PARTICLE_CONTROLLER.particleRenderRadius, 
                valZ.Min() - _PARTICLE_CONTROLLER.particleRenderRadius
            );
            triangle.upperBound = new(
                valX.Max() + _PARTICLE_CONTROLLER.particleRenderRadius, 
                valY.Max() + _PARTICLE_CONTROLLER.particleRenderRadius, 
                valZ.Max() + _PARTICLE_CONTROLLER.particleRenderRadius
            );
            // Form float3 versions of v1, v2, and v3
            v1f = new ObstacleVertex();
            v1f.obstacleIndex = obsIndex;
            v1f.position = new(v1.x, v1.y, v1.z);
            v2f = new ObstacleVertex();
            v2f.obstacleIndex = obsIndex;
            v2f.position = new(v2.x, v2.y, v2.z);
            v3f = new ObstacleVertex();
            v3f.obstacleIndex = obsIndex;
            v3f.position = new (v3.x, v3.y, v3.z);
            // Add to `_vertices`
            vertices[ts[t]] = v1f;
            vertices[ts[t+1]] = v2f;
            vertices[ts[t+2]] = v3f;
            // Add to `vertexIndices` of current triangle
            triangle.vertexIndices = new(ts[t],ts[t+1],ts[t+2]);
            // Calculate centroid based on average of v1,v2,v3
            Vector3 c = (v1 + v2 + v3) / 3f;
            //triangle.c = (v1f.position + v2f.position + v3f.position) / 3f;
            triangle.c = new(c.x, c.y, c.z);
            // Calculate normal based on v1,v2,v3
            Vector3 norm = Vector3.Cross(v2 - v1, v3 - v1);
            Vector3 normDir = norm.normalized;
            triangle.n = new(normDir.x, normDir.y, normDir.z);
            triangle.d = Vector3.Dot(-1f*c, normDir);
            // Calculate normals of each edge
            Vector3 v1v2 = v2 - v1;
            Vector3 v1c = c - v1;
            Vector3 v1v2n = -Vector3.Normalize(v1c - (Vector3.Dot(v1c,v1v2)/Vector3.Dot(v1v2,v1v2))*v1v2);
            Vector3 v2v3 = v3 - v2;
            Vector3 v2c = c - v2;
            Vector3 v2v3n = -Vector3.Normalize(v2c - (Vector3.Dot(v2c,v2v3)/Vector3.Dot(v2v3,v2v3))*v2v3);
            Vector3 v1v3 = v3-v1;
            Vector3 v1v3n = -Vector3.Normalize(v1c - (Vector3.Dot(v1c,v1v3)/Vector3.Dot(v1v3,v1v3))*v1v3);
            triangle.v1v2 = new(v1v2.x, v1v2.y, v1v2.z);
            triangle.v1v3 = new(v1v3.x, v1v3.y, v1v3.z);
            triangle.v2v3 = new(v2v3.x, v2v3.y, v2v3.z);
            triangle.v1v2n = new(v1v2n.x, v1v2n.y, v1v2n.z);
            triangle.v1v3n = new(v1v3n.x, v1v3n.y, v1v3n.z);
            triangle.v2v3n = new(v2v3n.x, v2v3n.y, v2v3n.z);
            // Calculate normals for vertices, based on edge directions
            // IMPORTANT: the commented-out implementation works best with 2D MESHES (ex. simple quads)...
            /*
            // First, for v1
            if (!vertexNormalDict.ContainsKey(vs[ts[t]]))   vertexNormalDict.Add( vs[ts[t]], new VertexNormal() );
            vertexNormalDict[vs[ts[t]]].Add( ts[t], -v1v2);
            vertexNormalDict[vs[ts[t]]].Add( ts[t], -v1v3);
            // Second, for v2
            if (!vertexNormalDict.ContainsKey(vs[ts[t+1]])) vertexNormalDict.Add( vs[ts[t+1]], new VertexNormal() );
            vertexNormalDict[vs[ts[t+1]]].Add( ts[t+1], v1v2 );
            vertexNormalDict[vs[ts[t+1]]].Add( ts[t+1], -v2v3 );
            // Lastly, for v3
            if (!vertexNormalDict.ContainsKey(vs[ts[t+2]])) vertexNormalDict.Add( vs[ts[t+2]], new VertexNormal() );
            vertexNormalDict[vs[ts[t+2]]].Add( ts[t+2], v2v3 );
            vertexNormalDict[vs[ts[t+2]]].Add( ts[t+2], v1v3 );
            */
            // This implementation version works best with 3D models!
            if (!vertexNormalDict.ContainsKey(vs[ts[t]]))   vertexNormalDict.Add( vs[ts[t]], new VertexNormal(ts[t], normDir) );
            else                                            vertexNormalDict[vs[ts[t]]].Add( ts[t], normDir );
            if (!vertexNormalDict.ContainsKey(vs[ts[t+1]])) vertexNormalDict.Add( vs[ts[t+1]], new VertexNormal(ts[t+1], normDir) );
            else                                            vertexNormalDict[vs[ts[t+1]]].Add( ts[t+1], normDir );
            if (!vertexNormalDict.ContainsKey(vs[ts[t+2]])) vertexNormalDict.Add( vs[ts[t+2]], new VertexNormal( ts[t+2], normDir ) );
            else                                            vertexNormalDict[vs[ts[t+2]]].Add( ts[t+2], normDir );

            // Add triangle to list of triangles we have
            triangle.obstacleIndex = obsIndex;
            triangles.Add(triangle);
        }

        foreach(KeyValuePair<Vector3,VertexNormal> kvp in vertexNormalDict) {
            Vector3 norm = kvp.Value.GetAverageNorm();
            Debug.Log(norm.ToString() + " : " + kvp.Value.vertexIndices.Count.ToString());
            //for(int i = 0; i < kvp.Value.vertexIndices.Count; i++)// in kvp.Value.vertexIndices) {
            foreach(int i in kvp.Value.vertexIndices)
            {
                vertices[i].normal = new(norm.x, norm.y, norm.z);
                //vertices[kvp.Value.vertexIndices[i]].normal = new(kvp.Value.norms[i].x, kvp.Value.norms[i].y, kvp.Value.norms[i].z);
            }
        }

        // We can get world-scale bounds via renderer.boudns
        var bounds = (obs.obstacle.GetComponent<MeshRenderer>() != null)
            ? obs.obstacle.GetComponent<MeshRenderer>().bounds
            : obs.obstacle.GetComponent<SkinnedMeshRenderer>().bounds;
        // To get 6-float bounds, we need to get bounds.min and bounds.max, which are Vector3s
        ObstacleBounds _bounds = new ObstacleBounds();
        // We make sure to limit the mins and maxes to within the simulation space
        
        _bounds.lowerBounds = new(
            Mathf.Max(_GRID.outerBounds[0], bounds.min.x - _GRID.gridCellSize), 
            Mathf.Max(_GRID.outerBounds[1], bounds.min.y - _GRID.gridCellSize), 
            Mathf.Max(_GRID.outerBounds[2], bounds.min.z - _GRID.gridCellSize)
        );
        _bounds.upperBounds = new(
            Mathf.Min(bounds.max.x + _GRID.gridCellSize, _GRID.outerBounds[3]), 
            Mathf.Min(bounds.max.y + _GRID.gridCellSize, _GRID.outerBounds[4]), 
            Mathf.Min(bounds.max.z + _GRID.gridCellSize, _GRID.outerBounds[5])
        );
        
        //_bounds.lowerBounds = new(bounds.min.x,bounds.min.y,bounds.min.z);
        //_bounds.upperBounds = new(bounds.max.x,bounds.max.y,bounds.max.z);
        // We check if the bounds are actually inside our simulation space
        // We can do so by checking if any of the 8 corners of the bounds are within the simulation space
        Vector3[] boundPositions = new Vector3[8] {
            bounds.min - new Vector3(_GRID.gridCellSize, _GRID.gridCellSize, _GRID.gridCellSize),
            new Vector3(bounds.min.x - _GRID.gridCellSize, bounds.min.y - _GRID.gridCellSize, bounds.max.z + _GRID.gridCellSize),
            new Vector3(bounds.min.x - _GRID.gridCellSize, bounds.max.y + _GRID.gridCellSize, bounds.min.z - _GRID.gridCellSize),
            new Vector3(bounds.max.x + _GRID.gridCellSize, bounds.min.y - _GRID.gridCellSize, bounds.min.z - _GRID.gridCellSize),
            new Vector3(bounds.min.x - _GRID.gridCellSize, bounds.max.y + _GRID.gridCellSize, bounds.max.z + _GRID.gridCellSize),
            new Vector3(bounds.max.x + _GRID.gridCellSize, bounds.min.y - _GRID.gridCellSize, bounds.max.z + _GRID.gridCellSize),
            new Vector3(bounds.max.x + _GRID.gridCellSize, bounds.max.y + _GRID.gridCellSize, bounds.min.z - _GRID.gridCellSize),
            bounds.max + new Vector3(_GRID.gridCellSize, _GRID.gridCellSize, _GRID.gridCellSize)
        };
        Vector3 center = (bounds.min + bounds.max)/2f;
        bool inSimulationSpace = (center.x > _GRID.outerBounds[0] && center.x < _GRID.outerBounds[3] && center.y > _GRID.outerBounds[1] && center.y < _GRID.outerBounds[4] && center.z > _GRID.outerBounds[2] && center.z < _GRID.outerBounds[5]);
        foreach(Vector3 bp in boundPositions) {
            if (inSimulationSpace) break;
            inSimulationSpace = 
                inSimulationSpace || 
                (bp.x > _GRID.outerBounds[0] && bp.x < _GRID.outerBounds[3] && bp.y > _GRID.outerBounds[1] && bp.y < _GRID.outerBounds[4] && bp.z > _GRID.outerBounds[2] && bp.z < _GRID.outerBounds[5]);
        }
        _bounds.isActive = (obs.obstacle.gameObject.activeInHierarchy && inSimulationSpace) ? 1 : 0;
        _bounds.frictionCoefficient = obs.frictionCoefficient;

        // Update obstacle
        obs.vertexIndices[1] = vs.Length;
        obs.triangleIndices[1] = triangles.Count;
        obs.vertices = vertices;
        obs.triangles = triangles;
        obs.bounds = _bounds;
        obs.index = obsIndex;
        obs.prevFriction = obs.frictionCoefficient;
        obs.obstacle.hasChanged = false;

        // If the obstacle list is already initialized, then we update that too
        if (_initialized) UpdateObstacleInCombined(obs);
    }

    private int _CLEAR_PROJECTIONS_KERNEL;
    private int _CALCULATE_PROJECTIONS_KERNEL;
    private void InitializeKernels() {
        _CLEAR_PROJECTIONS_KERNEL = _SHADER.FindKernel("ClearProjections");
        _CALCULATE_PROJECTIONS_KERNEL = _SHADER.FindKernel("CalculateProjections");
    }

    private ComputeBuffer _VERTICES_BUFFER;
    private ComputeBuffer _TRIANGLES_BUFFER;
    private ComputeBuffer _BOUNDS_BUFFER;
    private ComputeBuffer _VERTEX_INDEX_LIMITS_BUFFER;
    private ComputeBuffer _TRIANGLE_INDEX_LIMITS_BUFFER;
    private ComputeBuffer _DEBUG_BUFFER;
    // This function only initializes the buffers. It doesn't populate them just yet (that's done later after we initialize them)
    private void InitializeShader() {

        // == STEP #1: Combine data from _OBSTACLES into a usable form
        List<ObstacleVertex> VERTICES = new List<ObstacleVertex>();
        List<ObstacleTriangle> TRIANGLES = new List<ObstacleTriangle>();
        List<ObstacleBounds> BOUNDS = new List<ObstacleBounds>();
        List<int2> VERTEX_INDEX_LIMITS = new List<int2>();
        List<int2> TRIANGLE_INDEX_LIMITS = new List<int2>();
        _numVertices = 0;
        _numTriangles = 0;
        
        foreach(Obstacle obs in _OBSTACLES) {
            // We'll update the current obstacle with the index limits for their portion of _VERTICES and _TRIANGLES
            obs.vertexIndices[0] = _numVertices;
            obs.triangleIndices[0] = _numTriangles;
            // We'll be adding the vertices | triangles | bounds to _VERTICES | _TRIANGLES | _BOUNDS, respectively
            VERTICES.AddRange(obs.vertices);
            TRIANGLES.AddRange(obs.triangles);
            VERTEX_INDEX_LIMITS.Add(obs.vertexIndices);
            TRIANGLE_INDEX_LIMITS.Add(obs.triangleIndices);
            BOUNDS.Add(obs.bounds);
            // Update _numVertices and _numTriangles. No need to update any bounds number since that's tied specifically to index # of obstacle itself
            _numVertices += obs.vertexIndices[1];
            _numTriangles += obs.triangleIndices[1];
        }

        // == STEP #2: Initialize buffers
        if (_VERTICES_BUFFER != null) _VERTICES_BUFFER.Release();
        _VERTICES_BUFFER = new ComputeBuffer(VERTICES.Count, sizeof(int) + sizeof(float)*6);
        if (_TRIANGLES_BUFFER != null) _TRIANGLES_BUFFER.Release();
        _TRIANGLES_BUFFER = new ComputeBuffer(TRIANGLES.Count, sizeof(int)*4 + sizeof(float)*31);
        if (_BOUNDS_BUFFER != null) _BOUNDS_BUFFER.Release();
        _BOUNDS_BUFFER = new ComputeBuffer(BOUNDS.Count, sizeof(float)*7 + sizeof(int));
        if (_VERTEX_INDEX_LIMITS_BUFFER != null) _VERTEX_INDEX_LIMITS_BUFFER.Release();
        _VERTEX_INDEX_LIMITS_BUFFER = new ComputeBuffer(VERTEX_INDEX_LIMITS.Count, sizeof(int)*2);
        if (_TRIANGLE_INDEX_LIMITS_BUFFER != null) _TRIANGLE_INDEX_LIMITS_BUFFER.Release();
        _TRIANGLE_INDEX_LIMITS_BUFFER = new ComputeBuffer(TRIANGLE_INDEX_LIMITS.Count, sizeof(int)*2);
        if (_DEBUG_BUFFER != null) _DEBUG_BUFFER.Release();
        _DEBUG_BUFFER = new ComputeBuffer(_PARTICLE_CONTROLLER.numParticles, sizeof(int));

        // Set buffers inside relevant kernels
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_BOUNDS", _BOUNDS_BUFFER);
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_TRIANGLES", _TRIANGLES_BUFFER);
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_VERTICES", _VERTICES_BUFFER);
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_VERTEX_INDEX_LIMITS", _VERTEX_INDEX_LIMITS_BUFFER);
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_TRIANGLE_INDEX_LIMITS", _TRIANGLE_INDEX_LIMITS_BUFFER);
        _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_DEBUG_BUFFER", _DEBUG_BUFFER);

        // == STEP #3: Populate our buffers with the existing info, and initialize our shader variables
        _VERTICES_BUFFER.SetData(VERTICES.ToArray());
        _TRIANGLES_BUFFER.SetData(TRIANGLES.ToArray());
        _BOUNDS_BUFFER.SetData(BOUNDS.ToArray());
        _VERTEX_INDEX_LIMITS_BUFFER.SetData(VERTEX_INDEX_LIMITS.ToArray());
        _TRIANGLE_INDEX_LIMITS_BUFFER.SetData(TRIANGLE_INDEX_LIMITS.ToArray());
        
        _SHADER.SetInt("_numObstacles", _OBSTACLES.Count);
        _SHADER.SetInt("_numVertices", _numVertices);
        _SHADER.SetInt("_numTriangles", _numTriangles);
        _SHADER.SetInt("_numGridCells", _GRID.numGridCells);
        _SHADER.SetInts("_numCellsPerAxis", _GRID.numCellsPerAxis);
        _SHADER.SetFloat("_gridCellSize", _GRID.gridCellSize);
        _SHADER.SetFloat("_gridScalingX", _GRID.gridScaling[0]);
        _SHADER.SetFloat("_gridScalingY", _GRID.gridScaling[1]);
        _SHADER.SetFloat("_gridScalingZ", _GRID.gridScaling[2]);
    }










    private ComputeBuffer _NUM_CALLS_BUFFER;
    // This is called on the first frame of the simulation.
    // Check it: if the shader wasn't initialized (because maybe the Awake code wasn't completed due to an error), then we prevent the update coroutine from playing out.
    private void Start() {
        if (!_initialized) return;
        
        // If the _PARTICLE_CONTROLLER reference is not null, then we have a special task to do: 
        // 1. Getting its PARTICLE_BUFFER and EXTERNAL_FORCES_BUFFER into our own shader
        if (_PARTICLE_CONTROLLER != null) {
            _SHADER.SetInt("_numParticles", _PARTICLE_CONTROLLER.numParticles);
            _SHADER.SetFloat("_particleRenderRadius",_PARTICLE_CONTROLLER.particleRenderRadius);
            _SHADER.SetBuffer(_CLEAR_PROJECTIONS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
            _SHADER.SetBuffer(_CLEAR_PROJECTIONS_KERNEL, "_PROJECTIONS", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);
            _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
            _SHADER.SetBuffer(_CALCULATE_PROJECTIONS_KERNEL, "_PROJECTIONS", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);
        }

        if (updateMethod == UpdateOrCoroutine.Coroutine) StartCoroutine(UpdateCoroutine());
    }

    private void Update() {
        if (!_initialized) return;
        if (updateMethod == UpdateOrCoroutine.Update) {
            foreach(Obstacle obs in _OBSTACLES) {
                if (obs.isDynamic && (obs.obstacle.hasChanged || obs.prevFriction != obs.frictionCoefficient)) UpdateObstacle(obs, obs.index);
            }
        }
        if (_PARTICLE_CONTROLLER != null) {
            _SHADER.Dispatch(_CLEAR_PROJECTIONS_KERNEL, Mathf.CeilToInt((float)_PARTICLE_CONTROLLER.numParticles / 64f), 1, 1);
            _SHADER.Dispatch(_CALCULATE_PROJECTIONS_KERNEL, Mathf.CeilToInt((float)_PARTICLE_CONTROLLER.numParticles / 128f), 1, 1);
        }
        if (_toggle_debug) PrintDebugBuffer(1000);
    }

    // A separate coroutine that act as our "update" loop.
    // The use of a coroutine prevents the system from crashing due to the lag inherent in updating obstacles, perhaps some that are very big to begin with.
    public IEnumerator UpdateCoroutine() {
        while(true) {
            foreach(Obstacle obs in _OBSTACLES) {
                if (obs.isDynamic && obs.obstacle.hasChanged) UpdateObstacle(obs, obs.index);
                yield return null;
            }
            yield return new WaitForSeconds(0.1f);
        }
    }

    // This is a function that's called from `UpdateObstacle` during the runtime, never before the first frame of the simulation
    // This function essentially updates the buffers based on vertex and triangle indexes of each obstacle
    // This is done via `SetData`, which has 4 arguments in this case:
    // 1. The data to be updated
    // 2. The first element index in data to copy to the compute buffer.
    // 3. The first element index in compute buffer to receive the data.
    // 4. The number of elements to copy.
    public void UpdateObstacleInCombined(Obstacle obs) {
        // Update vertices
        _VERTICES_BUFFER.SetData(obs.vertices, 0, obs.vertexIndices[0], obs.vertexIndices[1]);

        // Update triangles
        _TRIANGLES_BUFFER.SetData(obs.triangles.ToArray(), 0, obs.triangleIndices[0], obs.triangleIndices[1]);
        
        // Update bounds
        ObstacleBounds[] BS = new ObstacleBounds[1] {obs.bounds};
        _BOUNDS_BUFFER.SetData(BS, 0, obs.index, 1);
    }

    // Update is called once per frame
    public void ManualUpdate() {
        Debug.Log("Performing a manual update");
        if (!_initialized) {
            Debug.LogError("ERROR: has not been initialized yet. Please manually initialize prior to manually updating");
            return;
        }
        foreach(Obstacle obs in _OBSTACLES) {
            if (obs.isDynamic && obs.obstacle.hasChanged) UpdateObstacle(obs, obs.index);
        }
    }

    private void OnDestroy() {
        ReleaseBuffers();
    }

    private void ReleaseBuffers() {
        if (_VERTICES_BUFFER != null) {
            _VERTICES_BUFFER.Release();
            _VERTICES_BUFFER = null;
        }
        if (_TRIANGLES_BUFFER != null) {
            _TRIANGLES_BUFFER.Release();
            _TRIANGLES_BUFFER = null;
        }
        if (_BOUNDS_BUFFER != null) {
            _BOUNDS_BUFFER.Release();
            _BOUNDS_BUFFER = null;
        }
        if (_VERTEX_INDEX_LIMITS_BUFFER != null) {
            _VERTEX_INDEX_LIMITS_BUFFER.Release();
            _VERTEX_INDEX_LIMITS_BUFFER = null;
        }
        if (_TRIANGLE_INDEX_LIMITS_BUFFER != null) {
            _TRIANGLE_INDEX_LIMITS_BUFFER.Release();
            _TRIANGLE_INDEX_LIMITS_BUFFER = null;
        }
        if (_DEBUG_BUFFER != null) {
            _DEBUG_BUFFER.Release();
            _DEBUG_BUFFER = null;
        }
    }

    public void FullReset() {
        ReleaseBuffers();
        _OBSTACLES = new List<Obstacle>();
        _initialized = false;
    }

    public void PrintDebug(int numRecords = 1000) {
        int[] debugRecords = new int[_PARTICLE_CONTROLLER.numParticles];
        _DEBUG_BUFFER.GetData(debugRecords);

    }

    private void PrintDebugBuffer(int debugSize = 500) {
        int[] temp = new int[debugSize];
        _DEBUG_BUFFER.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"DEBUF BUFFER:\n{top}\n{bottom}");
        temp = null;
    }

    public void PrintObstacleDetails(int obsIndex = 0) {
        if (obsIndex < 0 || obsIndex >= _OBSTACLES.Count) {
            Debug.LogError($"ERROR: Cannot get obstacle index {obsIndex} b/c it's out of bounds of the `_OBSTACLES` list");
            return;
        }
        Obstacle obstacle = _OBSTACLES[obsIndex];
        Mesh mesh = (obstacle.obstacle.GetComponent<MeshFilter>() != null) 
            ? obstacle.obstacle.GetComponent<MeshFilter>().sharedMesh
            : obstacle.obstacle.GetComponent<SkinnedMeshRenderer>().sharedMesh;
        var vs = mesh.vertices;
        var ts = mesh.triangles;
        
        string vs_x_str = "";
        string vs_y_str = "";
        string vs_z_str = "";
        foreach(Vector3 v in vs) {
            Vector3 vw = obstacle.obstacle.TransformPoint(v);
            vs_x_str += $"{vw.x},";
            vs_y_str += $"{vw.y},";
            vs_z_str += $"{vw.z},";
        }
        vs_x_str = vs_x_str.Remove(vs_x_str.Length-1, 1);
        vs_y_str = vs_y_str.Remove(vs_y_str.Length-1, 1);
        vs_z_str = vs_z_str.Remove(vs_z_str.Length-1, 1);

        string ts_0_str = "";
        string ts_1_str = "";
        string ts_2_str = "";
        for(int i = 0; i < ts.Length; i+=3) {
            ts_0_str += $"{ts[i]},";
            ts_1_str += $"{ts[i+1]},";
            ts_2_str += $"{ts[i+2]},";
        }
        ts_0_str = ts_0_str.Remove(ts_0_str.Length-1, 1);
        ts_1_str = ts_1_str.Remove(ts_1_str.Length-1, 1);
        ts_2_str = ts_2_str.Remove(ts_2_str.Length-1, 1);
        
        Debug.Log("Vertices:");
        Debug.Log($"X: {vs_x_str}");
        Debug.Log($"Y: {vs_y_str}");
        Debug.Log($"Z: {vs_z_str}");

        Debug.Log("Triangles");
        Debug.Log($"0: {ts_0_str}");
        Debug.Log($"1: {ts_1_str}");
        Debug.Log($"2: {ts_2_str}");
    }
}
