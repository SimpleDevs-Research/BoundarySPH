using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;
using Unity.Mathematics;
using UnityEditor;
using OP = ObstaclePrimitives.Structs;

public class MeshObsGPU : MonoBehaviour
{

    // This one stays in the CPU and simply stores basic info about the obstacle. 
    // This one is NOT sent over to the GPU in any way
    [System.Serializable]
    public class TestObstacle {
        [ReadOnly] public int obstacleID;
        public MeshObs obstacle;
        public float mass = 1f;
        [Range(0f,1f)] public float frictionCoefficient = 0f;
        public bool checkObstacleBounds = true;
        public bool checkTriangleBounds = true;
        public bool enable_external_forces = false;
        public bool alwaysUpdate = false;
        [HideInInspector] public float prevFriction = 0f;
        [HideInInspector] public bool prevCheckObstacleBounds = true;
        [HideInInspector] public bool prevCheckTriangleBounds = true;

        public bool show_vertex_positions = true;
        public bool show_vertex_normals = true;
        public bool show_vertex_forces = true;
        public bool show_face_normals = true;
        public bool show_3D_edge_normals = true;
        public bool show_2D_edge_normals = true;
        public bool show_centroids = true;
        public bool show_bounds = true;
        public bool show_triangle_bounds = true;

        public bool show_vertices => show_vertex_positions || show_vertex_normals || show_vertex_forces;
        public bool show_triangles => show_face_normals || show_2D_edge_normals || show_centroids || show_triangle_bounds;
        public bool show_gizmos => show_vertices || show_triangles;
    }
    
    public ComputeShader _SHADER;
    public ParticleController _PARTICLE_CONTROLLER = null;
    public BoidsController _BOIDS_CONTROLLER = null;
    public List<TestObstacle> obstacles;

    public List<OP.ObstacleStatic> obstacles_static;
    public List<OP.TriangleStatic> triangles_static;
    public List<OP.EdgeStatic> edges_static;

    public float particleRenderRadius = 1f;
    public List<Transform> particles = new List<Transform>();
    
    [SerializeField] private float _dt = 0.0165f;

    [SerializeField] private bool drawObstacleGizmos = false;
    [SerializeField] private bool drawProjectionGizmos = false;
    private bool drawGizmos => drawObstacleGizmos || drawProjectionGizmos;

    private int numObstacles, numVertices, numTriangles, numEdges;
    [HideInInspector] public int numParticles;
    
    public ComputeBuffer obstacles_static_buffer, obstacles_dynamic_buffer;
    public ComputeBuffer vertices_static_buffer, vertices_dynamic_buffer;
    public ComputeBuffer triangles_static_buffer, triangles_dynamic_buffer;
    public ComputeBuffer edges_static_buffer, edges_dynamic_buffer;
    public ComputeBuffer translational_forces_buffer, torque_forces_buffer;
    public ComputeBuffer particles_buffer, projections_buffer;

    private int resetVertexForcesKernel, updateVerticesKernel, updateTrianglesKernel, updateEdgesKernel, resetHasChangedKernel;
    private int resetTempProjectionsKernel, checkForProjectionKernel, combineForcesKernel;

    public OP.ObstacleDynamic[] obstacles_dynamic_array;
    public OP.VertexDynamic[] vertices_dynamic_array;
    public OP.TriangleDynamic[] triangles_dynamic_array; 
    public float3[] edges_dynamic_array; 
    public int3[] translational_forces_array, torque_forces_array;
    public OP.Projection[] projections_array;

    public bool alwaysUpdateTransforms = false;

    [SerializeField] private bool printDebugs = true;
    private bool _initialized = false;

    void OnDrawGizmos() {
        
        if (!Application.isPlaying || !_initialized) return;
        if (!drawGizmos) return;

        obstacles_dynamic_buffer.GetData(obstacles_dynamic_array);
        triangles_dynamic_buffer.GetData(triangles_dynamic_array);
        vertices_dynamic_buffer.GetData(vertices_dynamic_array);
          
        if (drawObstacleGizmos) {
            OP.VertexDynamic v_dynamic;
            for(int i = 0; i < obstacles_static.Count; i++) {
                if (!obstacles[i].show_gizmos) continue;
                OP.ObstacleStatic o_static = obstacles_static[i];
                OP.ObstacleDynamic o_dynamic = obstacles_dynamic_array[i];
                // Render vertices
                if (obstacles[i].show_vertices) {
                    for(uint vi = o_static.vs[0]; vi < o_static.vs[0] + o_static.vs[1]; vi++) {
                        v_dynamic = vertices_dynamic_array[(int)vi];
                        Vector3 vp = new Vector3(v_dynamic.position[0], v_dynamic.position[1], v_dynamic.position[2]);
                        Vector3 vn = new Vector3(v_dynamic.normal[0],v_dynamic.normal[1],v_dynamic.normal[2]).normalized;
                        if (obstacles[i].show_vertex_positions) {
                            Gizmos.color = Color.red;
                            Gizmos.DrawSphere(vp,0.1f);
                        }
                        if (obstacles[i].show_vertex_normals) {
                            Handles.color = Color.blue;
                            Handles.DrawLine(vp, vp + vn*5f, 3);
                        }
                        if (obstacles[i].show_vertex_forces) {
                            // We show vertex forces just because...
                            Vector3 vf = new Vector3(v_dynamic.force[0], v_dynamic.force[1], v_dynamic.force[2]);
                            Handles.color = new Vector4(0f,0f,0f,0.5f);
                            Handles.DrawLine(vp,vp+vf,2);
                        }
                        
                    }
                }
                // Render triangles
                if (obstacles[i].show_triangles) {
                    for(uint ti = o_static.ts[0]; ti < o_static.ts[0] + o_static.ts[1]; ti++) {
                        // Draw the center
                        if (obstacles[i].show_centroids) {
                            Gizmos.color = Color.red;
                            Gizmos.DrawSphere(triangles_dynamic_array[ti].center,0.75f);
                        }
                        // Draw the 3D face normal
                        if (obstacles[i].show_face_normals) {
                            Handles.color = Color.blue;
                            Handles.DrawLine(
                                triangles_dynamic_array[ti].center,
                                triangles_dynamic_array[ti].center + triangles_dynamic_array[ti].normal*10f, 
                                3
                            );
                        }
                        // Draw the 2D edge normals
                        if (obstacles[i].show_2D_edge_normals) {
                            Handles.color = Color.yellow;
                            Handles.DrawLine(
                                triangles_dynamic_array[ti].center,
                                triangles_dynamic_array[ti].center + triangles_dynamic_array[ti].v1v2n*3f,
                                3
                            );
                            Handles.DrawLine(
                                triangles_dynamic_array[ti].center,
                                triangles_dynamic_array[ti].center + triangles_dynamic_array[ti].v2v3n*3f,
                                3
                            );
                            Handles.DrawLine(
                                triangles_dynamic_array[ti].center,
                                triangles_dynamic_array[ti].center + triangles_dynamic_array[ti].v1v3n*3f,
                                3
                            );
                        }
                        // Show bounds
                        if (obstacles[i].show_triangle_bounds) {
                            Gizmos.color = Color.grey;
                            Gizmos.DrawWireCube(
                                (triangles_dynamic_array[ti].upperBound + triangles_dynamic_array[ti].lowerBound)/2f,
                                triangles_dynamic_array[ti].upperBound - triangles_dynamic_array[ti].lowerBound
                            );
                        }
                        
                    }
                }
                
                // Render the edge normals
                if (obstacles[i].show_3D_edge_normals) {
                    Handles.color = Color.blue;
                    for(uint ei = o_static.es[0]; ei < o_static.es[0] + o_static.es[1]; ei++) {
                        o_static = obstacles_static[(int)(edges_static[(int)ei].obstacleIndex)];
                        float3 v1 = vertices_dynamic_array[(int)(o_static.vs[0] + edges_static[(int)ei].vertices[0])].position;
                        float3 v2 = vertices_dynamic_array[(int)(
                            (int)(
                                o_static.vs[0] + edges_static[(int)ei].vertices[1]
                            )
                        )].position;
                        float3 m = (v1+v2)/2f;
                        Vector3 midpoint = new Vector3(m[0],m[1],m[2]);
                        float3 n = edges_dynamic_array[ei];
                        Vector3 normal = new Vector3(n[0],n[1],n[2]);
                        Handles.DrawLine(midpoint, midpoint + normal, 3);
                    }
                }
                
                // Render bounds
                if (obstacles[i].show_bounds) {
                    Gizmos.color = Color.black;
                    Gizmos.DrawWireCube(
                        (o_dynamic.upperBound + o_dynamic.lowerBound)/2f, 
                        o_dynamic.upperBound - o_dynamic.lowerBound
                    );
                }
                
            }
            
        }
        
        if (drawProjectionGizmos) {
            projections_buffer.GetData(projections_array);
            OP.Particle[] debug_particles = new OP.Particle[particles_buffer.count];
            particles_buffer.GetData(debug_particles);
            for(int i = 0; i < numParticles; i++) {
                float3 n = new(
                    (float)projections_array[i].normal[0] / 1024f,
                    (float)projections_array[i].normal[1] / 1024f,
                    (float)projections_array[i].normal[2] / 1024f
                );
                Handles.color = Color.red;
                Handles.DrawLine(
                    debug_particles[i].position, 
                    debug_particles[i].position + n, 3);
                /*
                int tri_index = (int)projections_array[i].triangleID;
                if (tri_index == numTriangles) continue;
                OP.TriangleStatic t_static = triangles_static[tri_index];
                OP.TriangleDynamic t_dynamic = triangles_dynamic_array[tri_index];
                OP.ObstacleStatic o_static = obstacles_static[(int)t_static.obstacleIndex];

                // Rendering the projection
                Gizmos.color = Color.red;
                Gizmos.DrawSphere(projections_array[i].projection,particleRenderRadius*0.5f);
                // Rendering the closest point
                Gizmos.color = new Vector4(0f,0f,1f,0.25f);
                Gizmos.DrawSphere(projections_array[i].position,particleRenderRadius);
                // Rendering the normal vector
                Gizmos.color = Color.blue;
                Vector3 n = new Vector3((float)projections_array[i].normal[0] / 1024f, (float)projections_array[i].normal[1] / 1024f, (float)projections_array[i].normal[2] / 1024f);
                Gizmos.DrawRay(projections_array[i].position, n);
                
                // Rendering the external force felt by the particle
                Handles.color = Color.yellow;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].external_force, 2);
                // Rendering the particle force exerted by the particle
                Handles.color = Color.red;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].particle_force, 2);
                // Rendering the resulting force exerted by the obstacle onto the particle
                float3 pForce = projections_array[i].particle_force;
                float3 pN = new(
                    (float)projections_array[i].normal[0] / 1024f,
                    (float)projections_array[i].normal[1] / 1024f,
                    (float)projections_array[i].normal[2] / 1024f 
                );
                float3 externalForce = projections_array[i].external_force + (pForce - 1.5f * Unity.Mathematics.math.dot(pForce, pN) * pN);
                Handles.color = Color.black;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + externalForce, 2);
                // Calculating the resulting force
                Handles.color = Color.blue;
                float3 finalForce = pForce + externalForce;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + finalForce, 10);

                // Rendering the two edges associated and 2D and 3D normal vectors. Yellow = first edge, grey = second edge
                /*
                Handles.color = Color.yellow;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].e1,3);
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].e1_3DN,2);
                Handles.DrawLine(t_dynamic.center, t_dynamic.center + projections_array[i].e1_2DN,2);
                Handles.color = Color.grey;
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].e2,3);
                Handles.DrawLine(projections_array[i].position, projections_array[i].position + projections_array[i].e2_3DN,2);
                Handles.DrawLine(t_dynamic.center, t_dynamic.center + projections_array[i].e2_2DN,2);
                */
                /*
                int v1i = (int)(o_static.vs[0] + t_static.vertices[0]);
                int v2i = (int)(o_static.vs[0] + t_static.vertices[1]);
                int v3i = (int)(o_static.vs[0] + t_static.vertices[2]);
                */
                //Handles.color = Color.red;
                //Handles.DrawLine(
                //    vertices_dynamic_array[v1i].position,
                //    vertices_dynamic_array[v2i].position,
                //    3
                //);
                //Handles.DrawLine(
                //    vertices_dynamic_array[v1i].position,
                //    vertices_dynamic_array[v3i].position,
                //    3
                //);
                //Handles.DrawLine(
                //    vertices_dynamic_array[v2i].position,
                //    vertices_dynamic_array[v3i].position,
                //    3
                //);
                
                /*
                Gizmos.color = Color.yellow;
                Gizmos.DrawCube(
                    vertices_dynamic_array[v1i].position,
                    new Vector3(0.25f,0.25f,0.25f)
                );
                Gizmos.DrawCube(
                    vertices_dynamic_array[v2i].position,
                    new Vector3(0.25f,0.25f,0.25f)
                );
                Gizmos.DrawCube(
                    vertices_dynamic_array[v3i].position,
                    new Vector3(0.25f,0.25f,0.25f)
                );

                Handles.color = Color.black;
                Handles.DrawLine(projections_array[i].position, t_dynamic.center);
                */
            }
            
        }
    }

    // We run the setup on Start, not Awake, becuase ParticleController initializes on Awake
    void Start() {
        // Only run if we actually have a shader
        if (_SHADER == null) {
            if (printDebugs) Debug.LogError("MeshObsGPU - ERROR: Cannot preprocess obstacles due to missing reference to a compute shader");
            return;
        }
        // Don't run either if we don't have any particles
        if (_PARTICLE_CONTROLLER == null && particles.Count == 0) {
            if (printDebugs) Debug.LogError("MeshObsGPU - ERROR: Cannot preprocess obstacles due to missing particles");
            return;
        }
        // Check if we have a boids controller. If we do, we add them to our list of obstacles
        if (_BOIDS_CONTROLLER != null && _BOIDS_CONTROLLER.numBoids > 0) {
            Debug.Log("Adding Boids");
            obstacles.AddRange(_BOIDS_CONTROLLER.boids);
        } 
        PreprocessObstacles();
        UpdateObstacles(true);
        _initialized = true;
    }   

    void Update() {
        // Only run if we actually have a shader
        if (_SHADER == null) {
            if (printDebugs) Debug.LogError("MeshObsGPU - ERROR: Cannot update obstacles due to missing reference to a compute shader");
            return;
        }
        // Don't run either if we don't have any particles
        if (numParticles == 0) {
            if (printDebugs) Debug.LogError("MeshObsGPU - ERROR: Cannot update obstacles due to missing particles");
            return;
        }
        UpdateObstacles(alwaysUpdateTransforms);
        UpdateParticlePositions();
        // Reset projections
        _SHADER.Dispatch(resetTempProjectionsKernel, Mathf.CeilToInt((float)numParticles / 64f), 1, 1);
        // Get those projections
        _SHADER.Dispatch(
            checkForProjectionKernel, 
            Mathf.CeilToInt((float)numParticles / 32f), 
            Mathf.CeilToInt((float)numTriangles / 16f), 
            1
        );
        // Combine projection forces
        //_SHADER.Dispatch(combineForcesKernel, Mathf.CeilToInt((float)numParticles / 64f), 1, 1);

        translational_forces_buffer.GetData(translational_forces_array);
        torque_forces_buffer.GetData(torque_forces_array);
    }

    void FixedUpdate() {
        // Finally, we update each obstacle by getting the projections and checking if any obstacles should be updated, assuming they have a rigidbody
        //projections_buffer.GetData(projections_array);
        //translational_forces_buffer.GetData(translational_forces_array);
        //torque_forces_buffer.GetData(torque_forces_array);
        /*
        for(int i = 0; i < numObstacles; i++) {
            if (obstacles_static[i].has_rb == (int)0) continue;
            Vector3 translation_force = new Vector3(
                (float)(translational_forces_array[i][0] / 1024f),
                (float)(translational_forces_array[i][1] / 1024f),
                (float)(translational_forces_array[i][2] / 1024f)
            );
            Vector3 torque_force = new Vector3(
                (float)(torque_forces_array[i][0] / 1024f),
                (float)(torque_forces_array[i][1] / 1024f),
                (float)(torque_forces_array[i][2] / 1024f)
            );
            obstacles[i].obstacle.GetComponent<Rigidbody>().AddForce(translation_force);
            obstacles[i].obstacle.GetComponent<Rigidbody>().AddTorque(torque_force);
        }
        */
        /*
        for(int i = 0; i < numParticles; i++) {
            // Grab the projection
            OP.Projection p = projections_array[i];
            // Exit early if the projection's triangleID is == numTriangles
            if (p.triangleID == numTriangles) continue;
            // Get the associated triangle
            OP.TriangleStatic tri = triangles_static[(int)p.triangleID];
            // Get the associated obstacle
            TestObstacle obs = obstacles[(int)tri.obstacleIndex];
            // exit early if the associated obstacle doesn't have a rigidbody
            if (obs.obstacle.GetComponent<Rigidbody>() == null) continue;
            // If all else, apply the particle force onto the rigidbody at the position
            obs.obstacle.GetComponent<Rigidbody>().AddForceAtPosition(p.particle_force, p.position);
        }
        */
    }

    public void PreprocessObstacles() {

        // Initialize our lists 
        obstacles_static = new List<OP.ObstacleStatic>();
        List<OP.ObstacleDynamic> obstacles_dynamic = new List<OP.ObstacleDynamic>();
        List<OP.VertexStatic> vertices_static = new List<OP.VertexStatic>();
        //List<OP.VertexDynamic> vertices_dynamic = new List<OP.VertexDynamic>();
        List<OP.VertexDynamic> vertices_dynamic = new List<OP.VertexDynamic>();
        //List<OP.TriangleStatic> triangles_static = new List<OP.TriangleStatic>();
        triangles_static = new List<OP.TriangleStatic>();
        List<OP.TriangleDynamic> triangles_dynamic = new List<OP.TriangleDynamic>();
        edges_static = new List<OP.EdgeStatic>();
        List<float3> edges_dynamic = new List<float3>();
        
        // Preprocess our obstacles
        for(int i = 0; i < obstacles.Count; i++) {
            TestObstacle obstacle = obstacles[i];
            PreprocessObstacle(
                obstacle,
                i,
                ref obstacles_static, ref obstacles_dynamic,
                ref vertices_static, ref vertices_dynamic,
                ref triangles_static, ref triangles_dynamic,
                ref edges_static, ref edges_dynamic
            );
            obstacle.obstacle.position_transform.hasChanged = true;
            obstacle.obstacleID = i;
        }

        // Update our variables in the shader
        numObstacles = obstacles_static.Count;
        numVertices = vertices_static.Count;
        numTriangles = triangles_static.Count;
        numEdges = edges_static.Count;
        numParticles = (_PARTICLE_CONTROLLER != null) ? _PARTICLE_CONTROLLER.numParticles : particles.Count;
        if (_PARTICLE_CONTROLLER != null) particleRenderRadius = _PARTICLE_CONTROLLER.h;
        if (_PARTICLE_CONTROLLER != null) _dt = _PARTICLE_CONTROLLER.dt;
        _SHADER.SetInt("numObstacles", numObstacles);
        _SHADER.SetInt("numVertices", numVertices);
        _SHADER.SetInt("numTriangles", numTriangles);
        _SHADER.SetInt("numEdges", numEdges);
        _SHADER.SetInt("numParticles",numParticles);
        _SHADER.SetFloat("particleRenderRadius",particleRenderRadius);
        _SHADER.SetFloat("dt",_dt);

        // Grab our kernel references
        resetVertexForcesKernel = _SHADER.FindKernel("ResetVertexForces");
        updateVerticesKernel = _SHADER.FindKernel("UpdateVertices");
        updateEdgesKernel = _SHADER.FindKernel("UpdateEdges");
        updateTrianglesKernel = _SHADER.FindKernel("UpdateTriangles");
        resetHasChangedKernel = _SHADER.FindKernel("ResetHasChanged");
        checkForProjectionKernel = _SHADER.FindKernel("CheckForProjection2");
        resetTempProjectionsKernel = _SHADER.FindKernel("ResetTempProjections");
        combineForcesKernel = _SHADER.FindKernel("CombineForces");

        // Initialize our buffers
        obstacles_static_buffer = new ComputeBuffer(obstacles_static.Count, sizeof(uint)*9 + sizeof(float));
        obstacles_dynamic_buffer = new ComputeBuffer(obstacles_dynamic.Count, sizeof(uint)*4 + sizeof(float)*20);
        vertices_static_buffer = new ComputeBuffer(vertices_static.Count, sizeof(uint) + sizeof(float)*6);
        vertices_dynamic_buffer = new ComputeBuffer(vertices_dynamic.Count, sizeof(uint) + sizeof(float)*9);
        triangles_static_buffer = new ComputeBuffer(triangles_static.Count, sizeof(uint)*7 + sizeof(float)*9);
        triangles_dynamic_buffer = new ComputeBuffer(triangles_dynamic.Count, sizeof(uint)*7 + sizeof(float)*22);
        edges_static_buffer = new ComputeBuffer(edges_static.Count, sizeof(uint)*5 + sizeof(float)*6);
        edges_dynamic_buffer = new ComputeBuffer(edges_dynamic.Count, sizeof(float)*3);
        translational_forces_buffer = new ComputeBuffer(obstacles_static.Count, sizeof(int)*3);
        torque_forces_buffer = new ComputeBuffer(obstacles_static.Count, sizeof(int)*3);
        particles_buffer = (_PARTICLE_CONTROLLER != null) ? _PARTICLE_CONTROLLER.PARTICLES_BUFFER : new ComputeBuffer(numParticles, sizeof(float)*6);
        projections_buffer = (_PARTICLE_CONTROLLER != null) ? _PARTICLE_CONTROLLER.PROJECTIONS_BUFFER : new ComputeBuffer(numParticles, sizeof(uint) + sizeof(int)*8 + sizeof(float)*27);

        // Update the buffers
        obstacles_static_buffer.SetData(obstacles_static.ToArray());
        obstacles_dynamic_buffer.SetData(obstacles_dynamic.ToArray());
        vertices_static_buffer.SetData(vertices_static.ToArray());
        vertices_dynamic_buffer.SetData(vertices_dynamic.ToArray());
        triangles_static_buffer.SetData(triangles_static.ToArray());
        triangles_dynamic_buffer.SetData(triangles_dynamic.ToArray());
        edges_static_buffer.SetData(edges_static.ToArray());
        edges_dynamic_buffer.SetData(edges_dynamic.ToArray());
        UpdateParticlePositions();
        if (_PARTICLE_CONTROLLER != null) {
            projections_array = new OP.Projection[numParticles];
            projections_buffer.SetData(projections_array);
        }
        translational_forces_array = new int3[numObstacles];
        torque_forces_array = new int3[numObstacles];
        for(int i = 0; i < obstacles_static.Count; i++) {
            translational_forces_array[i] = new(0,0,0);
            torque_forces_array[i] = new(0,0,0);
        }
        translational_forces_buffer.SetData(translational_forces_array);
        torque_forces_buffer.SetData(torque_forces_array);

        // Associate our buffers with our kernels
        // Reset Vertex Forces
        _SHADER.SetBuffer(resetVertexForcesKernel, "vertices_dynamic", vertices_dynamic_buffer);
        // Update Vertices
        _SHADER.SetBuffer(updateVerticesKernel,"vertices_static", vertices_static_buffer);
        _SHADER.SetBuffer(updateVerticesKernel,"vertices_dynamic", vertices_dynamic_buffer);
        _SHADER.SetBuffer(updateVerticesKernel,"obstacles_static", obstacles_static_buffer);
        _SHADER.SetBuffer(updateVerticesKernel,"obstacles_dynamic", obstacles_dynamic_buffer);
        // Update Edges
        _SHADER.SetBuffer(updateEdgesKernel,"edges_static", edges_static_buffer);
        _SHADER.SetBuffer(updateEdgesKernel,"edges_dynamic", edges_dynamic_buffer);
        _SHADER.SetBuffer(updateEdgesKernel,"obstacles_dynamic", obstacles_dynamic_buffer);
        // Update Triangles
        _SHADER.SetBuffer(updateTrianglesKernel,"triangles_static",triangles_static_buffer);
        _SHADER.SetBuffer(updateTrianglesKernel,"triangles_dynamic",triangles_dynamic_buffer);
        _SHADER.SetBuffer(updateTrianglesKernel,"obstacles_dynamic",obstacles_dynamic_buffer);
        _SHADER.SetBuffer(updateTrianglesKernel,"obstacles_static",obstacles_static_buffer);
        _SHADER.SetBuffer(updateTrianglesKernel,"vertices_dynamic",vertices_dynamic_buffer);
        // Reset Has Changed
        _SHADER.SetBuffer(resetHasChangedKernel,"obstacles_dynamic",obstacles_dynamic_buffer);
        // Reset Temp Projections
        _SHADER.SetBuffer(resetTempProjectionsKernel,"projections", projections_buffer);
        _SHADER.SetBuffer(resetTempProjectionsKernel,"particles",particles_buffer);
        _SHADER.SetBuffer(resetTempProjectionsKernel,"translational_forces",translational_forces_buffer);
        _SHADER.SetBuffer(resetTempProjectionsKernel,"torque_forces",torque_forces_buffer);
        // Check For Projections
        _SHADER.SetBuffer(checkForProjectionKernel,"particles",particles_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"projections",projections_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"triangles_dynamic",triangles_dynamic_buffer);
        //_SHADER.SetBuffer(checkForProjectionKernel,"triangles_static",triangles_static_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"obstacles_static",obstacles_static_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"obstacles_dynamic",obstacles_dynamic_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"vertices_dynamic",vertices_dynamic_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"edges_dynamic",edges_dynamic_buffer);
        _SHADER.SetBuffer(checkForProjectionKernel,"translational_forces",translational_forces_buffer);
        // Combine Forces Kernel
        _SHADER.SetBuffer(combineForcesKernel,"projections",projections_buffer);
        _SHADER.SetBuffer(combineForcesKernel,"triangles_static",triangles_static_buffer);
        _SHADER.SetBuffer(combineForcesKernel,"obstacles_static",obstacles_static_buffer);
        _SHADER.SetBuffer(combineForcesKernel,"obstacles_dynamic",obstacles_dynamic_buffer);
        _SHADER.SetBuffer(combineForcesKernel,"translational_forces",translational_forces_buffer);
        _SHADER.SetBuffer(combineForcesKernel,"torque_forces",torque_forces_buffer);

        // Finally, prepare our global arrays
        obstacles_dynamic_array = new OP.ObstacleDynamic[numObstacles];
        vertices_dynamic_array = new OP.VertexDynamic[numVertices];
        triangles_dynamic_array = new OP.TriangleDynamic[numTriangles];
        projections_array = new OP.Projection[numParticles];
        edges_dynamic_array = new float3[numEdges];

        if (_BOIDS_CONTROLLER != null) _BOIDS_CONTROLLER.SetExternalForcesBuffer(translational_forces_buffer);
    }

    void PreprocessObstacle(
        TestObstacle obstacle, 
        int index,
        ref List<OP.ObstacleStatic> obstacles_static, ref List<OP.ObstacleDynamic> obstacles_dynamic,
        ref List<OP.VertexStatic> vertices_static, ref List<OP.VertexDynamic> vertices_dynamic,
        ref List<OP.TriangleStatic> triangles_static, ref List<OP.TriangleDynamic> triangles_dynamic,
        ref List<OP.EdgeStatic> edges_static, ref List<float3> edges_dynamic
    ) {
        
        // Before anything, we need to extract the mesh data! 
        // We also need to extract the vertex anbd triangle data from that mesh
        Mesh mesh = obstacle.obstacle.GetMesh();
        var vs = obstacle.obstacle.vertices;
        var ts = obstacle.obstacle.triangles;
        obstacle.prevFriction = obstacle.frictionCoefficient;
        obstacle.prevCheckObstacleBounds = obstacle.checkObstacleBounds;
        obstacle.prevCheckTriangleBounds = obstacle.checkTriangleBounds;

        // We initialize the details for the obstacles themselves
        OP.ObstacleStatic o_static = new OP.ObstacleStatic();
        OP.ObstacleDynamic o_dynamic = new OP.ObstacleDynamic();
        o_static.index = (uint)index;
        o_dynamic.index = (uint)index;
        o_static.mass = obstacle.mass;
        o_static.has_rb = (obstacle.enable_external_forces) ? (uint)1 : (uint)0;
        o_static.has_smr = (obstacle.obstacle.hasSkinnedMeshFilter) ? (uint)1 : (uint)0;
        o_dynamic.frictionCoefficient = obstacle.frictionCoefficient;
        o_dynamic.checkObstacleBounds = (obstacle.checkObstacleBounds) ? (uint)1 : (uint)0;
        o_dynamic.checkTriangleBounds = (obstacle.checkTriangleBounds) ? (uint)1 : (uint)0;
        o_dynamic.hasChanged = 1;

        // ===== GENERATING VERTICES DATA ===== //

        // We need to condense the vertices into a smaller list
        // We do this by creating:
        //  1) `fixed_vs` : List<float3> = A list of vertices simialr to `vs` but without duplicates
        //  2) `vs_map`: int[] = A int-int mapping between `vs` and `fixed_vs`
        //  3) `vs_static` : List<OP.VertexStatic> = the static details of each vertex in `fixed_vs`
        //  4) `vs_dynamic` : List<OP.VertexDynamic> = the dynamic details of each vertex in `fixed_vs`
        //  5) `lowerBound` : float3 = the lower bounds of the obstacle
        //  6) `upperBound` : float3 = the upper bounds of the obstacle 
        //  7) `vertex_normals` : List<float3> = stores normal vectors for vertices. Will be added at the end.
        obstacle.obstacle.fixed_vs = new List<float3>();
        obstacle.obstacle.vs_map = new uint[vs.Length];
        List<OP.VertexStatic> vs_static = new List<OP.VertexStatic>();
        List<OP.VertexDynamic> vs_dynamic = new List<OP.VertexDynamic>();
        //float3 lowerBound = new(0f,0f,0f);
        //float3 upperBound = new(0f,0f,0f);
        List<float3> vertex_normals = new List<float3>();

        // We can now iterate through each verex in `vs` and filter out duplicates
        int mi;
        float3 localPosition;
        for(int i = 0; i < vs.Length; i++) {
            // We save the current vertex's local position
            mi = obstacle.obstacle.fixed_vs.IndexOf(vs[i]);
            // If the mapped index is -1, then there's no entry in `fixed_vs`
            if (mi == -1) {
                // Grab the current index, and add the local position to `fixed_vs`
                mi = obstacle.obstacle.fixed_vs.Count;
                obstacle.obstacle.fixed_vs.Add(vs[i]);
                // Initialize the corresponding entriers in `vs_static` and `vs_dynamic`
                OP.VertexStatic v_static = new OP.VertexStatic();
                OP.VertexDynamic v_dynamic = new OP.VertexDynamic();
                // We set the obstacle index to the current index
                v_static.obstacleIndex = (uint)index;
                v_dynamic.obstacleIndex = (uint)index;
                // We set the local position of the vertex
                localPosition = new(vs[i].x, vs[i].y, vs[i].z);
                if (Mathf.Abs(localPosition[0]) < 0.000000000001f) localPosition[0] = 0f;
                if (Mathf.Abs(localPosition[1]) < 0.000000000001f) localPosition[1] = 0f;
                if (Mathf.Abs(localPosition[2]) < 0.000000000001f) localPosition[2] = 0f;
                v_static.localPosition = localPosition;
                // We set the ending point to the local position, currently. When determining the world-space normal in the update loop, we'll use this to determine the normal direction
                //v_static.localNormalEnd = localPosition;
                vertex_normals.Add(localPosition);
                //v_dynamic.normal = new(0f,0f,0f);
                // Add `v_static` into `vs_static` and `v_dynamic` into `vs_dynamic`
                vs_static.Add(v_static);
                vs_dynamic.Add(v_dynamic);
                // Update the lower and upper bounds
                //lowerBound[0] = Mathf.Min(lowerBound[0],localPosition[0]);
                //lowerBound[1] = Mathf.Min(lowerBound[1],localPosition[1]);
                //lowerBound[2] = Mathf.Min(lowerBound[2],localPosition[2]);
                //upperBound[0] = Mathf.Max(upperBound[0],localPosition[0]);
                //upperBound[1] = Mathf.Max(upperBound[1],localPosition[1]);
                //upperBound[2] = Mathf.Max(upperBound[2],localPosition[2]);
            }
            // Now, we merely have to update `vs_map`
            obstacle.obstacle.vs_map[i] = (uint)mi;
        }
        // We update `o_static` with number of filtered vertices
        // we know [0] because it's merely the current count of `vertices_static`
        o_static.vs = new((uint)vertices_static.Count,(uint)obstacle.obstacle.fixed_vs.Count);

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
        Dictionary<uint2,int> es_static_map = new Dictionary<uint2,int>();

        // Need to iterate through each triangle
        uint ti;
        for(uint i = 0; i < (uint)ts.Length; i+=3) {
            // Index for current triangle
            ti = i/3;
            // Generate a new triangle at the approriate index for both `ts_static` and `ts_dynamic`
            ts_static[ti] = new OP.TriangleStatic();
            ts_dynamic[ti] = new OP.TriangleDynamic();
            ts_static[ti].obstacleIndex = (uint)index;
            ts_dynamic[ti].obstacleIndex = (uint)index;

            // We get the corresponding indices of `vs_fixed` for this triangle. We need to store them in `ts_static[ti]`
            // In case it's confusing, know that `ts[i]` is an integer index for the original entry in `vs`.
            //  We have to get the corresponding index in `fixed_vs`. Hence `vs_map[ts[i]]`
            ts_static[(int)ti].vertices = new(
                obstacle.obstacle.vs_map[ts[(int)i]], 
                obstacle.obstacle.vs_map[ts[(int)i+1]], 
                obstacle.obstacle.vs_map[ts[(int)i+2]]
            );
            ts_dynamic[(int)ti].vertices = ts_static[(int)ti].vertices;
            
            // We need to calculate the angles for each vertex. To do so, we need to reference the local positions of each vertex.
            float3 v1f = vs_static[(int)ts_static[(int)ti].vertices[0]].localPosition;
            float3 v2f = vs_static[(int)ts_static[(int)ti].vertices[1]].localPosition;
            float3 v3f = vs_static[(int)ts_static[(int)ti].vertices[2]].localPosition;
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
            //if (localNormal.magnitude == 0) {
            //    Debug.Log($"{ts[(int)i]}|{ts[(int)i+1]}|{ts[(int)i+2]}\t{vs[ts[(int)i]]}|{vs[ts[(int)i+1]]}|{vs[ts[(int)i+2]]}\t{v1}, {v2}, {v3}");
            //}
            float3 localNormalF = new(localNormal.x,localNormal.y,localNormal.z);
            ts_static[ti].localNormal = localNormalF;

            // We need to generate the normal vectors for our vertices too. We do this by, for each vertex, adding the normal weighed by the influence angle
            //vs_static[(int)ts_static[(int)ti].vertices[0]].localNormalEnd += localNormalF * ts_static[(int)ti].angles[0];
            vertex_normals[(int)ts_static[(int)ti].vertices[0]] += localNormalF * ts_static[(int)ti].angles[0];
            //vs_static[(int)ts_static[(int)ti].vertices[1]].localNormalEnd += localNormalF * ts_static[(int)ti].angles[1];
            vertex_normals[(int)ts_static[(int)ti].vertices[1]] += localNormalF * ts_static[(int)ti].angles[1];
            //vs_static[(int)ts_static[(int)ti].vertices[2]].localNormalEnd += localNormalF * ts_static[(int)ti].angles[2];
            vertex_normals[(int)ts_static[(int)ti].vertices[2]] += localNormalF * ts_static[(int)ti].angles[2];
            
            // Now we have to update `es_static` and `es_dynamic`
            // The honestly crummy thing is that int2 can't be used as a key... if it's a list
            // So we used a dictionary `es_static_map` to temporarily link int2 pairs to edges in `es_static`

            int es_map_index_v1v2, es_map_index_v1v3, es_map_index_v2v3;
            // v1v2
            uint2 edge1_index = (ts_static[(int)ti].vertices[0] < ts_static[(int)ti].vertices[1])
                ? new(ts_static[(int)ti].vertices[0], ts_static[(int)ti].vertices[1])
                : new(ts_static[(int)ti].vertices[1], ts_static[(int)ti].vertices[0]);
            // If `v1v2` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge1_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v1v2 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = (uint)index;
                e_static.vertices = edge1_index;
                e_static.triangles = new(ti,(uint)ts.Length);
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
            uint2 edge2_index = (ts_static[(int)ti].vertices[0] < ts_static[(int)ti].vertices[2])
                ? new(ts_static[(int)ti].vertices[0], ts_static[(int)ti].vertices[2])
                : new(ts_static[(int)ti].vertices[2], ts_static[(int)ti].vertices[0]);
            // If `v1v3` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge2_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v1v3 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = (uint)index;
                e_static.vertices = edge2_index;
                e_static.triangles = new(ti,(uint)ts.Length);
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
            uint2 edge3_index = (ts_static[(int)ti].vertices[1] < ts_static[(int)ti].vertices[2])
                ? new(ts_static[(int)ti].vertices[1], ts_static[(int)ti].vertices[2])
                : new(ts_static[(int)ti].vertices[2], ts_static[(int)ti].vertices[1]);
            // If `v2v3` is not inside `es_static_map`, then it's also not inside `es_static`.
            if (!es_static_map.ContainsKey(edge3_index)) {
                // We create new corresponding entries in both `es_static`, `es_dynamic`, and `es_static_map`
                es_map_index_v2v3 = es_static.Count;
                OP.EdgeStatic e_static = new OP.EdgeStatic();
                e_static.obstacleIndex = (uint)index;
                e_static.vertices = edge3_index;
                e_static.triangles = new(ti,(uint)ts.Length);
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
            ts_static[(int)ti].edges = new((uint)es_map_index_v1v2,(uint)es_map_index_v1v3,(uint)es_map_index_v2v3);
            // We populate the ts_dynamic with our edges too
            ts_dynamic[(int)ti].edges = ts_static[(int)ti].edges;
        }

        // At this point, we have to attribute the vertex normals to their respective localNormalEnds in each item in vs_static
        for(int i = 0; i < vs_static.Count; i++) {
            OP.VertexStatic v_static = new OP.VertexStatic();
            v_static.obstacleIndex = vs_static[i].obstacleIndex;
            v_static.localPosition = vs_static[i].localPosition;
            v_static.localNormal = Unity.Mathematics.math.normalize(vertex_normals[i]);
            vs_static[i] = v_static;
        }
        vertices_static.AddRange(vs_static);
        vertices_dynamic.AddRange(vs_dynamic);

        // After all that, we need to update `o_static` on the updates to triangles and edges
        // [0] is the current count of `triangles_static` and `edges_static`, respectively
        o_static.ts = new((uint)triangles_static.Count, (uint)ts_static.Length);
        o_static.es = new((uint)edges_static.Count, (uint)es_static.Count);
        
        // We update the global triangles array
        triangles_static.AddRange(ts_static);
        triangles_dynamic.AddRange(ts_dynamic);
        
        // We have to iterate through all edges, find their midpoints, and their localNormalEnd`s
        float3 midpoint,normal;
        for(int ei = 0; ei < es_static.Count; ei++) {
            uint2 vertices = es_static[ei].vertices;
            midpoint = (vs_static[(int)vertices[0]].localPosition + vs_static[(int)vertices[1]].localPosition)/2f;
            // We calculate the normal by approximating the normals of each triangle assoicated with this edge. We have to be careful about edges with just one triangle (open edges).
            uint2 triangles = es_static[ei].triangles;
            normal = ts_static[triangles[0]].localNormal;
            if (triangles[1] != (uint)ts.Length) normal += (ts_static[triangles[1]].localNormal);
            // We can now update the midpoint and localNormalEnd of this edge
            OP.EdgeStatic e_static = new OP.EdgeStatic();
            e_static.obstacleIndex = es_static[ei].obstacleIndex;
            e_static.vertices = es_static[ei].vertices;
            e_static.triangles = es_static[ei].triangles;
            e_static.midpoint = midpoint;
            e_static.localNormal = Unity.Mathematics.math.normalize(normal);
            es_static[ei] = e_static;
        }

        // And we add the edges to our global list of edges
        edges_static.AddRange(es_static);
        edges_dynamic.AddRange(es_dynamic);

        // Finally, we append `o_static` and `o_dynamic` to `obstacles_static` and `obstacles_dynamic`, respectively
        obstacles_static.Add(o_static);
        obstacles_dynamic.Add(o_dynamic);
    }

    public void UpdateObstacles(bool forceUpdate = false) {
        MeshObs mo;
        bool needsUpdating = false;

        // We always reset the vertex forces
        _SHADER.Dispatch(resetVertexForcesKernel, Mathf.CeilToInt((float)numVertices/64f),1,1);

        for(int i = 0; i < obstacles.Count; i++) {
            // Check if the transform has been changed in any way
            mo = obstacles[i].obstacle; 
            if (
                    obstacles[i].alwaysUpdate
                    || mo.position_transform.hasChanged 
                    || obstacles[i].prevFriction != obstacles[i].frictionCoefficient 
                    || obstacles[i].prevCheckObstacleBounds != obstacles[i].checkObstacleBounds
                    || obstacles[i].prevCheckTriangleBounds != obstacles[i].checkTriangleBounds
                    || forceUpdate
            ) {
                //  We need to update the associated TestObjectDynamic
                obstacles_dynamic_array[i].index = (uint)i;
                obstacles_dynamic_array[i].position = new(mo.position_transform.position.x, mo.position_transform.position.y, mo.position_transform.position.z);
                obstacles_dynamic_array[i].rotation = new(mo.position_transform.rotation.x, mo.position_transform.rotation.y, mo.position_transform.rotation.z, mo.position_transform.rotation.w);
                obstacles_dynamic_array[i].scale = new(mo.position_transform.lossyScale.x, mo.position_transform.lossyScale.y, mo.position_transform.lossyScale.z);
                obstacles_dynamic_array[i].lowerBound = obstacles_dynamic_array[i].position;
                obstacles_dynamic_array[i].upperBound = obstacles_dynamic_array[i].position;
                //Rigidbody rb = obstacles[i].obstacle.rb;
                //if (rb != null) obstacles_dynamic_array[i].centerOfMass = new(rb.centerOfMass.x, rb.centerOfMass.y, rb.centerOfMass.z);
                obstacles_dynamic_array[i].frictionCoefficient = obstacles[i].frictionCoefficient;
                obstacles_dynamic_array[i].checkObstacleBounds = (obstacles[i].checkObstacleBounds) ? (uint)1 : (uint)0;
                obstacles_dynamic_array[i].checkTriangleBounds = (obstacles[i].checkTriangleBounds) ? (uint)1 : (uint)0;
                obstacles_dynamic_array[i].hasChanged = 1;
                obstacles[i].prevFriction = obstacles[i].frictionCoefficient;
                obstacles[i].prevCheckTriangleBounds = obstacles[i].checkTriangleBounds;
                if (obstacles[i].obstacle.hasSkinnedMeshFilter) {
                    // We have to grab the vertices manually
                    float3[] vs = obstacles[i].obstacle.GetWorldVertices();
                    if (vs != null) {
                        // We can update the vertices in the VS buffer
                        OP.VertexDynamic[] vd = new OP.VertexDynamic[numVertices];
                        vertices_dynamic_buffer.GetData(vd);
                        for(int vi = 0; vi < vs.Length; vi++) {
                            int vi2 = vi + (int)obstacles_static[i].vs[0];
                            vd[vi2].force = obstacles_static[i].mass * (vs[vi] - vd[vi2].position)/_dt;
                            vd[vi2].position = vs[vi];
                        }
                        vertices_dynamic_buffer.SetData(vd);
                    }
                }
                needsUpdating = true;
            }
        }

        if (needsUpdating) {
            if (printDebugs) {
                if (forceUpdate) Debug.Log("Forced Update!");
                else Debug.Log("UPDATING!");
            }
            // Before anything, we have to update our Obstacles Dynamic buffer
            obstacles_dynamic_buffer.SetData(obstacles_dynamic_array);

            // First, let's update the vertices
            _SHADER.Dispatch(updateVerticesKernel, Mathf.CeilToInt((float)numObstacles / 16f), 1, 1);
            // Second, let's reset the edge normals
            _SHADER.Dispatch(updateEdgesKernel, Mathf.CeilToInt((float)numEdges / 64f), 1, 1);
            // Thirdly, let's update the triangles
            _SHADER.Dispatch(updateTrianglesKernel, Mathf.CeilToInt((float)numTriangles / 64f), 1, 1);
            // Finally reset
            _SHADER.Dispatch(resetHasChangedKernel, Mathf.CeilToInt((float)numObstacles / 64f), 1, 1);
            for(int i = 0; i < obstacles.Count; i++) {
                obstacles[i].obstacle.position_transform.hasChanged = false;
            }
        }
    }

    void UpdateParticlePositions() {
        // We only run this if _PARTICLE_CONTROLLER is null
        if (_PARTICLE_CONTROLLER != null) return;
        OP.Particle[] particles = new OP.Particle[numParticles];
        for(int i = 0; i < numParticles; i++) {
            particles[i].position = particles[i].position;
        }
        if (particles_buffer != null) particles_buffer.SetData(particles);
    }

    void OnDestroy() {
        if (obstacles_static_buffer != null) obstacles_static_buffer.Release();
        if (obstacles_dynamic_buffer != null) obstacles_dynamic_buffer.Release();
        if (vertices_static_buffer != null) vertices_static_buffer.Release();
        if (vertices_dynamic_buffer != null) vertices_dynamic_buffer.Release();
        if (triangles_static_buffer != null) triangles_static_buffer.Release();
        if (triangles_dynamic_buffer != null) triangles_dynamic_buffer.Release();
        if (edges_static_buffer != null) edges_static_buffer.Release();
        if (edges_dynamic_buffer != null) edges_dynamic_buffer.Release();
        if (translational_forces_buffer != null) translational_forces_buffer.Release();
        if (torque_forces_buffer != null) torque_forces_buffer.Release();
        if (_PARTICLE_CONTROLLER == null && particles_buffer != null) particles_buffer.Release();
        if (_PARTICLE_CONTROLLER == null && projections_buffer != null) projections_buffer.Release();
    }
    
    // https://forum.unity.com/threads/whats-the-math-behind-transform-transformpoint.107401/
    public static float3 LocalPointToWorldPoint(float3 pos, float4 rot, float3 scale, float3 localPoint) {
        float3 s = new(localPoint[0] * scale[0], localPoint[1] * scale[1], localPoint[2] * scale[2]);
        return RotMultVec3(rot, s) + pos;
    }
    // https://answers.unity.com/questions/372371/multiply-quaternion-by-vector3-how-is-done.html
    public static float3 RotMultVec3(float4 quat, float3 vec){
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
        return new(
            (1f - (num5 + num6)) * vec[0] + (num7 - num12) * vec[1] + (num8 + num11) * vec[2],
            (num7 + num12) * vec[0] + (1f - (num4 + num6)) * vec[1] + (num9 - num10) * vec[2],
            (num8 - num11) * vec[0] + (num9 + num10) * vec[1] + (1f - (num4 + num5)) * vec[2]
        );
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
