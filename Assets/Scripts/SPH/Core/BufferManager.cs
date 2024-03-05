using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;

public class BufferManager : MonoBehaviour
{
    [Header("=== CONTROLLERS and RENDERERS ===")]
    public ParticleController _PARTICLE_CONTROLLER;
    public BoidsController _BOIDS_CONTROLLER;
    public MeshObsGPU _OBSTACLES_CONTROLLER;
    public PressureRenderer _PRESSURE_RENDERER;
    public ContinuousFlow _CF;
    public RecordingManager _RM;

    [Header("=== PARTICLE-related BUFFERS === ")]
    public ComputeBuffer PARTICLES_BUFFER;                      // Stores particles' position data
    public ComputeBuffer PARTICLES_VELOCITIES_BUFFER;           // Stores particles' velocity data
    public ComputeBuffer PARTICLES_DENSITIES_BUFFER;            // Stores particles' density data
    public ComputeBuffer PARTICLES_PRESSURES_BUFFER;            // Stores particles' pressure data (within their kernel)
    public ComputeBuffer PARTICLES_PRESSURE_FORCES_BUFFER;      // Stores the pressure force experienced BY a particle based on its kernel radius
    public ComputeBuffer PARTICLES_VISCOSITY_FORCES_BUFFER;     // Stores the viscosity force experienced BY a particle based on its kernel radius
    public ComputeBuffer PARTICLES_EXTERNAL_FORCES_BUFFER;      // Stores the external forces experienced BY a particle due to forces such as gravity

    [Header("=== GRID-related BUFFERS")]
    public ComputeBuffer PARTICLES_GRID_BUFFER;                 // The grid representation for particle position partitioning and hashing
    public ComputeBuffer PARTICLES_PRESSURE_GRID_BUFFER;        // The grid representation for particle pressure partitioning and hashing
    public ComputeBuffer MESHOBS_GRID_BUFFER;                   // The grid representation for all obstacles (that are boids)

    [Header("=== MESH OBSTACLE-related BUFFERS")]
    public ComputeBuffer MESHOBS_REARRANGED_OBSTACLES_BUFFER;   // Stores the rearranged indices of obstacles
    public ComputeBuffer MESHOBS_OBSTACLES_STATIC_BUFFER;       // Stores the static characteristics of mesh obstacles
    public ComputeBuffer MESHOBS_OBSTACLES_DYNAMIC_BUFFER;      // Stores the dynamic characteristics of mesh obstacles
    public ComputeBuffer MESHOBS_TRIANGLES_STATIC_BUFFER;       // Stores the static characteristics of mesh triangles
    public ComputeBuffer MESHOBS_TRIANGLES_DYNAMIC_BUFFER;      // Stores the dynamic characteristics of mesh triangles
    public ComputeBuffer MESHOBS_VERTICES_STATIC_BUFFER;        // Stores the static characteristics of mesh vertices
    public ComputeBuffer MESHOBS_VERTICES_DYNAMIC_BUFFER;       // Stores the dynamic characteristics of mesh vertices
    public ComputeBuffer MESHOBS_EDGES_STATIC_BUFFER;           // Stores the static characteristics of mesh vertices
    public ComputeBuffer MESHOBS_EDGES_DYNAMIC_BUFFER;          // Stores the dynamic characteristics of mesh edges
    public ComputeBuffer MESHOBS_TRANSLATION_FORCES_BUFFER;     // Stores the translational forces experienced BY an obstacle due to forces such as gravity or pressure fields
    public ComputeBuffer MESHOBS_TORQUE_FORCES_BUFFER;          // Stores the rotational forces experienced BY an obstacle due to forces such as gravity or pressure fields
    public ComputeBuffer MESHOBS_VELOCITIES_BUFFER;             // Stores the current velocity of the obstacle

    [Header("=== DEBUG SETTINGS ===")]
    [SerializeField] private bool _verbose = true;

    [Header("= PARTICLES-related ARRAYS =")]
    [SerializeField] private OP.Particle[] _particles_array;
    public OP.Particle[] particles_array => _particles_array;
    [SerializeField] private int[] _particles_grid_array;
    public int[] particles_grid_array => _particles_grid_array;
    [SerializeField] private float3[] _particles_velocities_array;
    public float3[] particles_velocities_array => _particles_velocities_array;
    [SerializeField] private float[] _particles_densities_array;
    public float[] particles_densities_array => _particles_densities_array;
    [SerializeField] private float[] _particles_pressures_array;
    public float[] particles_pressures_array => _particles_pressures_array;
    [SerializeField] private OP.Projection[] _particles_external_forces_array;
    public OP.Projection[] particles_external_forces_array => _particles_external_forces_array;
    
    [Header("= MESHOBS-related ARRAYS =")]
    [SerializeField] private int[] _obstacles_grid_array;
    public int[] obstacles_grid_array => _obstacles_grid_array;
    [SerializeField] private uint[] _rearranged_obstacles_array;
    [SerializeField] private OP.ObstacleStatic[] _obstacles_static_array;
    public OP.ObstacleStatic[] obstacles_static_array => _obstacles_static_array;
    [SerializeField] private OP.ObstacleDynamic[] _obstacles_dynamic_array;
    public OP.ObstacleDynamic[] obstacles_dynamic_array => _obstacles_dynamic_array;
    [SerializeField] private OP.TriangleStatic[] _triangles_static_array;
    public OP.TriangleStatic[] triangles_static_array => _triangles_static_array;
    [SerializeField] private OP.TriangleDynamic[] _triangles_dynamic_array;
    public OP.TriangleDynamic[] triangles_dynamic_array => _triangles_dynamic_array;
    [SerializeField] private OP.VertexStatic[] _vertices_static_array;
    public OP.VertexStatic[] vertices_static_array => _vertices_static_array;
    [SerializeField] private OP.VertexDynamic[] _vertices_dynamic_array;
    public OP.VertexDynamic[] vertices_dynamic_array => _vertices_dynamic_array;
    [SerializeField] private OP.EdgeStatic[] _edges_static_array;
    public OP.EdgeStatic[] edges_static_array => _edges_static_array;
    [SerializeField] private float3[] _edges_dynamic_array;
    public float3[] edges_dynamic_array => _edges_dynamic_array;
    [SerializeField] private int3[] _translation_forces_array;
    public int3[] translation_forces_array => _translation_forces_array;

    public void InitializeParticleBuffers(int numParticles = 1, int numGridCells = 1) {
        PARTICLES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*6+sizeof(int));
        PARTICLES_GRID_BUFFER = new ComputeBuffer(numGridCells, sizeof(int));

        PARTICLES_VELOCITIES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_DENSITIES_BUFFER = new ComputeBuffer(numParticles, sizeof(float));
        PARTICLES_PRESSURES_BUFFER = new ComputeBuffer(numParticles, sizeof(float));

        PARTICLES_PRESSURE_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_VISCOSITY_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_EXTERNAL_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(uint) + sizeof(int) + sizeof(float)*35);

        _particles_array = new OP.Particle[Mathf.Min(_particles_array.Length,numParticles)];
        _particles_grid_array = new int[Mathf.Min(_particles_grid_array.Length, numGridCells)];
        _particles_velocities_array = new float3[Mathf.Min(_particles_velocities_array.Length, numParticles)];
        _particles_densities_array = new float[Mathf.Min(_particles_densities_array.Length, numParticles)];
        _particles_pressures_array = new float[Mathf.Min(_particles_pressures_array.Length, numParticles)];
        _particles_external_forces_array = new OP.Projection[Mathf.Min(_particles_external_forces_array.Length, numParticles)];

        if (_verbose) Debug.Log("[BUFFER MANAGER] Particle buffers initialized!");
    }

    public void InitializeMeshObsBuffers(int numGridCells = 1, int numObstacles = 1, int numTriangles = 1, int numVertices = 1, int numEdges = 1) {
        MESHOBS_GRID_BUFFER = new ComputeBuffer(numGridCells, sizeof(int));
        MESHOBS_REARRANGED_OBSTACLES_BUFFER = new ComputeBuffer(numObstacles, sizeof(uint));

        MESHOBS_OBSTACLES_STATIC_BUFFER = new ComputeBuffer(numObstacles, sizeof(uint)*9 + sizeof(float)*7);
        MESHOBS_OBSTACLES_DYNAMIC_BUFFER = new ComputeBuffer(numObstacles, sizeof(uint)*15 + sizeof(float)*22);
        MESHOBS_TRIANGLES_STATIC_BUFFER = new ComputeBuffer(numTriangles, sizeof(uint)*7 + sizeof(float)*9);
        MESHOBS_TRIANGLES_DYNAMIC_BUFFER = new ComputeBuffer(numTriangles, sizeof(uint)*7 + sizeof(float)*22);
        MESHOBS_VERTICES_STATIC_BUFFER = new ComputeBuffer(numVertices, sizeof(uint) + sizeof(float)*6);
        MESHOBS_VERTICES_DYNAMIC_BUFFER = new ComputeBuffer(numVertices, sizeof(uint) + sizeof(float)*9);
        MESHOBS_EDGES_STATIC_BUFFER = new ComputeBuffer(numEdges, sizeof(uint)*5 + sizeof(float)*6);
        MESHOBS_EDGES_DYNAMIC_BUFFER = new ComputeBuffer(numEdges, sizeof(float)*3);

        MESHOBS_TRANSLATION_FORCES_BUFFER = new ComputeBuffer(numObstacles, sizeof(int)*3);
        MESHOBS_TORQUE_FORCES_BUFFER = new ComputeBuffer(numObstacles, sizeof(int)*3);
        MESHOBS_VELOCITIES_BUFFER = new ComputeBuffer(numObstacles, sizeof(float)*3);

        _obstacles_grid_array = new int[Mathf.Min(_obstacles_grid_array.Length, numGridCells)];
        _rearranged_obstacles_array = new uint[Mathf.Min(_rearranged_obstacles_array.Length, numObstacles)];

        _obstacles_static_array = new OP.ObstacleStatic[Mathf.Min(_obstacles_static_array.Length, numObstacles)];
        _obstacles_dynamic_array = new OP.ObstacleDynamic[Mathf.Min(_obstacles_dynamic_array.Length, numObstacles)];
        _triangles_static_array = new OP.TriangleStatic[Mathf.Min(_triangles_static_array.Length, numTriangles)];
        _triangles_dynamic_array = new OP.TriangleDynamic[Mathf.Min(_triangles_dynamic_array.Length, numTriangles)];
        _vertices_dynamic_array = new OP.VertexDynamic[Mathf.Min(_vertices_dynamic_array.Length, numVertices)];
        _edges_dynamic_array = new float3[Mathf.Min(_edges_dynamic_array.Length, numEdges)];
        
        _translation_forces_array = new int3[Mathf.Min(_translation_forces_array.Length, numObstacles)];

        if (_verbose) Debug.Log("[BUFFER MANAGER] Mesh obs buffers initialized!");
    }

    void Start() {
        if (_PARTICLE_CONTROLLER != null) _PARTICLE_CONTROLLER.Initialize();
        if (_BOIDS_CONTROLLER != null) _BOIDS_CONTROLLER.Initialize();
        if (_OBSTACLES_CONTROLLER != null) _OBSTACLES_CONTROLLER.Initialize();
        if (_PRESSURE_RENDERER != null) _PRESSURE_RENDERER.Initialize();
        if (_CF != null) _CF.Initialize();
        if (_RM != null) _RM.Initialize();
    }

    void LateUpdate() {
        if (_particles_array.Length > 0) PARTICLES_BUFFER.GetData(_particles_array);
        if (_particles_grid_array.Length > 0) PARTICLES_GRID_BUFFER.GetData(_particles_grid_array);
        if (_particles_velocities_array.Length > 0) PARTICLES_VELOCITIES_BUFFER.GetData(_particles_velocities_array);
        if (_particles_densities_array.Length > 0) PARTICLES_DENSITIES_BUFFER.GetData(_particles_densities_array);
        if (_particles_pressures_array.Length > 0) PARTICLES_PRESSURES_BUFFER.GetData(_particles_pressures_array);
        if (_particles_external_forces_array.Length > 0) PARTICLES_EXTERNAL_FORCES_BUFFER.GetData(_particles_external_forces_array);

        if (_obstacles_grid_array.Length > 0) MESHOBS_GRID_BUFFER.GetData(_obstacles_grid_array);
        if (_rearranged_obstacles_array.Length > 0) MESHOBS_REARRANGED_OBSTACLES_BUFFER.GetData(_rearranged_obstacles_array);
        if (_obstacles_static_array.Length > 0) MESHOBS_OBSTACLES_STATIC_BUFFER.GetData(_obstacles_static_array);
        if (_obstacles_dynamic_array.Length > 0) MESHOBS_OBSTACLES_DYNAMIC_BUFFER.GetData(_obstacles_dynamic_array);
        if (_triangles_static_array.Length > 0) MESHOBS_TRIANGLES_STATIC_BUFFER.GetData(_triangles_static_array);
        if (_triangles_dynamic_array.Length > 0) MESHOBS_TRIANGLES_DYNAMIC_BUFFER.GetData(_triangles_dynamic_array);
        if (_vertices_static_array.Length > 0) MESHOBS_VERTICES_STATIC_BUFFER.GetData(_vertices_static_array);
        if (_vertices_dynamic_array.Length > 0) MESHOBS_VERTICES_DYNAMIC_BUFFER.GetData(_vertices_dynamic_array);
        if (_edges_static_array.Length > 0) MESHOBS_EDGES_STATIC_BUFFER.GetData(_edges_static_array);
        if (_edges_dynamic_array.Length > 0) MESHOBS_EDGES_DYNAMIC_BUFFER.GetData(_edges_dynamic_array);
        if (_translation_forces_array.Length > 0) MESHOBS_TRANSLATION_FORCES_BUFFER.GetData(_translation_forces_array);
    }

    public float dt => (_PARTICLE_CONTROLLER!=null) ? _PARTICLE_CONTROLLER.dt : 0f;

    void OnDestroy() {
        if (PARTICLES_BUFFER != null) PARTICLES_BUFFER.Release();
        if (PARTICLES_GRID_BUFFER != null) PARTICLES_GRID_BUFFER.Release();
        if (PARTICLES_VELOCITIES_BUFFER != null) PARTICLES_VELOCITIES_BUFFER.Release();
        if (PARTICLES_DENSITIES_BUFFER != null) PARTICLES_DENSITIES_BUFFER.Release();
        if (PARTICLES_PRESSURES_BUFFER != null) PARTICLES_PRESSURES_BUFFER.Release();
        if (PARTICLES_PRESSURE_GRID_BUFFER != null) PARTICLES_PRESSURE_GRID_BUFFER.Release();
        if (PARTICLES_PRESSURE_FORCES_BUFFER != null) PARTICLES_PRESSURE_FORCES_BUFFER.Release();
        if (PARTICLES_VISCOSITY_FORCES_BUFFER != null) PARTICLES_VISCOSITY_FORCES_BUFFER.Release();
        if (PARTICLES_EXTERNAL_FORCES_BUFFER != null) PARTICLES_EXTERNAL_FORCES_BUFFER.Release();

        if (MESHOBS_GRID_BUFFER != null) MESHOBS_GRID_BUFFER.Release();
        if (MESHOBS_REARRANGED_OBSTACLES_BUFFER != null) MESHOBS_REARRANGED_OBSTACLES_BUFFER.Release();
        if (MESHOBS_OBSTACLES_STATIC_BUFFER != null) MESHOBS_OBSTACLES_STATIC_BUFFER.Release();
        if (MESHOBS_OBSTACLES_DYNAMIC_BUFFER != null) MESHOBS_OBSTACLES_DYNAMIC_BUFFER.Release();
        if (MESHOBS_TRIANGLES_STATIC_BUFFER != null) MESHOBS_TRIANGLES_STATIC_BUFFER.Release();
        if (MESHOBS_TRIANGLES_DYNAMIC_BUFFER != null) MESHOBS_TRIANGLES_DYNAMIC_BUFFER.Release();
        if (MESHOBS_VERTICES_STATIC_BUFFER != null) MESHOBS_VERTICES_STATIC_BUFFER.Release();
        if (MESHOBS_VERTICES_DYNAMIC_BUFFER != null) MESHOBS_VERTICES_DYNAMIC_BUFFER.Release();
        if (MESHOBS_EDGES_STATIC_BUFFER != null) MESHOBS_EDGES_STATIC_BUFFER.Release();
        if (MESHOBS_EDGES_DYNAMIC_BUFFER != null) MESHOBS_EDGES_DYNAMIC_BUFFER.Release();
        if (MESHOBS_TRANSLATION_FORCES_BUFFER != null) MESHOBS_TRANSLATION_FORCES_BUFFER.Release();
        if (MESHOBS_TORQUE_FORCES_BUFFER != null) MESHOBS_TORQUE_FORCES_BUFFER.Release();
        if (MESHOBS_VELOCITIES_BUFFER != null) MESHOBS_VELOCITIES_BUFFER.Release();
    }
}

