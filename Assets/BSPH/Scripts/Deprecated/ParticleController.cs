using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;
using Unity.Mathematics;
using SerializableTypes;
using Helpers;
using WebTools;
using TMPro;
using OP = ObstaclePrimitives.Structs;
using UnityEditor;

public class ParticleController : MonoBehaviour
{

    [Header("== REFERENCES ==")]
    [SerializeField, Tooltip("Reference to the singular BufferManager component that handles all buffers in the system")]
    public BufferManager _BM;
    [SerializeField, Tooltip("Reference to a ParticleGrid component that acts as this controller's grid system")]
    public Grid _GRID = null;
    [SerializeField, Tooltip("Index of the section within the ParticleGrid component that we want to spawn particles inside. If set to -1, then we will use the number defined in this component and not particleGrid's.")]
    private int _SECTION_INDEX = -1;
    [SerializeField, Tooltip("The compute shader that handles all of our GPU calculations. Must-have, otherwise the script will not run")]
    private ComputeShader _SPH_Shader = null;
    [SerializeField, Tooltip("Reference to a Point Cloud Obstacle Manager. Only used if using a PointCloudObstacleManager system somewhere.")]
    private PointCloudObstacleManager _POINT_CLOUD_OBSTACLE_MANAGER = null;
    [SerializeField, Tooltip("Reference to a csvWriter that will record our findings for us.")]
    private CSVWriter recorder;

    [Header("== PARTICLE CONFIGURATIONS ==")]
    [SerializeField, Tooltip("How big (visually) should the particles be?")]
    private float _particleRenderSize = 0.5f;
    public float particleRenderSize => _particleRenderSize;
    public float particleRenderRadius => _particleRenderSize / 2f;
    [SerializeField, ReadOnlyInsp, Tooltip("The max number of particles possible within this controller's grid section")]
    private int _MAX_NUM_PARTICLES = 0;
    [SerializeField, Tooltip("The number of particles that need to be generated. If `sectionIndex` is defined, will default to using that section's # of particles instead. If `sectionIndex` is set to -1, will use this value by default")]
    private int _numParticles = 100;
    public int numParticles => _numParticles;
    [SerializeField, ReadOnlyInsp]
    private int[] _numParticlesPerAxis;
    [SerializeField, ReadOnlyInsp]
    private int _numParticlesPerGridCell;
    [SerializeField, ReadOnlyInsp]
    private float[] _spawnBounds;
    [SerializeField, Tooltip("The mass of each particle. We generally assume all particles have the same particle mass.")]
    private float _particleMass = 1f;


    [Header("== SPH CONFIGURATIONS ==")]
    [SerializeField, Tooltip("The delta time of the simulation. If set to any value < 0, then the simulation will default to using `Time.deltaTime`.")]
    private float _dt = 0.00825f;
    public float dt => _dt;
    [SerializeField, Tooltip("What's the gravitational force exerted on all particles?")]
    private float[] _g = {0f, -9.81f, 0f};
    public float[] g => _g;
    [SerializeField, Tooltip("`h` - the smoothing kernel radius used all across the SPH simulation")]
    private float _h = 1.5f;
    public float h => _h;
    [SerializeField, Tooltip("The `p_0` part of the SPH equations - the resting density the fluid")]
    private float _rest_density = 1000f;
    [SerializeField, Tooltip("`mu` - the viscosity coefficient of the fluid")]
    private float _mu = 0.5f;
    [SerializeField, Tooltip("The damping effect whenever a particle hits the bounds of the simulation")]
    private float _damping_effect = -0.5f;
    [SerializeField, Tooltip("`k` - technically the ideal gas constant, but some people (Bindel) refer this as the bulk modulus value.")]
    private float _k = 1f;
    public float k => _k;
    [SerializeField, Tooltip("`nk` - similar to `k` above, but for really close particle-particle interactions")]
    private float _nk = 1f;
    public float nk => _nk;
    [SerializeField, Tooltip("`c` - the speed of sound (10^-4 m/s, to model low Reynolds number flows). Used as a part of the `equation of state` to relate density to pressure. This is derived from the method mentioned in [https://www.sciencedirect.com/science/article/pii/S0045793021003145#b19], but the speed of sound value is derived from Morris et al. [https://www.sciencedirect.com/science/article/pii/S0021999197957764]")]
    private float _c = 5.77f;
    public float c => _c;
    private float _dt_passed = 0f;
    public float dt_passed => _dt_passed;
    private float _real_time_elapsed = 0f;
    private int _frames_elapsed = 0;
    private float _total_dt_elapsed = 0f;
    private float _total_real_time_elapsed = 0f;
    private int _total_frames_elapsed = 0;
    private float _fps = 0f;

    [SerializeField] private float _spawnDistanceBetweenParticles = 1.45f;
    public float spawnDistanceBetweenParticles => _spawnDistanceBetweenParticles;

    int size_property = Shader.PropertyToID("size");
    int num_particles_per_cell_property = Shader.PropertyToID("numParticlesPerCell");
    int render_touching_property = Shader.PropertyToID("renderTouching");
    int particle_buffer_property = Shader.PropertyToID("particle_buffer");
    int density_buffer_property = Shader.PropertyToID("density_buffer");
    int grid_cell_buffer_property = Shader.PropertyToID("grid_cell_buffer");
    int render_limits_property = Shader.PropertyToID("render_limits_buffer");

    int rest_density_property = Shader.PropertyToID("rest_density");
    int bulk_modulus_property = Shader.PropertyToID("bulk_modulus");
    int max_color_val_property = Shader.PropertyToID("max_color_val");
    int render_color_property = Shader.PropertyToID("render_color");

    public enum RenderType { Off, Particles, ParticlesTouchingObstacles, GridCells, Both }
    [Header("== DEBUG CONTROLS ==")]
    [SerializeField] private float _startDelay = 2f;
    public float startDelay => _startDelay;
    private float _timeStarted = 0f, _timePassed = 0f;
    [SerializeField] private RenderType _renderType = RenderType.Particles;
    [SerializeField] private float[] _renderLimits;

    private float t;

    [SerializeField] private bool render_colors = false;
    [SerializeField] private float max_color_value = 0.0f;
    [SerializeField] private bool gizmos_particles = false;
    [SerializeField] private Color gizmos_particle_color = Color.blue;
    [SerializeField] private bool gizmos_velocities = false;
    [SerializeField] private Color gizmos_velocity_color = Color.yellow;
    [SerializeField] private bool gizmos_accelerations = false;
    [SerializeField] private Color gizmos_accelerations_color = Color.white;
    [SerializeField] private bool gizmos_pressure_force = false;
    [SerializeField] private Color gizmos_pressure_force_color = Color.red;
    [SerializeField] private bool gizmos_viscosity_force = false;
    [SerializeField] private Color gizmos_viscosity_force_color = Color.green;
    [SerializeField] private bool gizmos_external_force = false;
    [SerializeField] private Color gizmos_external_force_color = Color.yellow;
    
    public enum RecordSettings { Off, CSV }
    public class Recording {
        public OP.Particle[] temp_particles;        
        public float3[] temp_velocities;
        public float[] temp_densities;
        public float3[] temp_pressures;
        public RecordingFrame frameData;
        public Recording(
                OP.Particle[] pa, float3[] ve, float[] de, float3[] pr, 
                int frame, float dt_passed_sec, float rl_passed_sec, float fps, 
                int n_obstacles, int n_vertices, int n_edges, int n_triangles
        ) {
            this.temp_particles = pa;
            this.temp_velocities = ve;
            this.temp_densities = de;
            this.temp_pressures = pr;
            this.frameData = new RecordingFrame(
                frame, dt_passed_sec, rl_passed_sec, fps,
                n_obstacles, n_vertices, n_edges, n_triangles 
            );
        }
    }
    public class RecordingFrame {
        public int frame;
        public float dt_passed_sec;
        public float rl_passed_sec;
        public float fps;
        public int n_obstacles, n_vertices, n_edges, n_triangles;
        public string fileName => $"nFrames_{frame}-dt_{dt_passed_sec}-realtime_{rl_passed_sec}-fps_{fps}";
        public RecordingFrame(
                int frame, float dt_passed_sec, float rl_passed_sec, float fps,
                int n_obstacles, int n_vertices, int n_edges, int n_triangles
        ) {
            this.frame = frame;
            this.dt_passed_sec = dt_passed_sec;
            this.rl_passed_sec = rl_passed_sec;
            this.fps = fps;
            this.n_obstacles = n_obstacles;
            this.n_vertices = n_vertices;
            this.n_edges = n_edges;
            this.n_triangles = n_triangles;
        }
    }

    [Header("== RECORDING CONFIGURATIONS ==")]
    [SerializeField, Tooltip("Public toggle to determine if we should be recording.")]
    private RecordSettings _record_statistics = RecordSettings.Off;
    [SerializeField, Tooltip("The duraction (in seconds) of the recording session")] 
    private float _record_interval = 1f;
    private IEnumerator _recordCoroutine = null;
    private Queue<Recording> _recordQueue = new Queue<Recording>();
    private Queue<RecordingFrame> _recordFrameQueue = new Queue<RecordingFrame>();
    [SerializeField] private RecordingCanvas textboxes = null;

    void OnDrawGizmos() {
        if (!Application.isPlaying) return;
        int minLength = Mathf.Min(_BM.particles_array.Length, _BM.particles_velocities_array.Length);
        minLength = Mathf.Min(minLength, _BM.particles_external_forces_array.Length);
        if (minLength == 0) return;
        for(int i = 0; i < minLength; i++) {
            Vector3 pos = new Vector3(_BM.particles_array[i].position[0], _BM.particles_array[i].position[1], _BM.particles_array[i].position[2]);
            if (gizmos_particles) {
                Gizmos.color = gizmos_particle_color;
                Gizmos.DrawSphere(
                    _BM.particles_array[i].position,
                    _h
                );
            }
            if (gizmos_velocities) {
                Handles.color = gizmos_velocity_color;
                Handles.DrawLine(
                    _BM.particles_array[i].position, 
                    _BM.particles_array[i].position + _BM.particles_velocities_array[i] * 10,
                    3
                );
            }
            if (gizmos_external_force) {
                Handles.color = gizmos_external_force_color;
                Handles.DrawLine(
                    _BM.particles_array[i].position, 
                    _BM.particles_array[i].position + _BM.particles_external_forces_array[i].external_force * 10,
                    3
                );
            }
        }
        /*
        if (!Application.isPlaying) return;
        // For now, just render particles
        OP.Particle[] temp_particles = new OP.Particle[_numParticles];
        _BM.PARTICLES_BUFFER.GetData(temp_particles);
        float3[] temp_velocities = new float3[_numParticles];
        _BM.PARTICLES_VELOCITIES_BUFFER.GetData(temp_velocities);
        float3[] temp_accelerations = new float3[_numParticles];
        FORCES_BUFFER.GetData(temp_accelerations);
        float3[] temp_pressure_forces = new float3[_numParticles];
        _BM.PARTICLES_PRESSURE_FORCES_BUFFER.GetData(temp_pressure_forces);
        float3[] temp_viscosity_forces = new float3[_numParticles];
        _BM.PARTICLES_VISCOSITY_FORCES_BUFFER.GetData(temp_viscosity_forces);
        float3[] temp_external_forces = new float3[_numParticles];
        //EXTERNAL_FORCES_BUFFER.GetData(temp_external_forces);

        for(int i = 0; i < _numParticles; i++) {
            Vector3 pos = new Vector3(temp_particles[i].position[0], temp_particles[i].position[1], temp_particles[i].position[2]);
            if (gizmos_particles) {
                Gizmos.color = gizmos_particle_color;
                Gizmos.DrawSphere(pos, particleRenderRadius);
            }
            if (gizmos_velocities) {
                Gizmos.color = gizmos_velocity_color;
                Gizmos.DrawRay(pos, new Vector3(temp_velocities[i][0], temp_velocities[i][1], temp_velocities[i][2]));
            }
            if (gizmos_accelerations) {
                Gizmos.color = gizmos_accelerations_color;
                Gizmos.DrawRay(pos, new Vector3(temp_accelerations[i][0], temp_accelerations[i][1], temp_accelerations[i][2]));
            }
            if (gizmos_pressure_force) {
                Gizmos.color = gizmos_pressure_force_color;
                Gizmos.DrawRay(pos, new Vector3(temp_pressure_forces[i][0], temp_pressure_forces[i][1], temp_pressure_forces[i][2]));
            }
            if (gizmos_viscosity_force) {
                Gizmos.color = gizmos_viscosity_force_color;
                Gizmos.DrawRay(pos, new Vector3(temp_viscosity_forces[i][0], temp_viscosity_forces[i][1], temp_viscosity_forces[i][2]));
            }
        }
        */
    }

    public void CalculateParticles() {
        //_spawnDistanceBetweenParticles = Mathf.Clamp(_spawnDistanceBetweenParticles,particleRenderRadius,_h);
        // GridCellSize is a global value across this grid. _particleRenderSize is unique to this controller
        int numParticlesPerCellAxis = Mathf.CeilToInt(_GRID.gridCellSize / _spawnDistanceBetweenParticles);
        // If we have a section index, we use that section's bounds. Otherwise, we use the grid's inner bounds
        Vector3 bs = (_SECTION_INDEX >= 0) ? _GRID.sections[_SECTION_INDEX].dimensionsV3 : _GRID.innerBoundsV3;
        // We store number of particles per axis; particle render size is unique to the controller
        _numParticlesPerAxis = new int[3] {
            Mathf.FloorToInt(bs[0] / _spawnDistanceBetweenParticles),
            Mathf.FloorToInt(bs[1] / _spawnDistanceBetweenParticles),
            Mathf.FloorToInt(bs[2] / _spawnDistanceBetweenParticles)
        };
        // The number of particles that would fit inside of the cell is defined already as a get function. We just need to apply it
        _MAX_NUM_PARTICLES = _numParticlesPerAxis[0] * _numParticlesPerAxis[1] * _numParticlesPerAxis[2];
        _numParticles = Mathf.Clamp(_numParticles,0,_MAX_NUM_PARTICLES);

        // Now we calculate based on grid dimensions!
        // The calculation for the # of particles is based on how many grid cells are below the line formed by _WATER_LEVEL_TRANSFORM's y-axis position
        _numParticlesPerGridCell = numParticlesPerCellAxis * _GRID.numCellsPerAxis[0] 
            + numParticlesPerCellAxis * _GRID.numCellsPerAxis[1]
            + numParticlesPerCellAxis * _GRID.numCellsPerAxis[2];
        _spawnBounds = (_SECTION_INDEX >= 0) ? _GRID.sections[_SECTION_INDEX].bounds : _GRID.innerBounds;
    }

    public void Initialize() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SPH_Shader == null) {
            Debug.LogError("SPH - ERROR: Cannot operate if either `GRID` or `SPH_SHADER` is set to `null`. Please define these references and restart the simulation.");
            return;
        }

        // Initialize Key Variables, many of which are from `grid`
        InitializeVariables();

        // We pass along key variable values to our shader too
        InitializeShaderVariables();

        // We initialize the kernels and the buffers used for the GPU
        InitializeKernels();
        InitializeBuffers();

        // We dispatch the initial functions for setup
        _SPH_Shader.Dispatch(_CLEAR_GRID, _NUM_BLOCKS_GRID,1,1);
        _SPH_Shader.Dispatch(_GENERATE_PARTICLES, _NUM_BLOCKS_PARTICLES,1,1);

        // If we are recording, we initialize a recording session
        if (_record_statistics != RecordSettings.Off) PrepareRecording();

        Debug.Log("ParticleController: Particles Initialized!");
        t = Time.time;
        _timeStarted = Time.time;

        if (textboxes != null) {
            textboxes.nParticles.gameObject.SetActive(true);
            textboxes.nObstacles.gameObject.SetActive(true);
            textboxes.obstacleDetails.gameObject.SetActive(true);
            textboxes.fps.gameObject.SetActive(true);
            textboxes.deltaTime.gameObject.SetActive(true);
            textboxes.realTime.gameObject.SetActive(true);
        }
    }

    private void PrepareRecording() {
        // Prepare the coroutine for the recordings
        _recordCoroutine = RecordCycle();
        StartCoroutine(_recordCoroutine);
        // Make our initial recording
        MakeRecord();
    }

    private int _BLOCK_SIZE = 512;
    private int _NUM_BLOCKS_GRID;
    private int _NUM_BLOCKS_PARTICLES;
    [SerializeField, ReadOnlyInsp] private int _NUM_BOUNDARY_PARTICLES;
    private void InitializeVariables() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SPH_Shader == null) return;

        // Determine block number for GPU threading
        _NUM_BLOCKS_GRID = Mathf.CeilToInt((float)_GRID.numGridCells / (float)_BLOCK_SIZE);
        _NUM_BLOCKS_PARTICLES = Mathf.CeilToInt((float)_numParticles / (float)_BLOCK_SIZE);

        // Determine the number of boundary particles
        _NUM_BOUNDARY_PARTICLES = (_POINT_CLOUD_OBSTACLE_MANAGER != null) ? _POINT_CLOUD_OBSTACLE_MANAGER.numBoundaryParticles : 0;

        // Debugs to show our current settings
        Debug.Log($"Number of particle grid cells per axis: ({_GRID.numCellsPerAxis[0]}, {_GRID.numCellsPerAxis[1]}, {_GRID.numCellsPerAxis[2]})");
        Debug.Log($"Total # of particle grid cells: {_GRID.numGridCells}");
        Debug.Log($"# Particles per grid cell: {_numParticlesPerGridCell}");
        Debug.Log($"Size of particle neighbors: {_GRID.numGridCells * _numParticlesPerGridCell}");
    }

    private void InitializeShaderVariables() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SPH_Shader == null) return;

        // == WORLD CONFIGURATIONS ==
        _SPH_Shader.SetFloat("gridCellSize", _GRID.gridCellSize);
        _SPH_Shader.SetInts("numCellsPerAxis", _GRID.numCellsPerAxis);
        _SPH_Shader.SetInt("total_number_of_cells", _GRID.numGridCells);
        //_SPH_Shader.SetFloats("origin", _GRID.originF);
        _SPH_Shader.SetFloat("epsilon", Mathf.Epsilon);
        _SPH_Shader.SetFloat("pi", Mathf.PI);
        _SPH_Shader.SetFloat("spawnDistanceBetweenParticles",_spawnDistanceBetweenParticles);
        /*
        float[] gridScale = new float[3]{
            _GRID.originF[0] - (_GRID.numCellsPerAxisF[0] * _GRID.gridCellSize)/2f,
            _GRID.originF[1] - (_GRID.numCellsPerAxisF[1] * _GRID.gridCellSize)/2f,
            _GRID.originF[2] - (_GRID.numCellsPerAxisF[2] * _GRID.gridCellSize)/2f
        };
        _SPH_Shader.SetFloats("gridScaling", gridScale);
        */
        _SPH_Shader.SetFloat("gridScalingX", _GRID.gridScaling[0]);
        _SPH_Shader.SetFloat("gridScalingY", _GRID.gridScaling[1]);
        _SPH_Shader.SetFloat("gridScalingZ", _GRID.gridScaling[2]);

        // == PARTICLE CONFIGURATIONS ==
        _SPH_Shader.SetInt("numParticles", _numParticles);
        _SPH_Shader.SetInts("numParticlesPerAxis", _numParticlesPerAxis);
        _SPH_Shader.SetInt("numParticlesPerGridCell", _numParticlesPerGridCell);
        _SPH_Shader.SetInt("numBoundaryParticles", _NUM_BOUNDARY_PARTICLES);

        // == GPU SETTINGS
        _SPH_Shader.SetInt("numBlocks", _NUM_BLOCKS_GRID);
        //_SPH_Shader.SetInt("randomSeed", Random.Range(0,int.MaxValue));

        // Update variables that may change over time due to modifying inspector values
        UpdateShaderVariables();
    }

    private void UpdateShaderVariables(bool updateDT = false) {
        // == PARTICLE CONFIGURATIONS ==
        _SPH_Shader.SetFloat("particleRenderRadius", particleRenderRadius);

        // == FLUID MECHANICS ==
        _SPH_Shader.SetFloat("smoothingRadius", _h);
        _SPH_Shader.SetFloat("radius2", _h * _h);
        _SPH_Shader.SetFloat("radius3", Mathf.Pow(_h,3));
        _SPH_Shader.SetFloat("radius5", Mathf.Pow(_h,5));
        _SPH_Shader.SetFloat("radius6", Mathf.Pow(_h,6));
        _SPH_Shader.SetFloat("radius8",  Mathf.Pow(_h,8));
        _SPH_Shader.SetFloat("radius9", Mathf.Pow(_h,9));

        _SPH_Shader.SetFloat("particleMass", _particleMass);
        _SPH_Shader.SetFloat("rest_density", _rest_density);
        _SPH_Shader.SetFloat("viscosity_coefficient", _mu);
        _SPH_Shader.SetFloat("damping", _damping_effect);
        _SPH_Shader.SetFloat("bulkModulus", _k);
        _SPH_Shader.SetFloat("nearBulkModulus", _nk);
        _SPH_Shader.SetFloats("g", _g);
        _SPH_Shader.SetFloat("c", _c);

        _renderLimits[0] = Mathf.Max(_renderLimits[0],_GRID.outerBounds[0]);
        _renderLimits[1] = Mathf.Max(_renderLimits[1],_GRID.outerBounds[1]);
        _renderLimits[2] = Mathf.Max(_renderLimits[2],_GRID.outerBounds[2]);
        _renderLimits[3] = Mathf.Min(_renderLimits[3],_GRID.outerBounds[3]);
        _renderLimits[4] = Mathf.Min(_renderLimits[4],_GRID.outerBounds[4]);
        _renderLimits[5] = Mathf.Min(_renderLimits[5],_GRID.outerBounds[5]);

        if (updateDT) {
            float t = (_dt >= 0) ? _dt : Time.deltaTime;
            _SPH_Shader.SetFloat("dt", t);
        }
    }

    private int _CLEAR_GRID;
    private int _GENERATE_PARTICLES;
    private int _UPDATE_GRID;
    private int _COMPUTE_DENSITY;
    private int _COMPUTE_PRESSURE_FORCES, _COMPUTE_VISCOSITY_FORCES;
    private int _INTEGRATE;
    //_INTEGRATE_DEBUG;
    private int _DAMPEN_BY_BOUNDS;
    private void InitializeKernels() {
        // Used on initialization. _CLEAR_GRID also is performed at the beginning of each update loop
        _CLEAR_GRID = _SPH_Shader.FindKernel("ClearGrid");
        _GENERATE_PARTICLES = _SPH_Shader.FindKernel("GenerateParticles");
        // These are run during each update loop
        _UPDATE_GRID = _SPH_Shader.FindKernel("UpdateGridCellCounts");
        _COMPUTE_DENSITY = _SPH_Shader.FindKernel("CV_ComputeDensity");
        _COMPUTE_PRESSURE_FORCES = _SPH_Shader.FindKernel("ComputePressureForces");
        _COMPUTE_VISCOSITY_FORCES = _SPH_Shader.FindKernel("ComputeViscosityForces");
        _INTEGRATE = _SPH_Shader.FindKernel("Integrate");
        //_INTEGRATE_DEBUG = _SPH_Shader.FindKernel("Integrate_Debug");
        _DAMPEN_BY_BOUNDS = _SPH_Shader.FindKernel("DampenByBounds");
    }

    public ComputeBuffer ARG_BUFFER;
    public ComputeBuffer CELL_LIMITS_BUFFER;
    public ComputeBuffer BOUNDS_BUFFER, SPAWN_BOUNDS_BUFFER;
    public ComputeBuffer PARTICLE_NEIGHBORS_BUFFER, PARTICLE_OFFSETS_BUFFER;
    public ComputeBuffer FORCES_BUFFER;
    public ComputeBuffer BOUNDARY_PARTICLES_BUFFER;

    private ComputeBuffer RENDER_LIMITS_BUFFER;

    private void InitializeBuffers() {
        uint[] arg = {_GRID.particle_mesh.GetIndexCount(0), (uint)(_numParticles), _GRID.particle_mesh.GetIndexStart(0), _GRID.particle_mesh.GetBaseVertex(0), 0};
        ARG_BUFFER = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        ARG_BUFFER.SetData(arg);

        // We call upon our BufferManager to initialize all buffers related to particles. Seems appropriate, given that the BufferManager is holding onto all the buffers for us.
        _BM.InitializeParticleBuffers(_numParticles + _NUM_BOUNDARY_PARTICLES, _GRID.numGridCells);
        //_BM.PARTICLES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float)*6+sizeof(int));
        //_BM.PARTICLES_GRID_BUFFER = new ComputeBuffer(_GRID.numGridCells, sizeof(int)*2 + sizeof(float)*3);
        //_BM.PARTICLES_DENSITIES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float));
        //_BM.PARTICLES_PRESSURE_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float));
        //_BM.PARTICLES_VELOCITIES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, 3 * sizeof(float));
        //_BM.PARTICLES_PRESSURE_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
        //_BM.PARTICLES_VISCOSITY_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
        //_BM.PARTICLES_EXTERNAL_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(uint) + sizeof(int)*8 + sizeof(float)*27);

        // We also have to initialize some custom compute buffers purely for our own purposes
        BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
        BOUNDS_BUFFER.SetData(_GRID.outerBounds);

        RENDER_LIMITS_BUFFER = new ComputeBuffer(6, sizeof(float));
        RENDER_LIMITS_BUFFER.SetData(_renderLimits);

        CELL_LIMITS_BUFFER = new ComputeBuffer(_GRID.numGridCells, sizeof(int)*9 + sizeof(float)*3);
        OP.CellLimits[] cellLimits = new OP.CellLimits[_GRID.numGridCells];
        for(int x = 0; x < _GRID.numCellsPerAxis[0]; x++) {
            for(int y = 0; y < _GRID.numCellsPerAxis[1]; y++) {
                for(int z = 0; z < _GRID.numCellsPerAxis[2]; z++) {
                    int i = GetProjectedGridIndexFromXYZ(x,y,z,_GRID.numCellsPerAxis);
                    cellLimits[i] = new OP.CellLimits();
                    cellLimits[i].id = new(x,y,z);
                    cellLimits[i].position = new(
                        _GRID.outerBounds[0] + (_GRID.gridCellSize * 0.5f) + (_GRID.gridCellSize * x),
                        _GRID.outerBounds[1] + (_GRID.gridCellSize * 0.5f) + (_GRID.gridCellSize * y),
                        _GRID.outerBounds[2] + (_GRID.gridCellSize * 0.5f) + (_GRID.gridCellSize * z)
                    );
                    cellLimits[i].lowerLimits = new(
                        Mathf.Max(0, x-1),
                        Mathf.Max(0, y-1),
                        Mathf.Max(0, z-1)
                    );
                    cellLimits[i].upperLimits = new(
                        Mathf.Min(x+1, _GRID.numCellsPerAxis[0]-1),
                        Mathf.Min(y+1, _GRID.numCellsPerAxis[1]-1),
                        Mathf.Min(z+1, _GRID.numCellsPerAxis[2]-1)
                    );
                }
            }
        }
        CELL_LIMITS_BUFFER.SetData(cellLimits);

        SPAWN_BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
        SPAWN_BOUNDS_BUFFER.SetData(_spawnBounds);

        if (_POINT_CLOUD_OBSTACLE_MANAGER != null && _NUM_BOUNDARY_PARTICLES > 0) {
            Debug.Log($"{_POINT_CLOUD_OBSTACLE_MANAGER.boundaryParticles.Count}");
            _BM.PARTICLES_BUFFER.SetData(_POINT_CLOUD_OBSTACLE_MANAGER.boundaryParticles.ToArray(), 0, _numParticles, _NUM_BOUNDARY_PARTICLES);
        }
        Debug.Log($"Number of Particles: {_numParticles + _NUM_BOUNDARY_PARTICLES}");
        PARTICLE_OFFSETS_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(int));
        PARTICLE_NEIGHBORS_BUFFER = new ComputeBuffer(_GRID.numGridCells * _numParticlesPerGridCell, sizeof(int));

        // We attempt to... "adjust" the density of the boundary particle ssuch tath they are able to really have an effect on fluid particles that hit them
        /*
        if (_POINT_CLOUD_OBSTACLE_MANAGER != null && _NUM_BOUNDARY_PARTICLES > 0) {
            float[] boundaryDensities = new float[_NUM_BOUNDARY_PARTICLES];
            for(int b = 0; b < _NUM_BOUNDARY_PARTICLES; b++) {
                boundaryDensities[b] = 100000000f;
            }
            _BM.PARTICLES_DENSITIES_BUFFER.SetData(boundaryDensities, 0, _numParticles, _NUM_BOUNDARY_PARTICLES);
        }
        */

        FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);

        // Setting the buffers

        _SPH_Shader.SetBuffer(_CLEAR_GRID, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_CLEAR_GRID, "bounds", BOUNDS_BUFFER);

        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "bounds", SPAWN_BOUNDS_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "densities", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "pressures", _BM.PARTICLES_PRESSURE_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "externalForces", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);

        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particleOffsets", PARTICLE_OFFSETS_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);

        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "cellLimits", CELL_LIMITS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "densities", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "near_densities", _BM.PARTICLES_NEAR_DENSITIES_BUFFER);
        //_SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        //_SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "pressure", _BM.PARTICLES_PRESSURE_BUFFER);

        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "cellLimits", CELL_LIMITS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "densities", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "near_densities", _BM.PARTICLES_NEAR_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_PRESSURE_FORCES, "pressures", _BM.PARTICLES_PRESSURE_BUFFER);

        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "cellLimits", CELL_LIMITS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "densities", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_VISCOSITY_FORCES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);

        _SPH_Shader.SetBuffer(_INTEGRATE, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "densities", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "externalForces", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);

        /*
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        */

        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "bounds", BOUNDS_BUFFER);
        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
    }

    // Update is called once per frame
    void Update() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SPH_Shader == null) return;
        //if (_record_statistics != RecordSettings.Off && !_recording_verified) return;
        // Update shader variables!
        UpdateShaderVariables(true);
        // Check if we can integrate particles, based on _startDelay and _timePassed
        if ((Time.time - _timeStarted) < _startDelay) return;
        // Update Particles!
        UpdateMullerParticles();
    }
    
    // =====================================
    // ===== MULLER'S COMPRESSIBLE SPH =====
    // =====================================

    public void UpdateMullerParticles() {
        //DebugBufferParticle("POSITIONS",1000,PARTICLES_BUFFER);
        // Reset the grid and num neighbors buffer, if we are debugging
        _SPH_Shader.Dispatch(
            _CLEAR_GRID, 
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[0] / 8f),
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[1] / 8f),
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[2] / 8f)
        );
        
        // Calculate how many particles are in each grid cell, and get each particle's offset for their particular grid cell
        _SPH_Shader.Dispatch(_UPDATE_GRID, Mathf.CeilToInt((float)(_numParticles+_NUM_BOUNDARY_PARTICLES) / 256f),1,1);
        
        // We perform the updates for density, pressure, and force/acceleration. 
        _SPH_Shader.Dispatch(_COMPUTE_DENSITY, Mathf.CeilToInt((float)(_numParticles+_NUM_BOUNDARY_PARTICLES) / 256f), 1, 1);
        _SPH_Shader.Dispatch(_COMPUTE_PRESSURE_FORCES, Mathf.CeilToInt((float)(_numParticles+_NUM_BOUNDARY_PARTICLES) / 256f),1,1);
        _SPH_Shader.Dispatch(_COMPUTE_VISCOSITY_FORCES, Mathf.CeilToInt((float)(_numParticles+_NUM_BOUNDARY_PARTICLES) / 256f),1,1);

        //_SPH_Shader.Dispatch(compute_external_acceleration_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);

        // Integrate over particles, update their positions after taking all force calcualtions into account
        // Note: We ONLY update the positions and velocities of the active fluid particles, not any boundary particles we might encounter.
        _SPH_Shader.Dispatch(_INTEGRATE, Mathf.CeilToInt((float)_numParticles / 256f), 1, 1);
        // Make sure that particles are within bounds, limit them if so
        //_SPH_Shader.Dispatch(_DAMPEN_BY_BOUNDS, Mathf.CeilToInt((float)_numParticles / 256f), 1, 1);

        // Update the time elapsed and frames elapsed
        _dt_passed += _dt;
        _real_time_elapsed += Time.deltaTime;
        _frames_elapsed += 1;
        _fps = (_real_time_elapsed > 0f) ? (float)_frames_elapsed / _real_time_elapsed : 0f;     // Calculate the FPS based on frames and time passed

        // Update visualization textboxes if they exist
        if (textboxes != null) {
            textboxes.nParticles.text = $"{_numParticles} Particles";
            textboxes.nObstacles.text = $"{_BM.MESHOBS_OBSTACLES_STATIC_BUFFER.count} Obstacles";
            textboxes.obstacleDetails.text = $"(V:{_BM.MESHOBS_VERTICES_STATIC_BUFFER.count} | E:{_BM.MESHOBS_EDGES_STATIC_BUFFER.count} | T:{_BM.MESHOBS_TRIANGLES_STATIC_BUFFER.count})";
            textboxes.fps.text = $"{_fps} FPS";
            textboxes.deltaTime.text = $"Simulation time: {_total_dt_elapsed + _dt_passed} sec.";
            textboxes.realTime.text = $"Real time: {_total_real_time_elapsed + _real_time_elapsed} sec.";
        }

        // If we're recording, record our session
        // Also note: if we're waiting for the recording to start, we won't actually record anything yet.
        if (_real_time_elapsed >= _record_interval) {
            // Here, we want to capture some important details
            // - We will instantiate a new file for every recording we capture. 
            // - We will time each capture on a time interval, meaning we will at least ensure that we won't lag the system.
            // - Each file will contain the following:
            //      1. the filename should contain the recorded timestamp (total elapsed time), the frame number itself, and the fps
            //      2. all subsequent lines should represent a particle each
            //      3. For each line/particle, we'll record their ID, position, velocity, density, and pressure value

            // If our time elapsed goes beyond our record_duration, we have to record
            MakeRecord();
        }

        /*
        int rest_density_property = Shader.PropertyToID("rest_density");
        int bulk_modulus_property = Shader.PropertyToID("bulk_modulus");
        int max_color_val_property = Shader.PropertyToID("max_color_val");
        int render_color_property = Shader.PropertyToID("render_color");
        */

        RENDER_LIMITS_BUFFER.SetData(_renderLimits);

        // In our materials buffer, let's set the necessary buffers for rendering particles and the like
        _GRID.particle_material.SetBuffer(particle_buffer_property, _BM.PARTICLES_BUFFER);
        _GRID.particle_material.SetBuffer(density_buffer_property, _BM.PARTICLES_DENSITIES_BUFFER);
        _GRID.particle_material.SetBuffer(render_limits_property, RENDER_LIMITS_BUFFER);
        // Same for variables
        _GRID.particle_material.SetFloat(size_property, particleRenderSize);
        //_GRID.particle_material.SetFloat(num_particles_per_cell_property, (float)_numParticlesPerGridCell);
        _GRID.particle_material.SetInt(render_color_property, render_colors ? 1 : 0);
        _GRID.particle_material.SetFloat(rest_density_property, _rest_density);
        _GRID.particle_material.SetFloat(bulk_modulus_property, _k);
        _GRID.particle_material.SetFloat(max_color_val_property, max_color_value);

        switch(_renderType) {
            case RenderType.Particles:
                Graphics.DrawMeshInstancedIndirect(
                    _GRID.particle_mesh, 
                    0, 
                    _GRID.particle_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
                break;
            /*
            case RenderType.ParticlesTouchingObstacles:
                _GRID.particle_material.SetFloat(num_particles_per_cell_property, (float)_numParticlesPerGridCell);
                _GRID.particle_material.SetFloat(render_denom_property, _renderDenom);
                _GRID.particle_material.SetInt(render_touching_property, 1);
                _GRID.particle_material.SetBuffer(render_limits_property, RENDER_LIMITS_BUFFER);
                Graphics.DrawMeshInstancedIndirect(
                    _GRID.particle_mesh, 
                    0, 
                    _GRID.particle_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
                break;
            case RenderType.GridCells:
                _GRID.grid_cell_material.SetFloat(size_property, _GRID.gridCellSize);
                _GRID.grid_cell_material.SetBuffer(render_limits_property, RENDER_LIMITS_BUFFER);
                _GRID.grid_cell_material.SetBuffer(grid_cell_buffer_property, _BM.PARTICLES_GRID_BUFFER);
                Graphics.DrawMeshInstancedIndirect(
                    _GRID.grid_cell_mesh, 
                    0, 
                    _GRID.grid_cell_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
                break;
            case RenderType.Both:
                _GRID.particle_material.SetFloat(size_property, particleRenderSize);
                _GRID.particle_material.SetFloat(num_particles_per_cell_property, (float)_numParticlesPerGridCell);
                Graphics.DrawMeshInstancedIndirect(
                    _GRID.particle_mesh, 
                    0, 
                    _GRID.particle_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
                _GRID.grid_cell_material.SetFloat(size_property, _GRID.gridCellSize);
                _GRID.grid_cell_material.SetBuffer(render_limits_property, RENDER_LIMITS_BUFFER);
                _GRID.grid_cell_material.SetBuffer(grid_cell_buffer_property, _BM.PARTICLES_GRID_BUFFER);
                Graphics.DrawMeshInstancedIndirect(
                    _GRID.grid_cell_mesh, 
                    0, 
                    _GRID.grid_cell_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
                break;
            */
            default:
                return;
        }

    }

    private void MakeRecord() {
        // STATS
        _total_dt_elapsed += _dt_passed;                   // Update the total time and frames elapsed
        _total_real_time_elapsed += _real_time_elapsed;
        _total_frames_elapsed += _frames_elapsed;
        _dt_passed = 0f;                                     // Reset _frames and _dt_passed
        _real_time_elapsed = 0f;
        _frames_elapsed = 0;

        if (_record_statistics != RecordSettings.CSV) return;

        // GET DATA
        OP.Particle[] temp_particles = new OP.Particle[_numParticles];
            _BM.PARTICLES_BUFFER.GetData(temp_particles);
        float3[] temp_velocities = new float3[_numParticles];
            _BM.PARTICLES_VELOCITIES_BUFFER.GetData(temp_velocities);
        float[] temp_densities = new float[_numParticles];
            _BM.PARTICLES_DENSITIES_BUFFER.GetData(temp_densities);
        float3[] temp_pressures = new float3[_numParticles];
            _BM.PARTICLES_PRESSURE_FORCES_BUFFER.GetData(temp_pressures);

        // Remember that each of these dynamic and static variations of the vertices, triangles, and edges all are the same size. So it doesn't matter if we get the number of them from either the static or dyanmic buffers
        int n_obstacles = (_BM.MESHOBS_OBSTACLES_STATIC_BUFFER != null) ? _BM.MESHOBS_OBSTACLES_STATIC_BUFFER.count : 0;
        int n_vertices = (_BM.MESHOBS_VERTICES_STATIC_BUFFER != null) ? _BM.MESHOBS_VERTICES_STATIC_BUFFER.count : 0;
        int n_edges = (_BM.MESHOBS_EDGES_STATIC_BUFFER != null) ? _BM.MESHOBS_EDGES_STATIC_BUFFER.count : 0;
        int n_triangles = (_BM.MESHOBS_TRIANGLES_STATIC_BUFFER != null) ? _BM.MESHOBS_TRIANGLES_STATIC_BUFFER.count : 0;
        
        // Create a new entry in a queue that holds all the recordings    
        _recordQueue.Enqueue(new Recording(
            temp_particles, temp_velocities, temp_densities, temp_pressures, 
            _total_frames_elapsed, _total_dt_elapsed, _total_real_time_elapsed, _fps,
            n_obstacles, n_vertices, n_edges, n_triangles
        ));
    }

    private IEnumerator RecordCycle() {
        while(true) {
            // Don't do anything if our queue is empty
            if ( _recordQueue.Count == 0) {
                yield return null;
                continue;
            }

            // Get the top record in our queue
            Recording r = _recordQueue.Dequeue();

            // Initialize the recorder
            recorder.fileName = r.frameData.fileName;
            recorder.Initialize();
            // Write our CSVs
            for(int i = 0; i < _numParticles; i++) {    
                recorder.AddPayload(i);                             // particle id
                recorder.AddPayload(r.temp_particles[i].position);    // particle position
                recorder.AddPayload(r.temp_velocities[i]);            // particle velocity
                recorder.AddPayload(r.temp_densities[i]);             // particle density
                recorder.AddPayload(r.temp_pressures[i]);             // particle pressure
                recorder.WriteLine();
            }

            // Flush and Disable CSV
            recorder.Disable();

            // Pass `r`'s RecordingFrame to the new queue
            _recordFrameQueue.Enqueue(r.frameData);

            // Yield return null
            yield return null;
        }
    }

    void OnDestroy() {
        ARG_BUFFER.Release();
        SPAWN_BOUNDS_BUFFER.Release();
        BOUNDS_BUFFER.Release();
        CELL_LIMITS_BUFFER.Release();
        PARTICLE_NEIGHBORS_BUFFER.Release();
        PARTICLE_OFFSETS_BUFFER.Release();
        FORCES_BUFFER.Release();
        RENDER_LIMITS_BUFFER.Release();
        
        // Make sure all records are processed
        if (_recordCoroutine != null) {
            while (_recordQueue.Count > 0) {}
            StopCoroutine(_recordCoroutine);
        }
        // Make sure all frame data are processed
        // Initialize the recorder for frame data
        recorder.fileName = "simulation_records";
        recorder.columns = new List<string> {"frame","dt_passed_sec","rl_passed_sec","fps","n_obstacles","n_vertices","n_edges","n_triangles","filename"};
        recorder.Initialize();
        while(_recordFrameQueue.Count > 0) {
            // Get the top record in our queue
            RecordingFrame r = _recordFrameQueue.Dequeue();        
            // Write to the CSV file    
            recorder.AddPayload(r.frame);                       // current frame
            recorder.AddPayload(r.dt_passed_sec);               // Amount of delta time passed since the start
            recorder.AddPayload(r.rl_passed_sec);               // Amount of real-world time passed since the start
            recorder.AddPayload(r.fps);                         // The fps at this frame
            recorder.AddPayload(r.n_obstacles);                 // The number of obstacles
            recorder.AddPayload(r.n_vertices);                  // The number of vertices
            recorder.AddPayload(r.n_edges);                     // The number of edges
            recorder.AddPayload(r.n_triangles);                 // The number of triangles
            recorder.AddPayload(r.fileName+".csv");             // The associated filename of this particular frame
            recorder.WriteLine();
        }
        // Flush and Disable CSV
        recorder.Disable();
    }

    private void OnlineSaveCallback(WebRequests.Response response) {
        Debug.Log($"Online Save Report Success: {response.success}: {response.response}");
    }

    private int GetProjectedGridIndexFromXYZ(int x, int y, int z, int[] numCellsPerAxis) {
        return x + (numCellsPerAxis[0] * y) + (numCellsPerAxis[0] * numCellsPerAxis[1] * z);
    }

    private void DebugBufferFloat(string debugText, int debugSize, ComputeBuffer b) {
        float[] temp = new float[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferFloat3(string debugText, int debugSize, ComputeBuffer b) {
        float3[] temp = new float3[debugSize];
        b.GetData(temp);
        string top = "", bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i]}\t|";
        }
        Debug.Log($"{debugText}:\n{top}\n{bottom}");
        temp = null;
    }
    private void DebugBufferParticle(string debugText, int debugSize, ComputeBuffer b) {
        OP.Particle[] temp = new OP.Particle[debugSize];
        b.GetData(temp);
        string top = "";
        string bottom = "";
        for(int i = 0; i < debugSize; i++) {
            top += $"{i}\t|";
            bottom += $"{temp[i].position}\t|";
        }
        Debug.Log($"{debugText} positions:\n{top}\n{bottom}");
        temp = null;
    }
}
