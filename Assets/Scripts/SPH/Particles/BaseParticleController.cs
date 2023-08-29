using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class BaseParticleController : MonoBehaviour
{
    [SerializeField] private Grid _GRID;
    [SerializeField] private int _section_index = -1;
    [SerializeField] private ComputeShader _SHADER;
    private int _num_blocks_grid, _num_blocks_particles;

    [Header("== Particle Generation ==")]
    [SerializeField] private float _particle_render_radius = 1f;
    public float particleRenderRadius => _particle_render_radius;
    public float particle_render_radius => _particle_render_radius;
    [SerializeField, Tooltip("How much space is needed between newly-generated particles?")] private float _spawn_distance_btw_particles = 1f;
    [SerializeField, ReadOnly, Tooltip("A READ-ONLY array that identifies how many particles can be spawned per axis")] private int[] _num_particles_per_axis;
    [SerializeField, ReadOnly] private int _num_particles_per_grid_cell;
    [SerializeField, ReadOnly] private float[] _spawn_bounds;
    [SerializeField, ReadOnly] private int _max_num_particles;
    [SerializeField] private int _num_particles;

    [Header("== SPH Variables ==")]
    [SerializeField] private float _particle_mass = 0.1f;
    [SerializeField] private float _smoothingKernel = 1f;
    public float smoothingKernel => _smoothingKernel;
    public float h => _smoothingKernel;
    [SerializeField] private float _rest_density = 1000f;
    [SerializeField] private float _viscosity_coefficient = 0.5f;
    public float viscosity_coefficient => _viscosity_coefficient;
    public float mu => _viscosity_coefficient;
    [SerializeField] private float _bulk_modulus = 1.25f;
    public float bulk_modulus => _bulk_modulus;
    public float k => _bulk_modulus;
    [SerializeField] private float[] _g = new float[3] {0f, -9.81f, 0f};
    public float[] g => _g;
    [SerializeField] private float _dt = 0.0165f;
    public float dt => _dt;

    [Header("== Shader kernels ==")]
    private int _CLEAR_GRID;
    private int _CLEAR_FORCES;
    private int _GENERATE_PARTICLES;
    private int _UPDATE_GRID;
    private int _COMPUTE_DENSITY;
    private int _COMPUTE_INTERNAL_FORCES;
    private int _INTEGRATE, _INTEGRATE_DEBUG;
    private int _DAMPEN_BY_BOUNDS;

    [Header("== Compute Buffers ==")]
    public ComputeBuffer ARG_BUFFER;
    public ComputeBuffer CELL_LIMITS_BUFFER;
    public ComputeBuffer BOUNDS_BUFFER, SPAWN_BOUNDS_BUFFER;
    public ComputeBuffer PARTICLE_NEIGHBORS_BUFFER, PARTICLE_OFFSETS_BUFFER;
    public ComputeBuffer FORCES_BUFFER;
    private ComputeBuffer GRID_CELL_RENDER_LIMITS_BUFFER;

    /// <summary>
    /// DESCRIPTION: This is called if the _GRID component is updated at any point. This specific function updates key variables that will be used later.
    /// INPUT: (none)
    /// OUTPUT: (none)
    /// STEPS:
    ///     1. Check to make sure that the _GRID reference is established
    ///     2. Make sure that the `_spawn_distance_btw_particles` variable is capped so that it's always bigger than 0
    ///     3. Determine the number of particles that can be spawned per axis
    ///     4. Determine what boundaries we are spawning inside, based on whether `_section_index` is -1 or a reference ID
    ///     5. Store number of particles to be spawned per axis
    ///     6. Update the maximum number of particles that can fit in the simulation space, as well as clamp the user-defined number of particles
    ///     7. Calculate how many particles can actually fit within a single grid cell
    ///     8. Identify the spawning boundaries for particle generation during startup
    /// </summary>
    public void GridUpdated() {
        ///     1. Check to make sure that the _GRID reference is established
        if (_GRID == null) {
            Debug.LogError("ERROR: GridUpdated() callback missing a `_GRID` reference. Please set this and update the grid once again.");
            return;
        }

        ///     2. Make sure that the `_spawn_distance_btw_particles` variable is capped so that it's always bigger than 0
        _spawn_distance_btw_particles = Mathf.Clamp(_spawn_distance_btw_particles,_particle_render_radius,_smoothingKernel);
        
        ///     3. Determine the number of particles that can be spawned per axis
        int numParticlesPerCellAxis = Mathf.CeilToInt(_GRID.gridCellSize / _spawn_distance_btw_particles);

        ///     4. Determine what boundaries we are spawning inside, based on whether `_section_index` is -1 or a reference ID
        Vector3 bs = (_section_index >= 0) 
            ? _GRID.sections[_section_index].dimensionsV3 
            : _GRID.innerBoundsV3;

        ///     5. Store number of particles to be spawned per axis
        _num_particles_per_axis = new int[3] {
            Mathf.FloorToInt(bs[0] / _spawn_distance_btw_particles),
            Mathf.FloorToInt(bs[1] / _spawn_distance_btw_particles),
            Mathf.FloorToInt(bs[2] / _spawn_distance_btw_particles)
        };

        ///     6. Update the maximum number of particles that can fit in the simulation space, as well as clamp the user-defined number of particles
        _max_num_particles = _num_particles_per_axis[0] * _num_particles_per_axis[1] * _num_particles_per_axis[2];
        _num_particles = Mathf.Clamp(_num_particles,1,_max_num_particles);

        ///     7. Calculate how many particles can actually fit within a single grid cell
        _num_particles_per_grid_cell = numParticlesPerCellAxis * _GRID.numCellsPerAxis[0] 
            + numParticlesPerCellAxis * _GRID.numCellsPerAxis[1]
            + numParticlesPerCellAxis * _GRID.numCellsPerAxis[2];

        ///     8. Identify the spawning boundaries for particle generation during startup
        _spawn_bounds = (_section_index >= 0) 
            ? _GRID.sections[_section_index].bounds 
            : _GRID.innerBounds;
    }
    
    /// <summary>
    /// This function is called to initialize the SPH particle system. This MUST be called prior to anything else
    /// INPUT: (none)
    /// OUTPUT: boolean - tells whatever called this function if the particle controller was successful in running
    /// STEPS:
    ///     1. End early if either the _GRID or _SHADER references are not set
    ///     2. Initialize key variables used only within this class component
    ///     3. Initialize key variables used in the _SHADER compute shader
    ///     4. Initialize the kernels that reference functions in _SHADER
    /// </summary>
    public bool Initialize() {
        ///     1. End early if either the _GRID or _SHADER references are not set
        if (_GRID == null || _SHADER == null) {
            Debug.LogError("ERROR: `GRID` or `SHADER` not set for the Particle Controller. Please add these references and run again.");
            return false;
        }

        ///     2. Initialize key variables used only within this class component
        InitializeCPUVariables();

        ///     3. Initialize key variables used in the _SHADER compute shader. This includes static and dynamic variables
        InitializeStaticGPUVariables();

        ///     4. Initialize the kernels that reference functions in _SHADER
        //InitializeKernels();

        return true;
    }

    /// <summary>
    /// DESCRIPTION: Initializes CPU variables that are required explicitly only within the BaseParticleController class and nowhere else.
    /// INPUT: (none)
    /// OUTPUT: (none)
    /// </summary>
    private void InitializeCPUVariables() {
        // Determine block number for GPU threading
        _num_blocks_grid = Mathf.CeilToInt((float)_GRID.numGridCells / 512f);
        _num_blocks_particles = Mathf.CeilToInt((float)_num_particles / 512f);
    }

    /// <summary>
    /// DESCRIPTION: Initializes static GPU variables that are required explicitly within the _SHADER compute shader.
    /// INPUTS: (none)
    /// OUTPUTS: (none)
    /// </summary>
    private void InitializeStaticGPUVariables() {
        // == WORLD CONFIGURATIONS ==
        _SHADER.SetFloat("gridCellSize", _GRID.gridCellSize);
        _SHADER.SetInts("numCellsPerAxis", _GRID.numCellsPerAxis);
        _SHADER.SetInt("total_number_of_cells", _GRID.numGridCells);
        _SHADER.SetFloat("epsilon", Mathf.Epsilon);
        _SHADER.SetFloat("pi", Mathf.PI);

        _SHADER.SetFloat("spawnDistanceBetweenParticles",_spawn_distance_btw_particles);
        _SHADER.SetFloat("gridScalingX", _GRID.gridScaling[0]);
        _SHADER.SetFloat("gridScalingY", _GRID.gridScaling[1]);
        _SHADER.SetFloat("gridScalingZ", _GRID.gridScaling[2]);

        // == PARTICLE CONFIGURATIONS ==
        _SHADER.SetInt("numParticles", _num_particles);
        _SHADER.SetInts("numParticlesPerAxis", _num_particles_per_axis);
        _SHADER.SetInt("numParticlesPerGridCell", _num_particles_per_grid_cell);

        // == GPU SETTINGS
        _SHADER.SetInt("numBlocks", _num_blocks_grid);

        // Update variables that may change over time due to modifying inspector values
        UpdateDynamicGPUVariables();
    }

    /// <summary>
    /// DESCRIPTION: Both initializes and updates dynamic GPU variables that are required explicitly within the _SHADER compute shader.
    /// INPUTS: bool = should the delta time be updated too? Default = TRUE
    /// OUTPUTS: (none)
    /// </summary>
    private void UpdateDynamicGPUVariables(bool updateDT = true) {
        // == PARTICLE CONFIGURATIONS ==
        _SHADER.SetFloat("particleRenderRadius", _particle_render_radius);

        // == FLUID MECHANICS ==
        _SHADER.SetFloat("smoothingRadius", _smoothingKernel);
        _SHADER.SetFloat("radius2", _smoothingKernel * _smoothingKernel);
        _SHADER.SetFloat("radius3", Mathf.Pow(_smoothingKernel,3));
        _SHADER.SetFloat("radius5", Mathf.Pow(_smoothingKernel,5));
        _SHADER.SetFloat("radius6", Mathf.Pow(_smoothingKernel,6));
        _SHADER.SetFloat("radius8", Mathf.Pow(_smoothingKernel,8));
        _SHADER.SetFloat("radius9", Mathf.Pow(_smoothingKernel,9));

        _SHADER.SetFloat("particleMass", _particle_mass);
        _SHADER.SetFloat("rest_density", _rest_density);
        _SHADER.SetFloat("viscosity_coefficient", _viscosity_coefficient);
        _SHADER.SetFloat("bulkModulus", _bulk_modulus);
        _SHADER.SetFloats("g", _g);

        if (updateDT) {
            float t = (_dt >= 0) ? _dt : Time.deltaTime;
            _SHADER.SetFloat("dt", t);
        }
    }

    /// <summary>
    /// DESCRIPTION: Initializes the pragma kernels that reference functions in the compute shader _SHADER
    /// INPUTS: (none)
    /// OUTPUTS: (none)
    /// </summary>
    private void InitializeKernels() {
        // Used on initialization. _CLEAR_GRID also is performed at the beginning of each update loop
        _CLEAR_GRID = _SHADER.FindKernel("ClearGrid");
        _CLEAR_FORCES = _SHADER.FindKernel("ClearForces");
        _GENERATE_PARTICLES = _SHADER.FindKernel("GenerateParticles");
        // These are run during each update loop
        _UPDATE_GRID = _SHADER.FindKernel("UpdateGridCellCounts");
        _COMPUTE_DENSITY = _SHADER.FindKernel("CV_ComputeDensity");
        _COMPUTE_INTERNAL_FORCES = _SHADER.FindKernel("ComputeInternalForces");
        _INTEGRATE = _SHADER.FindKernel("Integrate");
        _INTEGRATE_DEBUG = _SHADER.FindKernel("Integrate_Debug");
        _DAMPEN_BY_BOUNDS = _SHADER.FindKernel("DampenByBounds");
    }

    /*
    private void InitializeBuffers() {
        uint[] arg = {_GRID.particle_mesh.GetIndexCount(0), (uint)(_numParticles), _GRID.particle_mesh.GetIndexStart(0), _GRID.particle_mesh.GetBaseVertex(0), 0};
        ARG_BUFFER = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        ARG_BUFFER.SetData(arg);
        
        _BM.PARTICLES_GRID_BUFFER = new ComputeBuffer(_GRID.numGridCells, sizeof(int)*2 + sizeof(float)*3);
        BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
        BOUNDS_BUFFER.SetData(_GRID.outerBounds);
        GRID_CELL_RENDER_LIMITS_BUFFER = new ComputeBuffer(6, sizeof(float));
        _gridCellRenderLimits[0] = Mathf.Max(_gridCellRenderLimits[0],_GRID.outerBounds[0]);
        _gridCellRenderLimits[1] = Mathf.Max(_gridCellRenderLimits[1],_GRID.outerBounds[1]);
        _gridCellRenderLimits[2] = Mathf.Max(_gridCellRenderLimits[2],_GRID.outerBounds[2]);
        _gridCellRenderLimits[3] = Mathf.Min(_gridCellRenderLimits[3],_GRID.outerBounds[3]);
        _gridCellRenderLimits[4] = Mathf.Min(_gridCellRenderLimits[4],_GRID.outerBounds[4]);
        _gridCellRenderLimits[5] = Mathf.Min(_gridCellRenderLimits[5],_GRID.outerBounds[5]);
        GRID_CELL_RENDER_LIMITS_BUFFER.SetData(_gridCellRenderLimits);

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

        _BM.PARTICLES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float)*6+sizeof(int));
        debug_particles_array = new OP.Particle[_numParticles];
        if (_DEBUG_MODE) {
            // If we are in debug mode, we pre-emptively fill the particles buffer with the positions of our debug particles
            UpdateDebugParticles();
        }
        if (_POINT_CLOUD_OBSTACLE_MANAGER != null && _NUM_BOUNDARY_PARTICLES > 0) {
            Debug.Log($"{_POINT_CLOUD_OBSTACLE_MANAGER.boundaryParticles.Count}");
            _BM.PARTICLES_BUFFER.SetData(_POINT_CLOUD_OBSTACLE_MANAGER.boundaryParticles.ToArray(), 0, _numParticles, _NUM_BOUNDARY_PARTICLES);
        }
        Debug.Log($"Number of Particles: {_numParticles + _NUM_BOUNDARY_PARTICLES}");
        PARTICLE_OFFSETS_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(int));
        PARTICLE_NEIGHBORS_BUFFER = new ComputeBuffer(_GRID.numGridCells * _numParticlesPerGridCell, sizeof(int));

        _BM.PARTICLES_DENSITIES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float));
        _BM.PARTICLES_VELOCITIES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, 3 * sizeof(float));
        // We attempt to... "adjust" the density of the boundary particle ssuch tath they are able to really have an effect on fluid particles that hit them

        // We don't need to update the boundary aprticles for these
        _BM.PARTICLES_PRESSURE_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
        _BM.PARTICLES_VISCOSITY_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
        FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);

        _BM.PARTICLES_EXTERNAL_FORCES_BUFFER = new ComputeBuffer(_numParticles, sizeof(uint) + sizeof(int)*8 + sizeof(float)*27);
        OP.Projection[] projections = new OP.Projection[_numParticles];
        for(int i = 0; i < _numParticles; i++) {
            projections[i] = new OP.Projection();
            projections[i].projection = new(0f,0f,0f);
            projections[i].position = new(0f,0f,0f);
            projections[i].normal = new(0,0,0);
            projections[i].particle_force = new(0f,0f,0f);
            projections[i].external_force = new(0,0,0);
            projections[i].counter = 0;
            projections[i].frictionCoefficient = 0;
        }
        _BM.PARTICLES_EXTERNAL_FORCES_BUFFER.SetData(projections);

        // Setting the buffers

        _SHADER.SetBuffer(_CLEAR_GRID, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SHADER.SetBuffer(_CLEAR_GRID, "bounds", BOUNDS_BUFFER);

        _SHADER.SetBuffer(_CLEAR_FORCES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SHADER.SetBuffer(_CLEAR_FORCES, "force", FORCES_BUFFER);
        _SHADER.SetBuffer(_CLEAR_FORCES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);

        _SHADER.SetBuffer(_GENERATE_PARTICLES, "bounds", SPAWN_BOUNDS_BUFFER);
        _SHADER.SetBuffer(_GENERATE_PARTICLES, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_GENERATE_PARTICLES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SHADER.SetBuffer(_GENERATE_PARTICLES, "force", FORCES_BUFFER);
        _SHADER.SetBuffer(_GENERATE_PARTICLES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);

        _SHADER.SetBuffer(_UPDATE_GRID, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_UPDATE_GRID, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SHADER.SetBuffer(_UPDATE_GRID, "particleOffsets", PARTICLE_OFFSETS_BUFFER);
        _SHADER.SetBuffer(_UPDATE_GRID, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);

        _SHADER.SetBuffer(_COMPUTE_DENSITY, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_DENSITY, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_DENSITY, "cellLimits", CELL_LIMITS_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_DENSITY, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_DENSITY, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_DENSITY, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);

        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "cellLimits", CELL_LIMITS_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SHADER.SetBuffer(_COMPUTE_INTERNAL_FORCES, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);

        _SHADER.SetBuffer(_INTEGRATE, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "force", FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        //_SHADER.SetBuffer(_INTEGRATE, "bounds", BOUNDS_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE, "projections", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);

        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "force", FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SHADER.SetBuffer(_INTEGRATE_DEBUG, "density", _BM.PARTICLES_DENSITIES_BUFFER);

        _SHADER.SetBuffer(_DAMPEN_BY_BOUNDS, "bounds", BOUNDS_BUFFER);
        _SHADER.SetBuffer(_DAMPEN_BY_BOUNDS, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_DAMPEN_BY_BOUNDS, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
    }
    */

    /// <summary>
    /// This function generates particles within a given float[4] (for 2D) or float[6] (for 3D) space.
    /// To do this, it requires access to 
    /// </summary>
    public void GenerateParticles() {

    }
}
