using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;

public class Recording {
    public OP.Particle[] temp_particles;
    public float3[] temp_velocities;
    public float[] temp_velocity_magnitudes;
    public float[] temp_densities;
    public float[] temp_near_densities;
    public float[] temp_pressures;
    public float[] temp_near_pressures;
    public RecordingFrame frameData;
    public Recording(
            OP.Particle[] pa, float3[] ve, float[] vem, float[] de, float[] nde, float[] pr, float[] npr,
            int frame, float dt_passed_sec, float rl_passed_sec, float fps, 
            int n_obstacles, int n_vertices, int n_edges, int n_triangles
    ) {
        this.temp_particles = pa;
        this.temp_velocities = ve;
        this.temp_velocity_magnitudes = vem;
        this.temp_densities = de;
        this.temp_near_densities = nde;
        this.temp_pressures = pr;
        this.temp_near_pressures = npr;
        this.frameData = new RecordingFrame(
            frame, dt_passed_sec, rl_passed_sec, fps,
            n_obstacles, n_vertices, n_edges, n_triangles 
        );
    }

     public Recording(
            int frame, float dt_passed_sec, float rl_passed_sec, float fps, 
            int n_obstacles, int n_vertices, int n_edges, int n_triangles
    ) {
        this.frameData = new RecordingFrame(
            frame, dt_passed_sec, rl_passed_sec, fps,
            n_obstacles, n_vertices, n_edges, n_triangles 
        );
    }
}

public class RecordingFrame {
    public int frame;               // The simulation frame
    public float dt_passed_sec;     // The simulation time
    public float rl_passed_sec;     // The real-world time
    public float fps;               // FPS at this point in time
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

public class Simulation3D : MonoBehaviour
{
    [Header("== REFERENCES ==")]
    [SerializeField, Tooltip("Reference to the singular BufferManager component that handles all buffers in the system")]
    public BufferManager _BM;

    [SerializeField, Tooltip("Reference to the Particle Spawner we use to generate all our particles.")]
    public ParticleSpawner _SPAWNER;
    [SerializeField, Tooltip("The compute shader that performs all our SPH operations.")]
    private ComputeShader compute;

    [Header("== SPH CONFIGURATIONS ==")]
    [SerializeField, Tooltip("The delta time of the simulation. If set to any value < 0, then the simulation will default to using `Time.deltaTime`.")]
    private float _dt = 0.00825f;
    public float dt => _dt;
    [SerializeField, Tooltip("The mass of particles. In reality, pretty darn small. But we set to 1 here.")]
    private float _particleMass = 1f;
    [SerializeField, Tooltip("The external pressure - usually gravity.")]
    private Vector3 _external_force = new Vector3(0f, -9.8f, 0f);
    [SerializeField, Tooltip("The kernel radius for particle interactions. Referred to as `h` in mathematics of SPH.")]
    private float _kernel_radius = 0.25f;
    public float kernel_radius => _kernel_radius;
    public float h => _kernel_radius;
    [SerializeField, Tooltip("The `p_0` part of the SPH equations - the resting density the fluid. Also considered the 'target density' of the medium.")]
    private float _rest_density = 1000f;
    [SerializeField, Tooltip("The viscosity coefficient of the fluid. Also refered to in the SPH math as `mu`.")]
    private float _viscosity_influence = 0.5f;
    [SerializeField, Tooltip("Some refer to this as the ideal gas constant, others refer to it as the pressure influence, and some people (Bindel) refer this as the bulk modulus value. Traditionally labeled as `k` in SPH math.")]
    private float _pressure_influence = 1f;
    [SerializeField, Tooltip("`nk` - similar to pressure influence, but for really close particle-particle interactions")]
    private float _near_pressure_influence = 1f;
    [SerializeField, Tooltip("The damping effect whenever a particle hits the bounds of the simulation")]
    private float _damping_effect = -0.5f;
    [SerializeField, Tooltip("The number of particles used in the simulation. Don't touch!")]
    private int _numParticles = 0;
    public int numParticles => _numParticles;

    public enum RenderValue { Off, Velocity, Density, Pressure, Near_Pressure }
    [Header("== RENDERING & VISUALIZATION ==")]
    [SerializeField, Tooltip("How long should we delay the beginning of the simulation? This is important because we need to avoid the first frames of the unity runtime, which may be different than the rest of the simulation")]
    private float _delay = 1f;
    [SerializeField, Tooltip("How big should the particles be, visually?")]
    private float _render_size = 0.2f;
    [SerializeField, Tooltip("Should we color the particles based on the float buffer pass it, or do we just render as a blue sphere? toggle ON for coloration.")]
    private RenderValue _color_toggle = RenderValue.Off;
    [SerializeField, Tooltip("The lower limit that we normalize the values of each particle to whenm calculating what color to render.")]
    private float min_color_val = 0f;
    [SerializeField, Tooltip("The upper limit that we normalize the values of each particle to when calculating what color to render.")]
    private float max_color_val = 100f;
    [SerializeField, Tooltip("If you want to restrict which particles are rendered, where would you limit the range to be")]
    private float[] render_limits = { 0, 0, 0, 100, 100, 100 };
    [SerializeField, Tooltip("The mesh used to render each particle in the simulation. Usually just the default `Sphere` mesh from Unity.")]
    private Mesh _particle_mesh;
    [SerializeField, Tooltip("The material used to render the particles.")]
    private Material _particle_material;


    /*
    When it comes to timing, we want to keep in mind the following things:
    - The simulation may or may not be separate from the runtime of the simulation. However, it has its own time passed (either Time.fixedDeltaTime or Time.delta Time)
        and its own number of frames.
    - To calculate the timing, we also track the amount of real time that has passed. This is specifically done in the update loop using Time.deltaTime.
    - To calculate FPS, we track the number of frames passed in the update loop, and we record whenever the total amount of real time passed surpasses our intended interval
    */
    [Header("== FPS Calculation ==")]
    // The simulation time and frames
    [ReadOnly, SerializeField] private float _simulation_time_passed = 0f;
    [ReadOnly, SerializeField] private int _simulation_frames_passed = 0;
    // The real world time and frames
    [ReadOnly, SerializeField] private float _rl_time_passed = 0f;
    [ReadOnly, SerializeField] private int _rl_frames_passed = 0;
    // The variables required for FPS calculation
    [SerializeField, Tooltip("How much time we wait for each FPS update.")]  private float _fps_calculation_interval = 0.25f;
    [ReadOnly, SerializeField] private int _fps_frames_counter = 0;
    [ReadOnly, SerializeField] private float _fps_time_counter = 0f;
    [ReadOnly, SerializeField] private float _fps;


    public enum RecordingType { Off, All, FPS }
    [Header("== RECORDING ==")]
    [SerializeField, Tooltip("Reference to a csvWriter that will record our findings for us.")] private CSVWriter recorder;
    [SerializeField, Tooltip("Public toggle to determine if we should be recording.")] private RecordingType _record_statistics = RecordingType.Off;
    [SerializeField, Tooltip("The duraction (in seconds) of the recording session")]  private float _record_interval = 0.01f;
    [Tooltip("Amount of time passed since last recording update")] private float _record_interval_time = 0f;
    private IEnumerator _recordCoroutine = null;
    private Queue<Recording> _recordQueue = new Queue<Recording>();
    private Queue<RecordingFrame> _recordFrameQueue = new Queue<RecordingFrame>();
    [SerializeField] private RecordingCanvas textboxes = null;
    private bool is_recording_frame = false;


    // BUFFERS
    // We'll be using buffers primarily from our buffer manager. But some buffers we will be storing here.
    // The buffers we'll be using from the buffer manager are:
    // - particles buffer
    // - velocity buffer
    // - density buffer
    // - near_density buffer
    // - pressure buffer
    // - near_pressure buffer
    // The ones we will be maintaining here specifically are:
    // - ARG_BUFFER - for rendering particles
    // - predicted positions buffer
    // - spatial indices
    // - spatial offsets
    // The latter two are related to how we store data on the particles' positions and how we use neighbor search

    public ComputeBuffer ARG_BUFFER;
    public ComputeBuffer predictedPositionsBuffer;
    ComputeBuffer spatialIndices;
    ComputeBuffer spatialOffsets;
    ComputeBuffer renderLimitsBuffer;

    // KERNEL IDs
    int updateExternalForcesKernel;
    int updateSpatialHashKernel;
    int updateDensityKernel;
    int updatePressureKernel;
    int updateViscosityKernel;
    int updatePositionsKernel;

    // This is provided by Sebastian Lague, who actually inspired a lot of this code, admittedly.
    GPUSort gpuSort;

    // Shader Renderer kernels
    int renderer_size = Shader.PropertyToID("size");
    int renderer_particle_buffer = Shader.PropertyToID("particle_buffer");
    int renderer_float_buffer = Shader.PropertyToID("float_buffer");
    //int render_limits_property = Shader.PropertyToID("render_limits_buffer");
    int renderer_max_float = Shader.PropertyToID("max_color_val");
    int renderer_min_float = Shader.PropertyToID("min_color_val");
    int renderer_render_limits_buffer = Shader.PropertyToID("render_limits_buffer");
    int renderer_color_toggle = Shader.PropertyToID("color_toggle");

    // Other key variables
    ParticleSpawner.SpawnData spawnData;
    private int particle_threads;


    public void Initialize() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (compute == null) {
            Debug.LogError("SPH - ERROR: Cannot operate if we don't have a compute shader. Please define these references and restart the simulation.");
            return;
        }

        // Get the positions and the velocities of our particles
        // Issue: particles themselves actaully a struct, which contains not just the particle position but also the particle's acceleration (i.e. `force`) and a uint `render`.
        // Fix: let the spawner handle this step. `particles` will be the struct array we use to populate our particles buffer
        spawnData = _SPAWNER.GetSpawnData();
        _numParticles = spawnData.particles.Length;

        // Call the buffer manager to initialize particle buffers
        _BM.InitializeParticleBuffers(_numParticles);
        // Initialize Kernel IDs
        InitializeKernelIDs();
        // Now, populate buffers and connect them to the compute shader
        InitializeBuffers(spawnData);
        // Pass key variables to the compute shader
        InitializeShaderVariables(_dt);

        // If we are recording, we initialize a recording session
        if (_record_statistics != RecordingType.Off) {
            _recordCoroutine = RecordSim();
            StartCoroutine(_recordCoroutine);
        }

        if (textboxes != null) {
            textboxes.nParticles.gameObject.SetActive(true);
            textboxes.nObstacles.gameObject.SetActive(true);
            textboxes.obstacleDetails.gameObject.SetActive(true);
            textboxes.fps.gameObject.SetActive(true);
            textboxes.deltaTime.gameObject.SetActive(true);
            textboxes.realTime.gameObject.SetActive(true);
        }
    }

    private void InitializeKernelIDs() {
        updateExternalForcesKernel = compute.FindKernel("UpdateExternalForces");
        updateSpatialHashKernel = compute.FindKernel("UpdateSpatialHash");
        updateDensityKernel = compute.FindKernel("UpdateDensity");
        updatePressureKernel = compute.FindKernel("UpdatePressure");
        updateViscosityKernel = compute.FindKernel("UpdateViscosity");
        updatePositionsKernel = compute.FindKernel("UpdatePositions");
    }

    private void InitializeBuffers(ParticleSpawner.SpawnData spawnData) {

        // Create buffers
        // Remember: we'll be using buffers from our buffer manager. However, we need to define these buffers because these are unique to the setup here.
        predictedPositionsBuffer = ComputeHelper.CreateStructuredBuffer<float3>(_numParticles);
        spatialIndices = ComputeHelper.CreateStructuredBuffer<uint3>(_numParticles);
        spatialOffsets = ComputeHelper.CreateStructuredBuffer<uint>(_numParticles);
        renderLimitsBuffer = ComputeHelper.CreateStructuredBuffer<float>(6);

        // Set buffer data from the spawner data to populate our particle positions and velocities.
        // We'll also update our predicted positiosn too
        _BM.PARTICLES_BUFFER.SetData(spawnData.particles);
        _BM.PARTICLES_VELOCITIES_BUFFER.SetData(spawnData.velocities);
        predictedPositionsBuffer.SetData(spawnData.positions);
        renderLimitsBuffer.SetData(render_limits);

        // Now, we have to start attaching buffers to kernels.
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_BUFFER, "particles", updateExternalForcesKernel, updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, predictedPositionsBuffer, "predictedPositions", updateExternalForcesKernel, updateSpatialHashKernel, updateDensityKernel, 
                                                                                            updatePressureKernel, updateViscosityKernel, updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, spatialIndices, "spatialIndices", updateSpatialHashKernel, updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, spatialOffsets, "spatialOffsets", updateSpatialHashKernel, updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_DENSITIES_BUFFER, "densities", updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_NEAR_DENSITIES_BUFFER, "near_densities", updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_PRESSURE_BUFFER, "pressures", updatePressureKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_NEAR_PRESSURE_BUFFER, "near_pressures", updatePressureKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_VELOCITIES_BUFFER, "velocities", updateExternalForcesKernel, updatePressureKernel, updateViscosityKernel, updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_VELOCITY_MAGNITUDES_BUFFER, "velocity_magnitudes", updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_EXTERNAL_FORCES_BUFFER, "obstacleForces", updatePositionsKernel);

        // Now, let gpuSort look at spatialindicse and spatialoffsets.
        gpuSort = new();
        gpuSort.SetBuffers(spatialIndices, spatialOffsets);

        // Let's take a quickie to calculate the number of threads needed to run each dispatch
        particle_threads = Mathf.CeilToInt((float)_numParticles/256f);

        // Lastly, initialize variables for our renderer
        uint[] arg = {_particle_mesh.GetIndexCount(0), (uint)(_numParticles), _particle_mesh.GetIndexStart(0), _particle_mesh.GetBaseVertex(0), 0};
        ARG_BUFFER = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        ARG_BUFFER.SetData(arg);
    }

    private void InitializeShaderVariables(float timestep) {
        // Note: this is NOT updated every frame
        compute.SetInt("numParticles", _numParticles);
        compute.SetFloat("epsilon", Mathf.Epsilon);
        compute.SetFloat("pi", Mathf.PI);
        // We call `UpdateShaderVariables()` because that's kind of important too.
        UpdateShaderVariables(timestep);
    }

    private void UpdateShaderVariables(float timestep) {
        // Note: this is called every frame, so you can update the variables in the inspector and they'll also be updated in the shader.
        compute.SetFloat("particleMass", _particleMass);
        compute.SetFloat("deltaTime", timestep);
        compute.SetVector("external_force", _external_force);   // gravity
        compute.SetFloat("smoothingRadius", _kernel_radius);      // h
        compute.SetFloat("rest_density", _rest_density);        // p_0
        compute.SetFloat("viscosity_influence", _viscosity_influence);  // mu
        compute.SetFloat("pressure_influence", _pressure_influence);    // bulk modulus
        compute.SetFloat("near_pressure_influence", _near_pressure_influence);
        compute.SetFloat("damping_effect", _damping_effect);
        renderLimitsBuffer.SetData(render_limits);
    }

    private float delay_time_passed = 0f;


    private void Update() {
        // // Run the simulation step if we are not using fixed timestamp
        SimulationStep(_dt);

        // Don't continue anymore if we haven't passed the amount of time delay
        if (delay_time_passed < _delay) return;

        // Update the real amount of time that passed
        _rl_time_passed += Time.deltaTime;
        _rl_frames_passed += 1;

        // Update FPS
        UpdateFPS();

        // We have to render the particles now, if we are sufficiently past the delay time
        // In our materials buffer, let's set the necessary buffers for rendering particles and the like
        _particle_material.SetBuffer(renderer_particle_buffer, _BM.PARTICLES_BUFFER);
        _particle_material.SetBuffer(renderer_render_limits_buffer, renderLimitsBuffer);
        _particle_material.SetFloat(renderer_size, _render_size);
        switch(_color_toggle) {
            case RenderValue.Density:
                _particle_material.SetBuffer(renderer_float_buffer, _BM.PARTICLES_DENSITIES_BUFFER);
                _particle_material.SetInt(renderer_color_toggle, 1);
                break;
            case RenderValue.Velocity:
                _particle_material.SetBuffer(renderer_float_buffer, _BM.PARTICLES_VELOCITY_MAGNITUDES_BUFFER);
                _particle_material.SetInt(renderer_color_toggle, 1);
                break;
            case RenderValue.Pressure:
                _particle_material.SetBuffer(renderer_float_buffer, _BM.PARTICLES_PRESSURE_BUFFER);
                _particle_material.SetInt(renderer_color_toggle, 1);
                break;
            case RenderValue.Near_Pressure:
                _particle_material.SetBuffer(renderer_float_buffer, _BM.PARTICLES_NEAR_PRESSURE_BUFFER);
                _particle_material.SetInt(renderer_color_toggle, 1);
                break;
            default:
                _particle_material.SetInt(renderer_color_toggle, 0);
                break;
        }
        _particle_material.SetFloat(renderer_min_float, min_color_val);
        _particle_material.SetFloat(renderer_max_float, max_color_val);

        Graphics.DrawMeshInstancedIndirect(
                    _particle_mesh, 
                    0, 
                    _particle_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
        
        // Update visualization textboxes if they exist
        if (textboxes != null) {
            textboxes.nParticles.text = $"{_numParticles} Particles";
            textboxes.nObstacles.text = $"{_BM.MESHOBS_OBSTACLES_STATIC_BUFFER.count} Obstacles";
            textboxes.obstacleDetails.text = $"(V:{_BM.MESHOBS_VERTICES_STATIC_BUFFER.count} | E:{_BM.MESHOBS_EDGES_STATIC_BUFFER.count} | T:{_BM.MESHOBS_TRIANGLES_STATIC_BUFFER.count})";
            textboxes.fps.text = $"{_fps} FPS";
            textboxes.deltaTime.text = $"Simulation time: { _simulation_time_passed } sec.";
            textboxes.realTime.text = $"Real time: { _rl_time_passed } sec.";
        }

        if (is_recording_frame) return;
        // If we're recording, record our session
        // We based recording time based on the provided timestep. So if we use a fixed timestep, our recordings will be based on the fixed timestep
        // Alternatively, if we wer called via Update and not FixedUpdate(), then we use the variable timestep;
        _record_interval_time += Time.deltaTime;
        if (_record_interval_time >= _record_interval) {
            // Here, we want to capture some important details
            // - We will instantiate a new file for every recording we capture. 
            // - We will time each capture on a time interval, meaning we will at least ensure that we won't lag the system.
            // - Each file will contain the following:
            //      1. the filename should contain the recorded timestamp (total elapsed time), the frame number itself, and the fps
            //      2. all subsequent lines should represent a particle each
            //      3. For each line/particle, we'll record their ID, position, velocity, density, and pressure value
            
            print("Making Record");
            MakeRecord();
            _record_interval_time = 0;
        }

    }

    private void SimulationStep(float timestep) {

        if (delay_time_passed < _delay) {
            delay_time_passed += timestep;
            return;
        }

        if (is_recording_frame) return;

        // Update parameters, including the timestep
        UpdateShaderVariables(timestep);

        // Apply external forces (aka GRAVITY) + make predictions of particle positions
        compute.Dispatch(updateExternalForcesKernel, particle_threads, 1, 1);
        // Reset spatial offsets, recalibrate spatial hashes
        compute.Dispatch(updateSpatialHashKernel, particle_threads, 1, 1);
        // Let the CPU and GPU team up to sort indices of particles in cells and optimize search later on
        gpuSort.SortAndCalculateOffsets();
        // Upate density of particles
        compute.Dispatch(updateDensityKernel, particle_threads, 1, 1);
        // Update pressure of particles
        compute.Dispatch(updatePressureKernel, particle_threads, 1, 1);
        // Update viscosity of particles
        compute.Dispatch(updateViscosityKernel, particle_threads, 1, 1);
        // Update the true positions of the particles
        compute.Dispatch(updatePositionsKernel, particle_threads, 1, 1);

        // Update the time elapsed and frames elapsed specific to the simulation
        _simulation_frames_passed += 1;
        _simulation_time_passed += timestep;
    }

    private void UpdateFPS() {
        _fps_frames_counter += 1;
        _fps_time_counter += Time.deltaTime;
        if (_fps_time_counter > _fps_calculation_interval) {
            _fps = _fps_frames_counter / _fps_time_counter;
            _fps_frames_counter = 0;
            _fps_time_counter = 0f;
        }
    }

    private IEnumerator RecordSim() {
        // Make our initial recording
        MakeRecord();

        // Perform the while loop to continuously record until the end of the simulation
        while(true) {
            // Don't do anything if our queue is empty
            if ( _recordQueue.Count == 0) {
                yield return null;
                continue;
            }

            // Get the top record in our queue
            Recording r = _recordQueue.Dequeue();

            // If we want to write JUST the FPS, we don't record the particle positions.
            if (_record_statistics == RecordingType.All) {
                // Initialize the recorder
                recorder.fileName = r.frameData.fileName;
                recorder.columns = new List<string> {
                    "particle_id",
                    "position_x", "position_y", "position_z",
                    "velocity_x", "velocity_y", "velocity_z", "velocity_mag",
                    "density", "near_density", "pressure", "near_pressure"
                };
                recorder.Initialize();
                // Write our CSVs
                for(int i = 0; i < _numParticles; i++) {    
                    recorder.AddPayload(i);                               // particle id
                    recorder.AddPayload(r.temp_particles[i].position);    // particle position
                    recorder.AddPayload(r.temp_velocities[i]);            // particle velocity
                    recorder.AddPayload(r.temp_velocity_magnitudes[i]);   // particle velocity magnitudes
                    recorder.AddPayload(r.temp_densities[i]);             // particle density
                    recorder.AddPayload(r.temp_near_densities[i]);        // particle near density
                    recorder.AddPayload(r.temp_pressures[i]);             // particle pressure
                    recorder.AddPayload(r.temp_near_pressures[i]);        // particle near pressures
                    recorder.WriteLine();
                }

                // Flush and Disable CSV
                recorder.Disable();
            }

            // Pass `r`'s RecordingFrame to the new queue
            _recordFrameQueue.Enqueue(r.frameData);

            // Yield return null
            yield return null;
        }
    }

    private void MakeRecord() {
        if (_record_statistics == RecordingType.Off) return;

        is_recording_frame = true;

         // Remember that each of these dynamic and static variations of the vertices, triangles, and edges all are the same size. So it doesn't matter if we get the number of them from either the static or dyanmic buffers
        int n_obstacles = (_BM.MESHOBS_OBSTACLES_STATIC_BUFFER != null) ? _BM.MESHOBS_OBSTACLES_STATIC_BUFFER.count : 0;
        int n_vertices = (_BM.MESHOBS_VERTICES_STATIC_BUFFER != null) ? _BM.MESHOBS_VERTICES_STATIC_BUFFER.count : 0;
        int n_edges = (_BM.MESHOBS_EDGES_STATIC_BUFFER != null) ? _BM.MESHOBS_EDGES_STATIC_BUFFER.count : 0;
        int n_triangles = (_BM.MESHOBS_TRIANGLES_STATIC_BUFFER != null) ? _BM.MESHOBS_TRIANGLES_STATIC_BUFFER.count : 0;

        // GET DATA, but only if required
        if (_record_statistics == RecordingType.All) {
            OP.Particle[] temp_particles = new OP.Particle[_numParticles];
                _BM.PARTICLES_BUFFER.GetData(temp_particles);
            float3[] temp_velocities = new float3[_numParticles];
                _BM.PARTICLES_VELOCITIES_BUFFER.GetData(temp_velocities);
            float[] temp_velocity_magnitudes = new float[_numParticles];
                _BM.PARTICLES_VELOCITY_MAGNITUDES_BUFFER.GetData(temp_velocity_magnitudes);
            float[] temp_densities = new float[_numParticles];
                _BM.PARTICLES_DENSITIES_BUFFER.GetData(temp_densities);
            float[] temp_near_densities = new float[_numParticles];
                _BM.PARTICLES_NEAR_DENSITIES_BUFFER.GetData(temp_near_densities);
            float[] temp_pressures = new float[_numParticles];
                _BM.PARTICLES_PRESSURE_BUFFER.GetData(temp_pressures);
            float[] temp_near_pressures = new float[_numParticles];
                _BM.PARTICLES_NEAR_PRESSURE_BUFFER.GetData(temp_near_pressures);

            // Create a new entry in a queue that holds all the recordings    
            _recordQueue.Enqueue(new Recording(
                temp_particles, 
                temp_velocities, temp_velocity_magnitudes, 
                temp_densities, temp_near_densities, 
                temp_pressures, temp_near_pressures, 
                _simulation_frames_passed, _simulation_time_passed, _rl_time_passed, _fps,
                n_obstacles, n_vertices, n_edges, n_triangles
            ));
        }
        // However, if we're just recording the fps, no need to do all the above
        else {
            _recordQueue.Enqueue(new Recording(
                _simulation_frames_passed, _simulation_time_passed, _rl_time_passed, _fps,
                n_obstacles, n_vertices, n_edges, n_triangles
            ));
        }
        

        is_recording_frame = false;
    }

    void OnDestroy() {
        // While buffer manager handles releasing buffers, we also have to release our own.
        ARG_BUFFER.Release();
        if (predictedPositionsBuffer != null) predictedPositionsBuffer.Release();
        if (spatialIndices != null) spatialIndices.Release();
        if (spatialOffsets != null) spatialOffsets.Release();
        renderLimitsBuffer.Release();
        
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

}
