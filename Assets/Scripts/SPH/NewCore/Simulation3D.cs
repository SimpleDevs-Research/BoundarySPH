using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

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

    [Header("== RENDERING & VISUALIZATION ==")]
    [SerializeField, Tooltip("How long should we delay the beginning of the simulation? This is important because we need to avoid the first frames of the unity runtime, which may be different than the rest of the simulation")]
    private float _delay = 1f;
    [SerializeField, Tooltip("How big should the particles be, visually?")]
    private float _render_size = 0.2f;
    [SerializeField, Tooltip("Should we color the particles based on the float buffer pass it, or do we just render as a blue sphere? toggle ON for coloration.")]
    private bool _color_toggle = true;
    [SerializeField, Tooltip("The upper limit that we normalize the values of each particle to when calculating what color to render.")]
    private float max_color_val = 100f;
    [SerializeField, Tooltip("The mesh used to render each particle in the simulation. Usually just the default `Sphere` mesh from Unity.")]
    private Mesh _particle_mesh;
    [SerializeField, Tooltip("The material used to render the particles.")]
    private Material _particle_material;

    [Header("== RECORDING ==")]
    [Tooltip("The amount of simulation time that has passed since the last recording.")]
    private float _dt_passed = 0f;
    [Tooltip("The amount of world time that has passed since the last recording.")]
    private float _real_time_elapsed = 0f;
    [Tooltip("The number of Unity frames that passed since the last recording.")]
    private int _frames_elapsed = 0;
    [Tooltip("The total amount of simulation time that has passed since the start of the simulation.")]
    private float _total_dt_elapsed = 0f;
    [Tooltip("The total amount of world time that has passed since the start of the simulation")]
    private float _total_real_time_elapsed = 0f;
    [Tooltip("The total number of Unity frames that passed since the start of the simulation")]
    private int _total_frames_elapsed = 0;
    [Tooltip("The measured FPS of the simulation.")]
    private float _fps = 0f;

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
        InitializeShaderVariables();

        /*
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
        */
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

        // Set buffer data from the spawner data to populate our particle positions and velocities.
        // We'll also update our predicted positiosn too
        _BM.PARTICLES_BUFFER.SetData(spawnData.particles);
        _BM.PARTICLES_VELOCITIES_BUFFER.SetData(spawnData.velocities);
        predictedPositionsBuffer.SetData(spawnData.positions);

        // Now, we have to start attaching buffers to kernels.
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_BUFFER, "particles", updateExternalForcesKernel, updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, predictedPositionsBuffer, "predictedPositions", updateExternalForcesKernel, updateSpatialHashKernel, updateDensityKernel, 
                                                                                            updatePressureKernel, updateViscosityKernel, updatePositionsKernel);
        ComputeHelper.SetBuffer(compute, spatialIndices, "spatialIndices", updateSpatialHashKernel, updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, spatialOffsets, "spatialOffsets", updateSpatialHashKernel, updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_DENSITIES_BUFFER, "densities", updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_NEAR_DENSITIES_BUFFER, "near_densities", updateDensityKernel, updatePressureKernel, updateViscosityKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_PRESSURE_BUFFER, "pressures", updatePressureKernel);
        ComputeHelper.SetBuffer(compute, _BM.PARTICLES_VELOCITIES_BUFFER, "velocities", updateExternalForcesKernel, updatePressureKernel, updateViscosityKernel, updatePositionsKernel);
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

    private void InitializeShaderVariables() {
        // Note: this is NOT updated every frame
        compute.SetInt("numParticles", _numParticles);
        compute.SetFloat("epsilon", Mathf.Epsilon);
        compute.SetFloat("pi", Mathf.PI);
        // We call `UpdateShaderVariables()` because that's kind of important too.
        UpdateShaderVariables();
    }

    private void UpdateShaderVariables() {
        // Note: this is called every frame, so you can update the variables in the inspector and they'll also be updated in the shader.
        compute.SetFloat("particleMass", _particleMass);
        compute.SetFloat("deltaTime", _dt);
        compute.SetVector("external_force", _external_force);   // gravity
        compute.SetFloat("smoothingRadius", _kernel_radius);      // h
        compute.SetFloat("rest_density", _rest_density);        // p_0
        compute.SetFloat("viscosity_influence", _viscosity_influence);  // mu
        compute.SetFloat("pressure_influence", _pressure_influence);    // bulk modulus
        compute.SetFloat("near_pressure_influence", _near_pressure_influence);
        compute.SetFloat("damping_effect", _damping_effect);
    }

    private float delay_time_passed = 0f;
    private void Update() {
        // Check if we're still delayed. Discontinue if we aren't passed the delay point yet.
        if (delay_time_passed < _delay) {
            delay_time_passed += Time.deltaTime;
            return;
        }

        // Update parameters
        UpdateShaderVariables();

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

        // We have to render the particles now
        // In our materials buffer, let's set the necessary buffers for rendering particles and the like
        _particle_material.SetBuffer(renderer_particle_buffer, _BM.PARTICLES_BUFFER);
        _particle_material.SetBuffer(renderer_float_buffer, _BM.PARTICLES_NEAR_PRESSURE_BUFFER);
        _particle_material.SetFloat(renderer_size, _render_size);
        _particle_material.SetInt(renderer_color_toggle, _color_toggle ? 1 : 0);
        _particle_material.SetFloat(renderer_max_float, max_color_val);

        Graphics.DrawMeshInstancedIndirect(
                    _particle_mesh, 
                    0, 
                    _particle_material, 
                    new Bounds(Vector3.zero, new Vector3(1000f, 1000f, 1000f)),
                    ARG_BUFFER, 
                    castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
                );
    }

    void OnDestroy() {
        // While buffer manager handles releasing buffers, we also have to release our own.
        if (predictedPositionsBuffer != null) predictedPositionsBuffer.Release();
        if (spatialIndices != null) spatialIndices.Release();
        if (spatialOffsets != null) spatialOffsets.Release();
    }

}
