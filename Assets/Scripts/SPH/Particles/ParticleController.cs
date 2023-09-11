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

    [Header("== PARTICLE CONFIGURATIONS ==")]
    [SerializeField, Tooltip("How big (visually) should the particles be?")]
    private float _particleRenderSize = 0.5f;
    public float particleRenderSize => _particleRenderSize;
    public float particleRenderRadius => _particleRenderSize / 2f;
    [SerializeField, ReadOnly, Tooltip("The max number of particles possible within this controller's grid section")]
    private int _MAX_NUM_PARTICLES = 0;
    [SerializeField, Tooltip("The number of particles that need to be generated. If `sectionIndex` is defined, will default to using that section's # of particles instead. If `sectionIndex` is set to -1, will use this value by default")]
    private int _numParticles = 100;
    public int numParticles => _numParticles;
    [SerializeField, ReadOnly]
    private int[] _numParticlesPerAxis;
    [SerializeField, ReadOnly]
    private int _numParticlesPerGridCell;
    [SerializeField, ReadOnly]
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
    [SerializeField, Tooltip("`c` - the speed of sound (10^-4 m/s, to model low Reynolds number flows). Used as a part of the `equation of state` to relate density to pressure. This is derived from the method mentioned in [https://www.sciencedirect.com/science/article/pii/S0045793021003145#b19], but the speed of sound value is derived from Morris et al. [https://www.sciencedirect.com/science/article/pii/S0021999197957764]")]
    private float _c = 5.77f;
    public float c => _c;
    private float _time_elapsed = 0f;
    private int _frames_elapsed = 0;

    [SerializeField] private float _spawnDistanceBetweenParticles = 1.45f;
    public float spawnDistanceBetweenParticles => _spawnDistanceBetweenParticles;

    int size_property = Shader.PropertyToID("size");
    int num_particles_per_cell_property = Shader.PropertyToID("numParticlesPerCell");
    int render_touching_property = Shader.PropertyToID("renderTouching");
    int particle_buffer_property = Shader.PropertyToID("particle_buffer");
    int pressure_buffer_property = Shader.PropertyToID("pressure_buffer");
    int velocity_buffer_property = Shader.PropertyToID("velocity_buffer");
    int grid_cell_buffer_property = Shader.PropertyToID("grid_cell_buffer");
    int render_limits_property = Shader.PropertyToID("render_limits_buffer");
    int render_denom_property = Shader.PropertyToID("velocity_denom");

    public enum RenderType { Off, Particles, ParticlesTouchingObstacles, GridCells, Both }
    [Header("== DEBUG CONTROLS ==")]
    [SerializeField] private RenderType _renderType = RenderType.Particles;
    [SerializeField] private float[] _renderLimits;
    [SerializeField, Range(0f,30f)] private float _renderDenom = 10f;

    private float t;

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
    
    public enum RecordSettings { Off, Online }

    [Header("== RECORDING CONFIGURATIONS ==")]
    [SerializeField, Tooltip("Public toggle to determine if we should be recording.")]
    private RecordSettings _record_statistics = RecordSettings.Off;
    [SerializeField, Tooltip("If recording online, then this identifies if we got a response."), ReadOnly] 
    private bool _recording_verified = false;
    [SerializeField] private bool _record_on_start = true;
    [SerializeField,ReadOnly] private bool _recording_started = false;

    [SerializeField, Tooltip("Files are saved as JSONs. Do not include the file extension in the name.")] 
    private string _record_name = "SPH_Findings";
    [SerializeField, Tooltip("The base URL that we send position data to.")]
    private string _record_positions_url = "";
    [SerializeField, Tooltip("The base URL that we send velocity data to.")]
    private string _record_velocities_url = "";
    [SerializeField, Tooltip("The base URL that we send density data to.")]
    private string _record_densities_url = "";
    [SerializeField, Tooltip("The base URL that we send pressure data to.")]
    private string _record_pressures_url = "";
    [SerializeField, Tooltip("The base URL that we send fps data to.")]
    private string _record_fps_url = "";
    [SerializeField, Tooltip("The duraction (in seconds) of the recording session")] 
    private float _record_duration = 10f;
    [SerializeField, Tooltip("How many records per batch? If Set to <= 0, will default to # particles")]
    private int _num_records_per_batch = 1000;
    [SerializeField, Tooltip("How long between each batch web request? If set to any value <= 0, will default to 0.01")]
    private float _time_between_requests = 0.01f;

    // This stores all recordings for offline recording sessions
    //private Recording[] _offline_recordings;
    // This stores a queue for all updates that need to be passed along to our server
    private List<OnlineRecordingBatch> _positions_batch_queue = new List<OnlineRecordingBatch>();
    private List<OnlineRecordingBatch> _velocities_batch_queue = new List<OnlineRecordingBatch>();
    private List<OnlineRecordingBatch> _densities_batch_queue = new List<OnlineRecordingBatch>();
    private List<OnlineRecordingBatch> _pressures_batch_queue = new List<OnlineRecordingBatch>();
    private List<OnlineRecordingBatch> _fps_batch_queue = new List<OnlineRecordingBatch>();
    
    private OnlineRecordingBatch _positions_current_batch = null;
    private OnlineRecordingBatch _velocities_current_batch = null;
    private OnlineRecordingBatch _densities_current_batch = null;
    private OnlineRecordingBatch _pressures_current_batch = null;
    //private OnlineRecordingBatch _fps_current_batch = null;

    [SerializeField] private TextMeshProUGUI _positions_request_textbox = null;
    [SerializeField] private TextMeshProUGUI _velocities_request_textbox = null;
    [SerializeField] private TextMeshProUGUI _densities_request_textbox = null;
    [SerializeField] private TextMeshProUGUI _pressures_request_textbox = null;
    //[SerializeField] private TextMeshProUGUI _fps_request_textbox = null;

    /*
    [SerializeField] private bool _record_online = true;
    [SerializeField, Tooltip("If recording online, this will be the HTTP url. If not, it will be the directory that you want to save the files in")] 
    private string _record_directory;
    [SerializeField, Tooltip("Simply a check if you want to save records or not.")] 
    private bool _record_statistics;
    
    private bool _record => _record_name.Length > 0 && _record_duration > 0f && _record_statistics;
    
    */

    [System.Serializable]
    public class Recording {
        public int particle_id;
        public List<SVector3> positions;
        public List<SVector3> velocities;
        public List<float> densities;
        public List<float> pressures;
        public Recording(int particle_id) {
            this.particle_id = particle_id;
            this.positions = new List<SVector3>();
            this.velocities = new List<SVector3>();
            this.densities = new List<float>();
            this.pressures = new List<float>();
        }
    }
    [System.Serializable]
    public class OnlineRecording {
        public int particle_id;
        public float timestamp;
        public int frame;
        public List<float> value;
        public OnlineRecording(int particle_id, float timestamp, int frame) {
            this.particle_id = particle_id;
            this.timestamp = timestamp;
            this.frame = frame;
            this.value = new List<float>();
        }
        public OnlineRecording(int particle_id, float timestamp, int frame, SVector3 value) {
            this.particle_id = particle_id;
            this.timestamp = timestamp;
            this.frame = frame;
            this.value = new List<float>();
            this.value.Add(value.x);
            this.value.Add(value.y);
            this.value.Add(value.z);
        }
    }
    [System.Serializable]
    public class OnlineRecordingBatch {
        public List<OnlineRecording> batch;
        public OnlineRecordingBatch() {
            batch = new List<OnlineRecording>();
        }
    }
    [System.Serializable]
    public class RecordingSession {
        public string session_name;
        public string timestamp;
        public float dt;
        public int num_particles;
        public float render_radius;
        public float h;
        public float rest_density;
        public float mu;
        public float k;
        public string metric;

        public float duration;
        public int frames;
        public int success;

        public RecordingSession(
                string session_name,
                string timestamp,
                float dt, 
                int num_particles, 
                float render_radius, 
                float h, 
                float rest_density,
                float mu,
                float k,
                string metric
        ) {
            this.session_name = session_name;
            this.timestamp = timestamp;
            this.dt = dt;
            this.num_particles = num_particles;
            this.render_radius = render_radius;            
            this.h = h;
            this.rest_density = rest_density;
            this.mu = mu;
            this.k = k;
            this.metric = metric;
        }
    }

    private RecordingSession _positions_recording_session = null;
    private RecordingSession _velocities_recording_session = null;
    private RecordingSession _densities_recording_session = null;
    private RecordingSession _pressures_recording_session = null;
    //private RecordingSession _fps_recording_session = null;
    private string _recording_save_name;

    private bool _recording_saved = false;

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
        _spawnDistanceBetweenParticles = Mathf.Clamp(_spawnDistanceBetweenParticles,particleRenderRadius,_h);
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
    }

    private void PrepareRecording() {
        // Firstly, if the directory is empty, then we cannot record. So we default to turning off
        if (
            _record_positions_url.Length == 0
            && _record_velocities_url.Length == 0
            && _record_densities_url.Length == 0
            && _record_pressures_url.Length == 0
            && _record_fps_url.Length == 0
        ) {
            _record_statistics = RecordSettings.Off;
            _recording_verified = true;
            return;
        }

        // Prepare # records per batch, if using online
        if (_num_records_per_batch <= 0) _num_records_per_batch = _numParticles;
        // Prepare time interval between requests, if using online
        if (_time_between_requests <= 0f) _time_between_requests = 0.01f;
        // Get the current timestamp
        string timestamp = System.DateTime.Now.ToString("yyyy-MM-dd_HH-mm-ss");

        // We need to let the system know if we want to start the recording or not upon the startign of the simulation.
        _recording_started = _record_on_start;
        
        // Create our recording sessions
        string session_data;
        if (_record_positions_url.Length > 0) {
            _positions_recording_session = new RecordingSession(
                _record_name+"_positions",
                timestamp, _dt,
                _numParticles, particleRenderRadius, 
                _h, _rest_density, _mu, _k,
                "positions"
            );
            // Convert our recording session into a string
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_positions_recording_session);
            // Create a new batch, which will act as our initial batch
            _positions_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_positions_url+"create",session_data,false,OnlineRecordingCallback));
        }
        if (_record_velocities_url.Length > 0) {
            _velocities_recording_session = new RecordingSession(
                _record_name+"_velocities",
                timestamp, _dt,
                _numParticles, particleRenderRadius, 
                _h, _rest_density, _mu, _k,
                "velocities"
            );
            // Convert our recording session into a string
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_velocities_recording_session);
            // Create a new batch, which will act as our initial batch
            _velocities_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_velocities_url+"create",session_data,false,OnlineRecordingCallback));
        }
        if (_record_densities_url.Length > 0) {
            _densities_recording_session = new RecordingSession(
                _record_name+"_densities",
                timestamp, _dt,
                _numParticles, particleRenderRadius, 
                _h, _rest_density, _mu, _k,
                "densities"
            );
            // Convert our recording session into a string
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_densities_recording_session);
            // Create a new batch, which will act as our initial batch
            _densities_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_densities_url+"create",session_data,false,OnlineRecordingCallback));
        }
        if (_record_pressures_url.Length > 0) {
            _pressures_recording_session = new RecordingSession(
                _record_name+"_pressures",
                timestamp, _dt,
                _numParticles, particleRenderRadius, 
                _h, _rest_density, _mu, _k,
                "pressures"
            );
            // Convert our recording session into a string
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_pressures_recording_session);
            // Create a new batch, which will act as our initial batch
            _pressures_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_pressures_url+"create",session_data,false,OnlineRecordingCallback));
        }

        /*
        if (_record_fps_url.Length > 0) {
            _fps_recording_session = new RecordingSession(
                _record_name+"_fps",
                timestamp, _dt,
                1, particleRenderRadius, 
                _h, _rest_density, _mu, _k,
                "fps"
            );
            // Convert our recording session into a string
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_pressures_recording_session);
            // Create a new batch, which will act as our initial batch
            _pressures_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_pressures_url+"create",session_data,false,OnlineRecordingCallback));
        }
        */
        
        /*
        // If we're recording online, we need to contact the server to check if we can record
        if (_record_statistics == RecordSettings.Online) {            
            // Firstly, convert our recording session into a string
            string session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_recording_session);
            // Create a new batch, which will act as our initial batch
            _online_current_batch = new OnlineRecordingBatch();
            // Next, we initialize the recording by calling to our API
            // The callback function `OnlineRecordingCallback` determines if we record during the update loop
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_directory+"create",session_data,false,OnlineRecordingCallback));
        }
        */
        /*
        // If we're recording offline, we need to initialize folders and the like
        else {
            // We first need to determine the directory we'll be saving inside
            string recording_save_dir = SaveSystemMethods.GetSaveLoadDirectory(_record_directory);
            // We create that directory if it doesn't exist yet
            SaveSystemMethods.CheckOrCreateDirectory(recording_save_dir);
            // We then Check or Create the recording_save_dir + _record_name directory. Delete an existing one if needed
            _recording_save_name = (_record_name.EndsWith("/")) 
                ? recording_save_dir + timestamp + "_" + _record_name 
                : recording_save_dir + timestamp + "_" + _record_name + "/";
            // We have to delete the old directory, if one exists. We won't stand for overlapping recordings
            SaveSystemMethods.DeleteDirectory(_recording_save_name);
            // We create the save directory
            SaveSystemMethods.CheckOrCreateDirectory(_recording_save_name);
            // Create separate folder inside of this current session folder for all particles
            SaveSystemMethods.CheckOrCreateDirectory(_recording_save_name+"particles/");
            // We need to create recordings for each particle
            _offline_recordings = new Recording[_numParticles];
            for(int i = 0; i < _numParticles; i++) _offline_recordings[i] = new Recording(i);
            // Add our first entry into the recording
            Record();
            // Toggle our boolean for record verification
            _recording_verified = true;
        }
        */
    }

    public void StartRecording() {
        _recording_started = true;
    }

    private void OnlineRecordingCallback(WebRequests.Response response) {
        if (!response.success) {
            _record_statistics = RecordSettings.Off;
            Debug.LogError(response.response);
        }
        _recording_verified = true;
        // We start the coroutine for the batch processing
        if (_record_positions_url.Length > 0) StartCoroutine(PositionsBatchRequester());
        if (_record_velocities_url.Length > 0) StartCoroutine(VelocitiesBatchRequester());
        if (_record_densities_url.Length > 0) StartCoroutine(DensitiesBatchRequester());
        if (_record_pressures_url.Length > 0) StartCoroutine(PressuresBatchRequester());
        // Record our first entry
        if (_recording_started) Record();
    }

    private int _BLOCK_SIZE = 512;
    private int _NUM_BLOCKS_GRID;
    private int _NUM_BLOCKS_PARTICLES;
    [SerializeField, ReadOnly] private int _NUM_BOUNDARY_PARTICLES;
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
    private int _COMPUTE_INTERNAL_FORCES;
    private int _INTEGRATE, _INTEGRATE_DEBUG;
    private int _DAMPEN_BY_BOUNDS;
    private void InitializeKernels() {
        // Used on initialization. _CLEAR_GRID also is performed at the beginning of each update loop
        _CLEAR_GRID = _SPH_Shader.FindKernel("ClearGrid");
        _GENERATE_PARTICLES = _SPH_Shader.FindKernel("GenerateParticles");
        // These are run during each update loop
        _UPDATE_GRID = _SPH_Shader.FindKernel("UpdateGridCellCounts");
        _COMPUTE_DENSITY = _SPH_Shader.FindKernel("CV_ComputeDensity");
        _COMPUTE_INTERNAL_FORCES = _SPH_Shader.FindKernel("ComputeInternalForces");
        _INTEGRATE = _SPH_Shader.FindKernel("Integrate");
        _INTEGRATE_DEBUG = _SPH_Shader.FindKernel("Integrate_Debug");
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
        //_BM.PARTICLES_PRESSURES_BUFFER = new ComputeBuffer(_numParticles + _NUM_BOUNDARY_PARTICLES, sizeof(float));
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
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "pressure", _BM.PARTICLES_PRESSURES_BUFFER);
        _SPH_Shader.SetBuffer(_GENERATE_PARTICLES, "externalForces", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);

        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particleOffsets", PARTICLE_OFFSETS_BUFFER);
        _SPH_Shader.SetBuffer(_UPDATE_GRID, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);

        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "cellLimits", CELL_LIMITS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_DENSITY, "pressure", _BM.PARTICLES_PRESSURES_BUFFER);

        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "grid", _BM.PARTICLES_GRID_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "cellLimits", CELL_LIMITS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "particleNeighbors", PARTICLE_NEIGHBORS_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_COMPUTE_INTERNAL_FORCES, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);

        _SPH_Shader.SetBuffer(_INTEGRATE, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE, "externalForces", _BM.PARTICLES_EXTERNAL_FORCES_BUFFER);

        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "force", FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "pressureForces", _BM.PARTICLES_PRESSURE_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "viscosityForces", _BM.PARTICLES_VISCOSITY_FORCES_BUFFER);
        _SPH_Shader.SetBuffer(_INTEGRATE_DEBUG, "density", _BM.PARTICLES_DENSITIES_BUFFER);

        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "bounds", BOUNDS_BUFFER);
        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "particles", _BM.PARTICLES_BUFFER);
        _SPH_Shader.SetBuffer(_DAMPEN_BY_BOUNDS, "velocity", _BM.PARTICLES_VELOCITIES_BUFFER);
    }

    // Update is called once per frame
    void Update() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SPH_Shader == null) return;
        if (_record_statistics != RecordSettings.Off && !_recording_verified) return;
        // Update shader variables!
        UpdateShaderVariables(true);
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
        _SPH_Shader.Dispatch(_COMPUTE_INTERNAL_FORCES, Mathf.CeilToInt((float)(_numParticles+_NUM_BOUNDARY_PARTICLES) / 256f),1,1);

        //_SPH_Shader.Dispatch(compute_external_acceleration_kernel, Mathf.CeilToInt((float)numParticles / (float)_BLOCK_SIZE), 1, 1);

        // Integrate over particles, update their positions after taking all force calcualtions into account
        // Note: We ONLY update the positions and velocities of the active fluid particles, not any boundary particles we might encounter.
        _SPH_Shader.Dispatch(_INTEGRATE, Mathf.CeilToInt((float)_numParticles / 256f), 1, 1);
        // Make sure that particles are within bounds, limit them if so
        //_SPH_Shader.Dispatch(_DAMPEN_BY_BOUNDS, Mathf.CeilToInt((float)_numParticles / 256f), 1, 1);

        // If we're recording, record our session
        // Also note: if we're waiting for the recording to start, we won't actually record anything yet.
        if (_record_statistics != RecordSettings.Off && _recording_started) {
            // Update time elapsed
            _time_elapsed += _dt;
            _frames_elapsed += 1;
            // Still record if we haven't saved yet
            if (!_recording_saved && _time_elapsed <= _record_duration) Record();
            // If we've reached our time limit, do indicate that we've saved!
            if (!_recording_saved && (_time_elapsed >= _record_duration || _time_elapsed + _dt > _record_duration)) StartCoroutine(SaveRecording());
        }

        RENDER_LIMITS_BUFFER.SetData(_renderLimits);
        switch(_renderType) {
            case RenderType.Particles:
                _GRID.particle_material.SetFloat(size_property, particleRenderSize);
                _GRID.particle_material.SetFloat(num_particles_per_cell_property, (float)_numParticlesPerGridCell);
                _GRID.particle_material.SetFloat(render_denom_property, _renderDenom);
                _GRID.particle_material.SetInt(render_touching_property, 0);
                _GRID.particle_material.SetBuffer(pressure_buffer_property, _BM.PARTICLES_PRESSURES_BUFFER);
                _GRID.particle_material.SetBuffer(velocity_buffer_property, _BM.PARTICLES_VELOCITIES_BUFFER);
                _GRID.particle_material.SetBuffer(particle_buffer_property, _BM.PARTICLES_BUFFER);
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
            case RenderType.ParticlesTouchingObstacles:
                _GRID.particle_material.SetFloat(size_property, particleRenderSize);
                _GRID.particle_material.SetFloat(num_particles_per_cell_property, (float)_numParticlesPerGridCell);
                _GRID.particle_material.SetFloat(render_denom_property, _renderDenom);
                _GRID.particle_material.SetInt(render_touching_property, 1);
                _GRID.particle_material.SetBuffer(pressure_buffer_property, _BM.PARTICLES_PRESSURES_BUFFER);
                _GRID.particle_material.SetBuffer(particle_buffer_property, _BM.PARTICLES_BUFFER);
                _GRID.particle_material.SetBuffer(velocity_buffer_property, _BM.PARTICLES_VELOCITIES_BUFFER);
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
                _GRID.particle_material.SetBuffer(pressure_buffer_property, _BM.PARTICLES_PRESSURES_BUFFER);
                _GRID.particle_material.SetBuffer(particle_buffer_property, _BM.PARTICLES_BUFFER);
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
            default:
                return;
        }

    }

    private void Record() {
        if (_record_statistics == RecordSettings.Off || _recording_saved) return;

        // Current timestamp and frame are captured in `_time_elapsed` and `_frames_elapsed`

        // We're recording online. In this case, we're working with batches of data
        // The batch we're working with is `_online_current_batch`.
        // At specific time intervals, we pass on 
        //      `_positions/velocities/densities/pressures_current_batch` to our 
        //      `_positions/velocities/densities/pressures_batch_queue` List, 
        //      which is looked at by a separate coroutine

        // OnlineRecordingBatch is merely a container for a list of type OnlineRecording
        // Our server is expecting the data in the following format:
        // - particle_id : int
        // - timestamp : float
        // - frame : int
        // - value : SVector3

        // Depending on what data we're recording, we'll do each separately
        OP.Particle[] temp_particles = new OP.Particle[0];
        float3[] temp_velocities = new float3[0];
        float[] temp_densities = new float[0];
        float3[] temp_pressures = new float3[0];
        if (_record_positions_url.Length > 0) {
            temp_particles = new OP.Particle[_numParticles];
            _BM.PARTICLES_BUFFER.GetData(temp_particles);
        }
        if (_record_velocities_url.Length > 0) {
            temp_velocities = new float3[_numParticles];
            _BM.PARTICLES_VELOCITIES_BUFFER.GetData(temp_velocities);
        }
        if (_record_densities_url.Length > 0) {
            temp_densities = new float[_numParticles];
            _BM.PARTICLES_DENSITIES_BUFFER.GetData(temp_densities);
        }
        if (_record_pressures_url.Length > 0) {
            temp_pressures = new float3[_numParticles];
            _BM.PARTICLES_PRESSURE_FORCES_BUFFER.GetData(temp_pressures);
        }

        // Initialize the data package and particle
        OnlineRecording r;
        OP.Particle p;
        for(int i = 0; i < _numParticles; i++) {    
            // Do we need to update positions?
            if (_record_positions_url.Length > 0) {
                r = new OnlineRecording(i,_time_elapsed,_frames_elapsed);
                // We need to update positions batch
                p = temp_particles[i];
                r.value = new List<float>() { p.position[0], p.position[1], p.position[2] };
                _positions_current_batch.batch.Add(r);
                // At this moment, if the # of records in our batch exceeds _num_records_per_batch, we pass the batch and create a new one
                if (_positions_current_batch.batch.Count >= _num_records_per_batch) {
                    // Enqueue the batch for processing
                    _positions_batch_queue.Add(_positions_current_batch);
                    // create a new batch
                    _positions_current_batch = new OnlineRecordingBatch();
                }
            }

            // Do we need to update velocities?
            if (_record_velocities_url.Length > 0) {
                r = new OnlineRecording(i,_time_elapsed,_frames_elapsed);
                // We need to update velocities batch
                r.value = new List<float>() { temp_velocities[i][0], temp_velocities[i][1], temp_velocities[i][2] };
                _velocities_current_batch.batch.Add(r);
                // At this moment, if the # of records in our batch exceeds _num_records_per_batch, we pass the batch and create a new one
                if (_velocities_current_batch.batch.Count >= _num_records_per_batch) {
                    // Enqueue the batch for processing
                    _velocities_batch_queue.Add(_velocities_current_batch);
                    // create a new batch
                    _velocities_current_batch = new OnlineRecordingBatch();
                }
            }

            // Do we need to update densities
            if (_record_densities_url.Length > 0) {
                r = new OnlineRecording(i,_time_elapsed,_frames_elapsed);
                // We need to update densities batch
                r.value = new List<float>() { temp_densities[i] };
                _densities_current_batch.batch.Add(r);
                // At this moment, if the # of records in our batch exceeds _num_records_per_batch, we pass the batch and create a new one
                if (_densities_current_batch.batch.Count >= _num_records_per_batch) {
                    // Enqueue the batch for processing
                    _densities_batch_queue.Add(_densities_current_batch);
                    // create a new batch
                    _densities_current_batch = new OnlineRecordingBatch();
                }
            }

            // Do we need to update pressures
            if (_record_pressures_url.Length > 0) {
                r = new OnlineRecording(i,_time_elapsed,_frames_elapsed);
                // We need to update densities batch
                r.value = new List<float>() { temp_pressures[i][0], temp_pressures[i][1], temp_pressures[i][2] };
                _pressures_current_batch.batch.Add(r);
                // At this moment, if the # of records in our batch exceeds _num_records_per_batch, we pass the batch and create a new one
                if (_pressures_current_batch.batch.Count >= _num_records_per_batch) {
                    // Enqueue the batch for processing
                    _pressures_batch_queue.Add(_pressures_current_batch);
                    // create a new batch
                    _pressures_current_batch = new OnlineRecordingBatch();
                }
            }
        }

        /*
        // We must determine best course of action depending on what kind of recording we are using
        if (_record_statistics == RecordSettings.Online) {
            // We're recording online. In this case, we're working with batches of data
            // The batch we're working with is `_online_current_batch`.
            // At specific time intervals, we pass on `_online_current_batch` to our `_positions_batch_queue` Queue, which is looked at by a separate coroutine
            // OnlineRecordingBatch is merely a container for a list of type OnlineRecording
            // Our server is expecting the data in the following format:
            // - particle_id : int
            // - timestamp : float
            // - frame : int
            // - position : SVector3
            // - velocity : SVector3
            // - density : float
            // - pressure : SVector3
            
            // Initialize the data package
            OnlineRecording r;
            for(int i = 0; i < _numParticles; i++) {
                Particle p = temp_particles[i];
                r = new OnlineRecording(
                    i,
                    _time_elapsed,
                    _frames_elapsed,
                    p.position,
                    temp_velocities[i],
                    temp_densities[i],
                    temp_pressures[i]
                );
                // push our current record into the current `online_current_batch`
                _online_current_batch.batch.Add(r);
                // At this moment, if the # of records in our batch exceeds _num_records_per_batch, we pass the batch and create a new one
                if (_online_current_batch.batch.Count >= _num_records_per_batch) {
                    // Enqueue the batch for processing
                    _positions_batch_queue.Add(_online_current_batch);
                    // create a new batch
                    _online_current_batch = new OnlineRecordingBatch();
                }
            }
        } else {
            // We're recording offline. This means we have a different process more optimized for local handling
            Recording r;
            for(int i = 0; i < _numParticles; i++) {
                r = _offline_recordings[i];
                // Add to position
                r.positions.Add(new SVector3(temp_particles[i].position[0],temp_particles[i].position[1],temp_particles[i].position[2]));
                // Add to velocities
                r.velocities.Add(new SVector3(temp_velocities[i][0],temp_velocities[i][1],temp_velocities[i][2]));
                // Add to densities
                r.densities.Add(temp_densities[i]);
                // Add to pressures
                r.pressures.Add((new SVector3(temp_pressures[i][0], temp_pressures[i][1], temp_pressures[i][2])).magnitude);
                // Set it back inside _recordings
                _offline_recordings[i] = r;
            }
            Debug.Log("Saving Recording!");
        }    
        */
    }

    private IEnumerator PositionsBatchRequester() {
        while(_record_statistics == RecordSettings.Online) {
            if (_positions_batch_queue.Count == 0) {
                yield return null;
                continue;
            }
            // We grab the first item in the queue
            OnlineRecordingBatch b = _positions_batch_queue[0];
            _positions_batch_queue.Remove(b);
            
            if (_positions_request_textbox != null) 
                _positions_request_textbox.text = $"Positions: {_positions_batch_queue.Count} Request Remaining";
            // We build the JSON for this batch
            string b_data = SaveSystemMethods.ConvertToJSON<OnlineRecordingBatch>(b);
            // We pass it to our server
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_positions_url+"update_batch",b_data,false));
            // We wait for a specific time to prevent overloading the server
            yield return new WaitForSeconds(_time_between_requests);
        }
    }
    private IEnumerator VelocitiesBatchRequester() {
        while(_record_statistics == RecordSettings.Online) {
            if (_velocities_batch_queue.Count == 0) {
                yield return null;
                continue;
            }
            // We grab the first item in the queue
            OnlineRecordingBatch b = _velocities_batch_queue[0];
            _velocities_batch_queue.Remove(b);
            
            if (_velocities_request_textbox != null) 
                _velocities_request_textbox.text = $"Velocities: {_velocities_batch_queue.Count} Request Remaining";
            // We build the JSON for this batch
            string b_data = SaveSystemMethods.ConvertToJSON<OnlineRecordingBatch>(b);
            // We pass it to our server
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_velocities_url+"update_batch",b_data,false));
            // We wait for a specific time to prevent overloading the server
            yield return new WaitForSeconds(_time_between_requests);
        }
    }
    private IEnumerator DensitiesBatchRequester() {
        while(_record_statistics == RecordSettings.Online) {
            if (_densities_batch_queue.Count == 0) {
                yield return null;
                continue;
            }
            // We grab the first item in the queue
            OnlineRecordingBatch b = _densities_batch_queue[0];
            _densities_batch_queue.Remove(b);
            
            if (_densities_request_textbox != null) 
                _densities_request_textbox.text = $"Densities: {_densities_batch_queue.Count} Request Remaining";
            // We build the JSON for this batch
            string b_data = SaveSystemMethods.ConvertToJSON<OnlineRecordingBatch>(b);
            // We pass it to our server
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_densities_url+"update_batch",b_data,false));
            // We wait for a specific time to prevent overloading the server
            yield return new WaitForSeconds(_time_between_requests);
        }
    }
    private IEnumerator PressuresBatchRequester() {
        while(_record_statistics == RecordSettings.Online) {
            if (_pressures_batch_queue.Count == 0) {
                yield return null;
                continue;
            }
            // We grab the first item in the queue
            OnlineRecordingBatch b = _pressures_batch_queue[0];
            _pressures_batch_queue.Remove(b);
            
            if (_pressures_request_textbox != null) 
                _pressures_request_textbox.text = $"Pressures: {_pressures_batch_queue.Count} Request Remaining";
            // We build the JSON for this batch
            string b_data = SaveSystemMethods.ConvertToJSON<OnlineRecordingBatch>(b);
            // We pass it to our server
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_pressures_url+"update_batch",b_data,false));
            // We wait for a specific time to prevent overloading the server
            yield return new WaitForSeconds(_time_between_requests);
        }
    }

    private IEnumerator SaveRecording() {
        // Turn off saving
        _recording_saved = true;
        
        // Pass on the current batch to our requesters.
        string session_data;
        if (_record_positions_url.Length > 0) {
            _positions_batch_queue.Add(_positions_current_batch);
            // Save duration-related data
            _positions_recording_session.duration = _time_elapsed;
            _positions_recording_session.frames = _frames_elapsed;
            _positions_recording_session.success = 1;
            // Convert and save session to JSON
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_positions_recording_session);
            // Update our API
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_positions_url+"update_session",session_data,false,OnlineSaveCallback));
            yield return null;
        }
        if (_record_velocities_url.Length > 0) {
            _velocities_batch_queue.Add(_velocities_current_batch);
            _velocities_recording_session.duration = _time_elapsed;
            _velocities_recording_session.frames = _frames_elapsed;
            _velocities_recording_session.success = 1;
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_velocities_recording_session);
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_velocities_url+"update_session",session_data,false,OnlineSaveCallback));
            yield return null;
        }
        if (_record_densities_url.Length > 0) {
            _densities_batch_queue.Add(_densities_current_batch);
            _densities_recording_session.duration = _time_elapsed;
            _densities_recording_session.frames = _frames_elapsed;
            _densities_recording_session.success = 1;
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_densities_recording_session);
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_densities_url+"update_session",session_data,false,OnlineSaveCallback));
            yield return null;
        }
        if (_record_pressures_url.Length > 0) {
            _pressures_batch_queue.Add(_pressures_current_batch);
            _pressures_recording_session.duration = _time_elapsed;
            _pressures_recording_session.frames = _frames_elapsed;
            _pressures_recording_session.success = 1;
            session_data = SaveSystemMethods.ConvertToJSON<RecordingSession>(_pressures_recording_session);
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_pressures_url+"update_session",session_data,false,OnlineSaveCallback));
            yield return null;
        }

        // Depending on if we're offline or online, decide the next course of action
        /*
        if (_record_statistics == RecordSettings.Online) {
            StartCoroutine(WebRequests.PostRequestWithJSON(_record_directory+"update_session",session_data,false,OnlineSaveCallback));
            yield return null;
        } else {
            SaveSystemMethods.SaveJSON(_recording_save_name+"session", session_data);
            yield return null;
            // Save particle data
            string record_json_data;
            int save_count = 0;
            for(int i = 0; i < _numParticles; i++) {
                record_json_data = SaveSystemMethods.ConvertToJSON<Recording>(_offline_recordings[i]);
                SaveSystemMethods.SaveJSON($"{_recording_save_name}particles/{i}",record_json_data);
                save_count += 1;
                if (save_count >= 10000) {
                    save_count = 0;
                    yield return null;
                }
            }
            Debug.Log("Saving Particle Data Recording!");
            yield return null;
        }
        */
    }

    /*
    private void Record() {
        if (_recording_saved || _recording_session == null) return;

        
        
        Recording r;
        for(int i = 0; i < _numParticles; i++) {
            r = _recordings[i];
            // Add to position
            r.positions.Add(new SVector3(temp_particles[i].position[0],temp_particles[i].position[1],temp_particles[i].position[2]));
            // Add to velocities
            r.velocities.Add(new SVector3(temp_velocities[i][0],temp_velocities[i][1],temp_velocities[i][2]));
            // Add to densities
            r.densities.Add(temp_densities[i]);
            // Add to pressures
            r.pressures.Add((new SVector3(temp_pressures[i][0], temp_pressures[i][1], temp_pressures[i][2])).magnitude);
            // Set it back inside _recordings
            _recordings[i] = r;
        }

        // If the time elapsed is equal to or greater than _record_duration, then we save our recording session into a JSON
        if (_time_elapsed >= _record_duration || _time_elapsed + _dt > _record_duration) {
            // Prevent future saves
            _recording_saved = true;
            // Convert and save particle data via Coroutine
            StartCoroutine(SaveRecordingCoroutine());
        } else {
            Debug.Log("Saving Recording!");
        }
    }
    */

    void OnDestroy() {
        ARG_BUFFER.Release();
        SPAWN_BOUNDS_BUFFER.Release();
        BOUNDS_BUFFER.Release();
        CELL_LIMITS_BUFFER.Release();
        PARTICLE_NEIGHBORS_BUFFER.Release();
        PARTICLE_OFFSETS_BUFFER.Release();
        FORCES_BUFFER.Release();
        RENDER_LIMITS_BUFFER.Release();
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
