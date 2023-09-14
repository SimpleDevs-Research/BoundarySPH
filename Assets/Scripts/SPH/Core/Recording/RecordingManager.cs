using System.IO;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using SaveMethods = Helpers.SaveSystemMethods;
using System.Linq;
using Unity.Mathematics;

public class RecordingManager : MonoBehaviour
{
    public enum RecStart {
        Off,
        At_Start,
        At_Particle_Controller_Start,
        Custom_Time,
        Manual
    }

    [Header("=== REFERENCES ===")]
    public BufferManager _BM = null;
    public ParticleController _PC = null;
    public ComputeShader _SHADER = null;

    [Header("=== RECORDING SETTINGS ===")]
    [SerializeField] private string _recordingName;
    [SerializeField] private string _recordingFilepath;
    [SerializeField] private RecStart _recordingStartSetting = RecStart.Off;
    [SerializeField] private float _delayBeforeStart = 0f;
    [SerializeField] private float _recordTimeBuffer = 0.1f;

    [Header("=== RECORDING STATISTICS ===")]
    [SerializeField, ReadOnly] private string _positions_filepath, _velocities_filepath;
    private float _timeStarted;
    private float _recordingTimeStarted;
    [SerializeField, ReadOnly] private float _timeElapsed;
    [SerializeField, ReadOnly] private int _numRecordsPassed = -1;
    [SerializeField] private float3[] particle_positions_array, particle_velocities_array;
    [SerializeField, ReadOnly] private bool _filespace_initialized = false;
    [SerializeField, ReadOnly] private bool _initialized = false;


    // === PRIVATE VARIABLES ===
    private int condenseParticleDataKernel;
    private ComputeBuffer particlePositionsHashedBuffer, particleVelocitiesHashedBuffer;
    private List<float3[]> position_raw_rows, velocity_raw_rows;
    private IEnumerator updateCoroutine;
    private StreamWriter positions_writer, velocities_writer;
 
    // This function must be called in order for the recording to actually start. It can be called via the editor in inspector, or called via `BufferManager`
    public void Initialize() {
        // Make sure our references are all set. Namely Shader and PC
        if (_BM == null || _PC == null || _SHADER == null) {
            Debug.LogError("[RecordingManager] Error: Cannot operate without reference to a BufferManager, ParticleController, and/or ComputeShader");
            return;
        }
        // Initialize the _timeStarted.
        // This variable does NOT mean the time the recording started. It indicates what the current time is at the start of the simulation.
        // The variable that tracks when the RECORDING starts is _recordingTimeStarted
        _timeStarted = Time.time;

        // We initialize our shader before we begin.
        InitializeShader();

        // We also need to initialize our file saving data
        InitializeFiles();

        // Depending on our starting settings, we will start recording:
        switch(_recordingStartSetting) {
            case RecStart.At_Start:
                // - at the beginning of the simulation, no holds barred
                StartRecording();
                break;
            case RecStart.At_Particle_Controller_Start:
                // - When the particle controller starts
                _delayBeforeStart = _PC.startDelay;
                break;
            case RecStart.Custom_Time:
                // - At a custom time set by the user
                break;
            default:
                // = The recorder is off and will not record unless "StartRecording()" is called externally.
                break;
        }
    }

    public void StartRecording() {
        _recordingTimeStarted = Time.time;
        _initialized = _filespace_initialized;
        if (_initialized) {
            updateCoroutine = FilesaveCoroutine();
            StartCoroutine(updateCoroutine);
        }
    }

    void Update() {
        // We have to check first if the recording is on or not.
        if (!_initialized) {
            // We only record if we are set to record when the particle controller starts or after a custom time.
            // Note that if we set "At_Start" at the settings, then this wouldn't be called
            // Also note that we won't do anything if htat setting was set to "Off"
            if (
                (_recordingStartSetting == RecStart.At_Particle_Controller_Start || _recordingStartSetting == RecStart.Custom_Time)
                && Time.time - _timeStarted >= _delayBeforeStart
            ) StartRecording();
            return;
        }

        // Update the elapsed time
        _timeElapsed = Time.time - _recordingTimeStarted;

        // We need to calculate the current time buffer we're on
        int curTimeBufferIndex;
        if (_recordTimeBuffer > 0f) {
            curTimeBufferIndex = Mathf.FloorToInt(_timeElapsed / _recordTimeBuffer);
        } else {
            curTimeBufferIndex = _numRecordsPassed + 1;
        }
        // We don't do anything if the current recorded time buffer index is the same as the current index
        if (curTimeBufferIndex <= _numRecordsPassed) return;

        // Assuming we got this far, we update _prevNumRecords
        _numRecordsPassed = curTimeBufferIndex;

        // We need to update our buffers
        _SHADER.Dispatch(condenseParticleDataKernel, Mathf.CeilToInt((float)_PC.numParticles / 64f), 1, 1);

        // We need to extract the data from our buffer
        particle_positions_array = new float3[_PC.numParticles];
        particle_velocities_array = new float3[_PC.numParticles];
        particlePositionsHashedBuffer.GetData(particle_positions_array);
        _BM.PARTICLES_VELOCITIES_BUFFER.GetData(particle_velocities_array);

        // push the current state of the array into our buffer space, that'll eventually be read by our coroutine
        position_raw_rows.Add(particle_positions_array);
        velocity_raw_rows.Add(particle_velocities_array);
    }

    private void InitializeShader() {
        InitializeShaderVariables();
        InitializeShaderKernels();
        InitializeShaderBuffers();
    }

    private void InitializeShaderVariables() {
        _SHADER.SetInt("numParticles", _PC.numParticles);
    }

    private void InitializeShaderKernels() {
        condenseParticleDataKernel = _SHADER.FindKernel("CondenseParticleData");
    }

    private void InitializeShaderBuffers() {
        // We need to initialize two buffers: one for particle positions, and another for particle velocities
        particlePositionsHashedBuffer = new ComputeBuffer(_PC.numParticles, sizeof(float)*3);
        //particleVelocitiesHashedBuffer = new ComputeBuffer(_PC.numParticles, sizeof(int));

        // Set buffers
        _SHADER.SetBuffer(condenseParticleDataKernel, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(condenseParticleDataKernel, "particle_velocities", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SHADER.SetBuffer(condenseParticleDataKernel, "particle_positions", particlePositionsHashedBuffer);
        //_SHADER.SetBuffer(condenseParticleDataKernel, "particle_positions_hashed", particlePositionsHashedBuffer);
        //_SHADER.SetBuffer(condenseParticleDataKernel, "particle_velocities_hashed", particleVelocitiesHashedBuffer);
    }

    private void InitializeFiles() {
        string finalFilepath = SaveMethods.GetSaveLoadDirectory(_recordingFilepath) + _recordingName + "/";
        if (SaveMethods.CheckOrCreateDirectory(finalFilepath)) {
            // At this point, we can initialize our writers for the particle positions and particle velocities
            // Initialize our filepaths
            _positions_filepath = finalFilepath + "positions.csv";
            _velocities_filepath = finalFilepath + "velocities.csv";
            // Initialize our header
            string[] headers = new string[_PC.numParticles + 1];
            headers[0] = "timestamp";
            for(int i = 1; i <= _PC.numParticles; i++) {
                headers[i] = (i-1).ToString();
            }
            string h = string.Join(",",headers);
            // Initialize our data lists
            position_raw_rows = new List<float3[]>();
            velocity_raw_rows = new List<float3[]>();
            // Initialize our StreamWriters
            positions_writer = new StreamWriter(_positions_filepath);
            velocities_writer = new StreamWriter(_velocities_filepath);
            // Write our headers into both writers
            positions_writer.WriteLine(h);
            velocities_writer.WriteLine(h);
            // Finally, indicate that we can write to files
            _filespace_initialized = true;
        }
    }

    private IEnumerator FilesaveCoroutine() {
        string newLine;
        // We can only do the stuff in this coroutine if we have initialized the filesave stuff.
        while(_filespace_initialized) {
            //string ct = _timeElapsed.ToString();
            string ct = _PC.dt_passed.ToString();
            if (position_raw_rows.Count > 0) {
                newLine = ct+"," + string.Join(", ", position_raw_rows[0].Select(i => i[0].ToString()+"|"+i[1]+"|"+i[2].ToString()).ToArray());
                //position_rows.Add(newLine);
                // Write to our streamwriter
                positions_writer.WriteLine(newLine);
                // Remove the first-most raw row
                position_raw_rows.RemoveAt(0);
            }
            if (velocity_raw_rows.Count > 0) {
                newLine = ct+"," + string.Join(", ", velocity_raw_rows[0].Select(i =>  i[0].ToString()+"|"+i[1]+"|"+i[2].ToString()).ToArray());
                //velocity_rows.Add(newLine);
                // Write to our streamwriter
                velocities_writer.WriteLine(newLine);
                // Remove the first-most raw row
                velocity_raw_rows.RemoveAt(0);
            }
            yield return null;
        }
        yield return null;
    }

    void OnDestroy() {
        if (updateCoroutine != null) {
            // We need to end the coroutine and therefore end the file writing, to prevent data corruption
            StopCoroutine(updateCoroutine);
            // We must flush and close opur writers
            positions_writer.Flush();
            positions_writer.Close();
            velocities_writer.Flush();
            velocities_writer.Close();
        }
        if (particlePositionsHashedBuffer != null) particlePositionsHashedBuffer.Release();
        if (particleVelocitiesHashedBuffer != null) particleVelocitiesHashedBuffer.Release();
    }
}
