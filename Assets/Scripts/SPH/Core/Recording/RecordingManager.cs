using System.IO;
using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using SaveMethods = Helpers.SaveSystemMethods;
using System.Linq;
using Unity.Mathematics;

[System.Serializable]
public class RecordingWriter {
    public StreamWriter writer;
    public ComputeBuffer bufferRef;
    public int numParticles;
    public RecordingWriter(string fp, ComputeBuffer buffer, int n) {
        writer = new StreamWriter(fp);
        bufferRef = buffer;
        numParticles = n;
    }

    public virtual void GetData() {}
    public virtual void WriteData(string ct) {}
    public virtual bool HasRawRows() { return true; }

    public void DestroyWriter() {
        writer.Flush();
        writer.Close();
    }
}

public class RecordingWriterFloat3 : RecordingWriter {
    public List<float3[]> raws;

    public RecordingWriterFloat3(string fp, ComputeBuffer buffer, int n, string header) : base(fp, buffer, n) {
        writer.WriteLine(header);
        raws = new List<float3[]>();
    }

    public override void GetData() {
        // Prep the new current data array
        float3[] current_data = new float3[numParticles];
        // Get the data from the buffer
        bufferRef.GetData(current_data);
        // Save the line
        raws.Add(current_data);
    }

    public override void WriteData(string ct) {
        if (raws.Count == 0) return;
        // Prepare the next line
        string newLine = ct+"," + string.Join(", ", raws[0].Select(i => i[0].ToString()+"|"+i[1]+"|"+i[2].ToString()).ToArray());
        // Write to our streamwriter
        writer.WriteLine(newLine);
        // Remove the first-most raw row
        raws.RemoveAt(0);
    }

    public override bool HasRawRows() { return raws.Count > 0; }
}

public class RecordingWriterFloat : RecordingWriter {
    public List<float[]> raws;

    public RecordingWriterFloat(string fp, ComputeBuffer buffer, int n, string header) : base(fp, buffer, n) {
        writer.WriteLine(header);
        raws = new List<float[]>();
    }

    public override void GetData() {
        // Prep the new current data array
        float[] current_data = new float[numParticles];
        // Get the data from the buffer
        bufferRef.GetData(current_data);
        // Save the line
        raws.Add(current_data);
    }

    public override void WriteData(string ct) {
        if (raws.Count == 0) return;
        // Prepare the next line
        string newLine = ct+"," + string.Join(", ", raws[0].Select(i =>  i.ToString()).ToArray());
        // Write to our streamwriter
        writer.WriteLine(newLine);
        // Remove the first-most raw row
        raws.RemoveAt(0);
    }

    public override bool HasRawRows() { return raws.Count > 0; }
}

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
    [SerializeField, ReadOnly] private string _positions_filepath, _velocities_filepath, _densities_filepath, _pressures_filepath;
    private float _timeStarted;
    private float _recordingTimeStarted;
    [SerializeField, ReadOnly] private float _timeElapsed;
    [SerializeField, ReadOnly] private int _numRecordsPassed = -1;
    [SerializeField, ReadOnly] private List<RecordingWriter> _writers;
    [SerializeField] private float3[] particle_positions_array, particle_velocities_array;
    [SerializeField] private float[] particle_densities_array, particle_pressures_array;
    [SerializeField, ReadOnly] private bool _filespace_initialized = false;
    [SerializeField, ReadOnly] private bool _initialized = false;


    // === PRIVATE VARIABLES ===
    private int condenseParticleDataKernel;
    private ComputeBuffer particlePositionsHashedBuffer, particleVelocitiesHashedBuffer;
    private List<float3[]> position_raw_rows, velocity_raw_rows;
    private List<float[]> density_raw_rows, pressure_raw_rows;
    private IEnumerator updateCoroutine;
    private StreamWriter positions_writer, velocities_writer, densities_writer, pressures_writer;
 
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
        foreach(RecordingWriter w in _writers) {
            w.GetData();
        }
        /*
        particle_positions_array = new float3[_PC.numParticles];
        particle_velocities_array = new float3[_PC.numParticles];
        particle_densities_array = new float[_PC.numParticles];
        particle_pressures_array = new float[_PC.numParticles];
        particlePositionsHashedBuffer.GetData(particle_positions_array);
        _BM.PARTICLES_VELOCITIES_BUFFER.GetData(particle_velocities_array);
        _BM.PARTICLES_DENSITIES_BUFFER.GetData(particle_densities_array);
        _BM.PARTICLES_PRESSURES_BUFFER.GetData(particle_pressures_array);

        // push the current state of the array into our buffer space, that'll eventually be read by our coroutine
        position_raw_rows.Add(particle_positions_array);
        velocity_raw_rows.Add(particle_velocities_array);
        density_raw_rows.Add(particle_densities_array);
        pressure_raw_rows.Add(particle_pressures_array);
        */
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
            _densities_filepath = finalFilepath + "densities.csv";
            _pressures_filepath = finalFilepath + "pressures.csv";
            // Initialize our header
            string[] headers = new string[_PC.numParticles + 1];
            headers[0] = "timestamp";
            for(int i = 1; i <= _PC.numParticles; i++) {
                headers[i] = (i-1).ToString();
            }
            string h = string.Join(",",headers);
            // Initialize our RecordingWriters
            RecordingWriterFloat3 positionsWriter = new RecordingWriterFloat3(_positions_filepath, particlePositionsHashedBuffer, _PC.numParticles, h);
            RecordingWriterFloat3 velocitiesWriter = new RecordingWriterFloat3(_velocities_filepath, _BM.PARTICLES_VELOCITIES_BUFFER, _PC.numParticles, h);
            RecordingWriterFloat densitiesWriter = new RecordingWriterFloat(_densities_filepath, _BM.PARTICLES_DENSITIES_BUFFER, _PC.numParticles, h);
            RecordingWriterFloat pressuresWriter = new RecordingWriterFloat(_pressures_filepath, _BM.PARTICLES_PRESSURES_BUFFER, _PC.numParticles, h);
            // Add our writers
            _writers.Add(positionsWriter);
            _writers.Add(velocitiesWriter);
            _writers.Add(densitiesWriter);
            _writers.Add(pressuresWriter);
            /*
            // Initialize our data lists
            position_raw_rows = new List<float3[]>();
            velocity_raw_rows = new List<float3[]>();
            density_raw_rows = new List<float[]>();
            pressure_raw_rows = new List<float[]>();
            // Initialize our StreamWriters
            positions_writer = new StreamWriter(_positions_filepath);
            velocities_writer = new StreamWriter(_velocities_filepath);
            densities_writer = new StreamWriter(_densities_filepath);
            pressures_writer = new StreamWriter(_pressures_filepath);
            // Write our headers into both writers
            positions_writer.WriteLine(h);
            velocities_writer.WriteLine(h);
            densities_writer.WriteLine(h);
            pressures_writer.WriteLine(h);
            */
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
            foreach(RecordingWriter w in _writers) {
                w.WriteData(ct);
            }
            /*
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
            if (density_raw_rows.Count > 0) {
                newLine = ct+"," + string.Join(", ", density_raw_rows[0].Select(i =>  i.ToString()).ToArray());
                //velocity_rows.Add(newLine);
                // Write to our streamwriter
                densities_writer.WriteLine(newLine);
                // Remove the first-most raw row
                density_raw_rows.RemoveAt(0);
            }
            if (pressure_raw_rows.Count > 0) {
                newLine = ct+"," + string.Join(", ", pressure_raw_rows[0].Select(i =>  i.ToString()).ToArray());
                //velocity_rows.Add(newLine);
                // Write to our streamwriter
                pressures_writer.WriteLine(newLine);
                // Remove the first-most raw row
                pressure_raw_rows.RemoveAt(0);
            }
            */
            yield return null;
        }
        yield return null;
    }

    void OnDestroy() {
        if (updateCoroutine != null) {
            // We need to end the coroutine and therefore end the file writing, to prevent data corruption
            bool stillHasRaws = false;
            do {
                stillHasRaws = false;
                foreach(RecordingWriter w in _writers) {
                    stillHasRaws = stillHasRaws || w.HasRawRows();
                }
            } while (stillHasRaws);
            //while(position_raw_rows.Count > 0 || velocity_raw_rows.Count > 0 || density_raw_rows.Count > 0 || pressure_raw_rows.Count > 0) {}
            StopCoroutine(updateCoroutine);
            // We must flush and close opur writers
            foreach(RecordingWriter w in _writers) {
                w.DestroyWriter();
            }
            /*
            positions_writer.Flush();
            positions_writer.Close();
            velocities_writer.Flush();
            velocities_writer.Close();
            densities_writer.Flush();
            densities_writer.Close();
            pressures_writer.Flush();
            pressures_writer.Close();
            */
        }
        if (particlePositionsHashedBuffer != null) particlePositionsHashedBuffer.Release();
        if (particleVelocitiesHashedBuffer != null) particleVelocitiesHashedBuffer.Release();
    }
}
