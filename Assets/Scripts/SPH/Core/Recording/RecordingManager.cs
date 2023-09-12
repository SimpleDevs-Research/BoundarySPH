using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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
    [SerializeField] private RecStart _recordingStartSetting = RecStart.Off;
    [SerializeField] private float _delayBeforeStart = 0f;

    [Header("=== RECORDING STATISTICS ===")]
    private float _timeStarted;
    private float _recordingTimeStarted;
    [SerializeField, ReadOnly] private float _timeElapsed;
    [SerializeField, ReadOnly] private bool _initialized = false;

    public void Initialize() {
        // Make sure our references are all set
        _timeStarted = Time.time;
        switch(_recordingStartSetting) {
            case RecStart.At_Start:
                StartRecording();
                break;
            case RecStart.At_Particle_Controller_Start:
                _delayBeforeStart = _PC.startDelay;
                break;
            case RecStart.Custom_Time:
                break;
            default:
                break;
        }
    }

    public void StartRecording() {
        _recordingTimeStarted = Time.time;
        _initialized = true;
    }

    void Update() {
        if (!_initialized) {
            if (_recordingStartSetting == RecStart.At_Particle_Controller_Start || _recordingStartSetting == RecStart.Custom_Time) {
                if (Time.time - _timeStarted >= _delayBeforeStart) StartRecording();
            }
            return;
        }
        _timeElapsed = Time.time - _recordingTimeStarted;

    }
}
