using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BoidsController : MonoBehaviour
{
    [Header("== REFERENCES ==")]
    [SerializeField, Tooltip("Reference to a ParticleGrid component that acts as this controller's grid system")]
    public ParticleGrid _GRID = null;
    [SerializeField, Tooltip("Index of the section within the ParticleGrid component that we want to spawn particles inside. If set to -1, then we will use the number defined in this component and not particleGrid's.")]
    private int _SECTION_INDEX = -1;
    [SerializeField, Tooltip("The compute shader that handles all of our GPU calculations.")]
    private ComputeShader _Shader = null;

    [Header("== BOIDS SETTINGS ==")]
    [SerializeField] private int _numBoids;
    public int numBoids => _numBoids;
    private int _CPU_LIMIT = 2048;

    [Header("== DEBUG TOOLS ==")]
    [SerializeField] private bool _verbose = true;

    [SerializeField, ReadOnly] private bool _useGPU;
    [SerializeField, ReadOnly] private float _visualRange;
    [SerializeField, ReadOnly] private int _numBlocks;

    void Awake() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `GRID` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }


        // Initialize some key variables central to Boid behavior
        InitializeVariables();
        
        // If we're using the GPU, 
        if (_useGPU) {
        }

    }

    private void InitializeVariables() {
        _useGPU = _Shader != null;              // Should we use the GPU?
        _visualRange = _GRID.gridCellSize;      // Get the visual range of boids based on grid cell size
        _numBlocks = Mathf.CeilToInt((float)_GRID.numGridCells / 64f);  // Initialize `numBlocks` for GPU's sake
        if (_verbose) Debug.Log($"Number of blocks: {_numBlocks}");

        // Limit our # of boids if only on the CPU
        if (!_useGPU && _numBoids > _CPU_LIMIT) {
            if (_verbose) Debug.Log($"WARNING - `numBoids` over the CPU limit. Will reduce down to clamp to {_CPU_LIMIT} boids");
            _numBoids = _CPU_LIMIT;
        }

        

    }

    // Update is called once per frame
    void Update() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `GRID` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }
        
    }
}
