using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;

public class PressureRenderer : MonoBehaviour
{
    [Header("== REFERENCES ==")]
    [SerializeField, Tooltip("Reference to the singular BufferManager component that handles all buffers in the system")]
    public BufferManager _BM = null;
    [SerializeField, Tooltip("Reference to a ParticleGrid component that acts as this controller's grid system")]
    public ParticleGrid _GRID = null;
    [SerializeField, Tooltip("Index of the section within the ParticleGrid component that we want to spawn particles inside. If set to -1, then we will use the number defined in this component and not particleGrid's.")]
    private int _SECTION_INDEX = -1;
    [SerializeField, Tooltip("The compute shader that handles all of our GPU calculations. Must-have, otherwise the script will not run")]
    private ComputeShader _SHADER = null;
    [SerializeField, Tooltip("The ParticleController that controls the SPH flow")]
    private ParticleController _PC = null;

    [Header("== CONFIGURATIONS ==")]
    [SerializeField] private float[] _gridCellRenderLimits;
    private int3 _NUM_BLOCKS_GRID;
    private int _NUM_BLOCKS_PARTICLES;

    [Header("== DEBUG ==")]
    [SerializeField] private bool _getPressures = false;
    [SerializeField] private int _numTempPressures;
    [ReadOnly, SerializeField] private float[] _tempPressures;
    [ReadOnly, SerializeField] private OP.GridCell[] _tempCells;
    //[ReadOnly, SerializeField] private int[] _tempParticles;

    void Start() {
        if (_GRID == null || _SHADER == null || _PC == null) {
            Debug.LogError("Pressure Renderer - ERROR: Cannot operate if either `GRID`, `SHADER`, or `PC` is set to `null`. Please define these references and restart the simulation.");
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
        _SHADER.Dispatch(_CLEAR_GRID, _NUM_BLOCKS_GRID[0], _NUM_BLOCKS_GRID[1], _NUM_BLOCKS_GRID[2]);

        Debug.Log("Pressure Renderer: Initialized!");
    }

    private void InitializeVariables() {
        // Determine block number for GPU threading
        _NUM_BLOCKS_GRID = new(
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[0] / 8f),
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[1] / 8f),
            Mathf.CeilToInt((float)_GRID.numCellsPerAxis[2] / 8f)
        );
        _NUM_BLOCKS_PARTICLES = Mathf.CeilToInt((float)_PC.numParticles / 256f);
        if (_getPressures) {
            _tempPressures = new float[_numTempPressures];
            _tempCells = new OP.GridCell[_numTempPressures];
        }
        //_tempParticles = new int[_numTempPressures];
    }

    private void InitializeShaderVariables() {
        // == WORLD CONFIGURATIONS ==
        _SHADER.SetFloat("gridCellSize", _GRID.gridCellSize);
        _SHADER.SetInts("numCellsPerAxis", _GRID.numCellsPerAxis);
        _SHADER.SetInt("total_number_of_cells", _GRID.numGridCells);

        _SHADER.SetFloat("gridScalingX", _GRID.gridScaling[0]);
        _SHADER.SetFloat("gridScalingY", _GRID.gridScaling[1]);
        _SHADER.SetFloat("gridScalingZ", _GRID.gridScaling[2]);

        // == PARTICLE CONFIGURATIONS ==
        _SHADER.SetInt("numParticles", _PC.numParticles);

        // Update variables that may change over time due to modifying inspector values
        UpdateShaderVariables();
    }

    private void UpdateShaderVariables() {
        _gridCellRenderLimits = (_SECTION_INDEX == -1) 
            ? _GRID.outerBounds
            : _GRID.sections[_SECTION_INDEX].bounds;
    
        _SHADER.SetFloat("bulkModulus", _PC.k);
    }

    private int _CLEAR_GRID, _UPDATE_PRESSURES, _CONDENSE_PRESSURES;
    private int size_property, grid_cell_buffer_property, grid_cell_pressures_property, grid_cell_render_limits_property;
    private void InitializeKernels() {
        // Used on initialization. _CLEAR_GRID also is performed at the beginning of each update loop
        _CLEAR_GRID = _SHADER.FindKernel("ClearGrid");
        // These are run during each update loop
        _UPDATE_PRESSURES = _SHADER.FindKernel("UpdatePressures");
        _CONDENSE_PRESSURES = _SHADER.FindKernel("CondensePressures");

        size_property = Shader.PropertyToID("size");
        grid_cell_buffer_property = Shader.PropertyToID("grid_cell_buffer");
        grid_cell_pressures_property = Shader.PropertyToID("pressures_buffer");
        grid_cell_render_limits_property = Shader.PropertyToID("render_limits_buffer");
    }

    private ComputeBuffer ARG_BUFFER;
    private ComputeBuffer PRESSURE_GRID_BUFFER;
    private ComputeBuffer BOUNDS_BUFFER, GRID_RENDER_LIMITS_BUFFER;
    //private ComputeBuffer TEMP_PARTICLES_BUFFER;
    private void InitializeBuffers() {
        // Initialize the buffers
        uint[] arg = {_GRID.grid_cell_mesh.GetIndexCount(0), (uint)(_GRID.numGridCells), _GRID.grid_cell_mesh.GetIndexStart(0), _GRID.grid_cell_mesh.GetBaseVertex(0), 0};
        ARG_BUFFER = new ComputeBuffer(1, arg.Length * sizeof(uint), ComputeBufferType.IndirectArguments);
        ARG_BUFFER.SetData(arg);
        
        PRESSURE_GRID_BUFFER = new ComputeBuffer(_GRID.numGridCells, sizeof(int)*2 + sizeof(float)*3);
        BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
        BOUNDS_BUFFER.SetData(_GRID.outerBounds);
        GRID_RENDER_LIMITS_BUFFER = new ComputeBuffer(6, sizeof(float));
        GRID_RENDER_LIMITS_BUFFER.SetData(_gridCellRenderLimits);
        _BM.PARTICLES_PRESSURES_BUFFER = new ComputeBuffer(_GRID.numGridCells, sizeof(float));

        //TEMP_PARTICLES_BUFFER = new ComputeBuffer(_PC.numParticles, sizeof(int));

        // Setting the buffers
        _SHADER.SetBuffer(_CLEAR_GRID, "grid", PRESSURE_GRID_BUFFER);
        _SHADER.SetBuffer(_CLEAR_GRID, "pressures", _BM.PARTICLES_PRESSURES_BUFFER);
        _SHADER.SetBuffer(_CLEAR_GRID, "bounds", BOUNDS_BUFFER);

        _SHADER.SetBuffer(_UPDATE_PRESSURES, "grid", PRESSURE_GRID_BUFFER);
        _SHADER.SetBuffer(_UPDATE_PRESSURES, "bounds", BOUNDS_BUFFER);
        _SHADER.SetBuffer(_UPDATE_PRESSURES, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(_UPDATE_PRESSURES, "density", _BM.PARTICLES_DENSITIES_BUFFER);
        //_SHADER.SetBuffer(_UPDATE_PRESSURES, "tempParticles", TEMP_PARTICLES_BUFFER);

        _SHADER.SetBuffer(_CONDENSE_PRESSURES, "grid", PRESSURE_GRID_BUFFER);
        _SHADER.SetBuffer(_CONDENSE_PRESSURES, "pressures", _BM.PARTICLES_PRESSURES_BUFFER);

        _GRID.grid_cell_material.SetBuffer(grid_cell_buffer_property, PRESSURE_GRID_BUFFER);
        _GRID.grid_cell_material.SetBuffer(grid_cell_render_limits_property, GRID_RENDER_LIMITS_BUFFER);
        _GRID.grid_cell_material.SetBuffer(grid_cell_pressures_property, _BM.PARTICLES_PRESSURES_BUFFER);
    }

    void Update() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null || _SHADER == null || _PC == null) return;
        // Update shader variables!
        UpdateShaderVariables();
        // Update Particles!
        UpdatePressures();

        if (_getPressures) {
            _BM.PARTICLES_PRESSURES_BUFFER.GetData(_tempPressures);
            PRESSURE_GRID_BUFFER.GetData(_tempCells);
            //TEMP_PARTICLES_BUFFER.GetData(_tempParticles);
        }
    }

    private void UpdatePressures() {
        _SHADER.Dispatch(_CLEAR_GRID, _NUM_BLOCKS_GRID[0], _NUM_BLOCKS_GRID[1], _NUM_BLOCKS_GRID[2]);
        _SHADER.Dispatch(_UPDATE_PRESSURES, _NUM_BLOCKS_PARTICLES,1,1);
        _SHADER.Dispatch(_CONDENSE_PRESSURES, _NUM_BLOCKS_GRID[0], _NUM_BLOCKS_GRID[1], _NUM_BLOCKS_GRID[2]);

        
        _GRID.grid_cell_material.SetFloat(size_property, _GRID.gridCellSize);
        GRID_RENDER_LIMITS_BUFFER.SetData(_gridCellRenderLimits);
        Graphics.DrawMeshInstancedIndirect(
            _GRID.grid_cell_mesh, 
            0, 
            _GRID.grid_cell_material, 
            new Bounds(Vector3.zero, new Vector3(100f, 100f, 100f)),
            ARG_BUFFER, 
            castShadows: UnityEngine.Rendering.ShadowCastingMode.Off
        );
        
    }

    void OnDestroy() {
        ARG_BUFFER.Release();
        PRESSURE_GRID_BUFFER.Release();
        GRID_RENDER_LIMITS_BUFFER.Release();
        BOUNDS_BUFFER.Release();

        //TEMP_PARTICLES_BUFFER.Release();
    }
}
