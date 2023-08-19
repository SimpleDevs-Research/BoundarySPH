using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
using Random = UnityEngine.Random;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;


public class BoidsController : MonoBehaviour
{

    [Header("== REFERENCES ==")]
    [SerializeField, Tooltip("Reference to a ParticleGrid component that acts as this controller's grid system")]
    public ParticleGrid _GRID = null;
    [SerializeField, Tooltip("Index of the section within the ParticleGrid component that we want to spawn particles inside. If set to -1, then we will use the number defined in this component and not particleGrid's.")]
    private int _SECTION_INDEX = -1;
    [SerializeField, Tooltip("The compute shader that handles all of our GPU calculations.")]
    private ComputeShader _SHADER = null;
    [SerializeField, Tooltip("Prefab for the boid")]
    private MeshObs _boidPrefab;

    [Header("== BOID SETTINGS ==")]
    [SerializeField] private bool _autoGenerateBoids = true;
    [SerializeField] private List<MeshObsGPU.TestObstacle> _boids = new List<MeshObsGPU.TestObstacle>();
    public List<MeshObsGPU.TestObstacle> boids => _boids;
    public int numBoids => _boids.Count;
    //public float boidSize = 0.05f;
    //public Color boidColor = Color.red;
    [SerializeField, ReadOnly] private float _visualRange;
    [SerializeField, Tooltip("Make sure this is smaller than the `visualRange`.")] private float _innerRange;
    [SerializeField] private float _maxSpeed = 1.5f;
    [SerializeField] private float _minSpeed = 0.5f;
    [SerializeField] private float _cohesionFactor = 1f;
    [SerializeField] private float _separationFactor = 30f;
    [SerializeField] private float _alignmentFactor = 5f;
    [SerializeField] private float _externalForceFactor = 1f;
    [SerializeField] private float _turnSpeed = 2f;
    [SerializeField] private float _meshTranslateSpeed = 1f;
    [SerializeField] private float _meshTurnSpeed = 1f;
    [SerializeField] private float _dt = -1f;
    [SerializeField] private bool _useGravity = true;
    [SerializeField] private bool _restrictX = false, _restrictY = false, _restrictZ = false;

    [SerializeField] private Transform _targetPosition = null;
    [SerializeField, Range(0f,1f)] private float _targetBias = 0.5f;
    [SerializeField] private float _targetRange = 1f;

    [SerializeField, ReadOnly] private OP.Boid[] _gpuBoids;
    [SerializeField, ReadOnly] private float3[] _gpuBoidVelocities;
    [SerializeField, ReadOnly] private float3[] _gpuBoidCurrentDirections;
    [SerializeField, ReadOnly] private float[] _gpuBoidDirectionDiffs;

    private int _CPU_LIMIT = 2048;

    [Tooltip("Stores the number of boids in each grid cell")]
    private int[] grid;

    [Header("== DEBUG TOOLS ==")]
    [SerializeField] private bool _verbose = true;

    [SerializeField, ReadOnly] private bool _useGPU;
    [SerializeField, ReadOnly] private int _numGridBlocks;
    [SerializeField, ReadOnly] private int _numBoidBlocks;

    [SerializeField] private bool showRadii = true;
    [SerializeField] private bool showDirection = true;
    [SerializeField] private bool showExternalForces = true;
    private bool showGizmos => showRadii || showDirection || showExternalForces;

    void OnDrawGizmos() {
        if (!Application.isPlaying || !showGizmos) return;
        int3[] external_forces = new int3[externalForcesBuffer.count];
        externalForcesBuffer.GetData(external_forces);
        for(int i = 0; i < numBoids; i++) {
            if (showRadii) {
                Gizmos.color = new Vector4(1f,0f,0f,0.25f);
                Gizmos.DrawSphere(_boids[i].obstacle.position_transform.position, _visualRange);
                Gizmos.color = new Vector4(0f,0f,1f,0.5f);
                Gizmos.DrawSphere(_boids[i].obstacle.position_transform.position, _innerRange);
            }
            if (showDirection) {
                Handles.color = Color.red;
                Handles.DrawLine(_gpuBoids[i].position, _gpuBoids[i].position + _gpuBoidVelocities[i] * 5f, 3);
            }
            if (showExternalForces) {
                if (_boids[i].obstacleID == -1) continue;
                int3 extInt3 = external_forces[_boids[i].obstacleID];
                float3 external_force = new(
                    (float)extInt3[0] / 1024f,
                    (float)extInt3[1] / 1024f,
                    (float)extInt3[2] / 1024f
                );
                Handles.color = Color.blue;
                Handles.DrawLine(_gpuBoids[i].position,_gpuBoids[i].position + external_force * 5f * _dt, 3);
                Handles.color = Color.black;
                Handles.DrawLine(_gpuBoids[i].position,_gpuBoids[i].position + (float3)new(0f,-9.81f,0f) * 5f * _dt, 3);
            }
        }
    }

    void Awake() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `GRID` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }

        if (_SHADER == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `SHADER` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }
        // We prime the boids to have an obstacle ID of -1. This is the catch in the `UpdateBoids` gpu kernel.
        // This will also be updated later after `Start()` is called in `MeshObsGPU`, if it is active in the scene.
        foreach(MeshObsGPU.TestObstacle obs in _boids) obs.obstacleID = -1;
        // Initialize some key variables central to Boid behavior
        InitializeGlobalVariables();
        // Initialize the GPU-specific parameters and buffers
        InitializeShaderKernels();
        InitializeShaderVariables();
        InitializeShaderBuffers();
        _SHADER.Dispatch(clearGridKernel, _numGridBlocks, 1, 1);
        // If we need to auto-generate the boids, we do so via the boolean check. Generating boids is done via the GPU
        if (_autoGenerateBoids) {
            _SHADER.Dispatch(generateBoidsKernel, _numBoidBlocks, 1, 1);
            // Reposition the relevant transforms inside of the `boids` array
            GenerateBoidTransforms();
        } 
        // In this case, we already have the boids set already. All that's needed is to update the boid transforms in the buffer
        else {
            SetBoidTransforms();
        }
        
    }

    private void InitializeGlobalVariables() {
        _visualRange = _GRID.gridCellSize;      // Get the visual range of boids based on grid cell size
        /*
        // Limit our # of boids if only on the CPU
        if (!_useGPU && numBoids > _CPU_LIMIT) {
            if (_verbose) Debug.Log($"WARNING - `numBoids` over the CPU limit. Will reduce down to clamp to {_CPU_LIMIT} boids");
            numBoids = _CPU_LIMIT;
        }
        */
        // Calcualte grid blocks and boid blocks
        _numGridBlocks = Mathf.CeilToInt((float)_GRID.numGridCells / 64f);  // Initialize `numGridBlocks` for GPU's sake
        _numBoidBlocks = Mathf.CeilToInt((float)numBoids / 64f);           // Initialize `numBoidBlocks` for GPU's sake
        if (_verbose) Debug.Log($"Number of Grid Blocks: {_numGridBlocks}\tNumber of Boid Blocks: {_numBoidBlocks}");
    }

    private int clearGridKernel, generateBoidsKernel;
    private int updateGridCellCountsKernel;
    private int prefixSumKernel, sumBlocksKernel, addSumsKernel;
    private int rearrangeBoidsKernel, updateBoidsKernel;
    private void InitializeShaderKernels() {
        clearGridKernel = _SHADER.FindKernel("ClearGrid");        
        generateBoidsKernel = _SHADER.FindKernel("GenerateBoids");
        updateGridCellCountsKernel = _SHADER.FindKernel("UpdateGridCellCounts");
        prefixSumKernel = _SHADER.FindKernel("PrefixSum");
        sumBlocksKernel = _SHADER.FindKernel("SumBlocks");
        addSumsKernel = _SHADER.FindKernel("AddSums");
        rearrangeBoidsKernel = _SHADER.FindKernel("RearrangeBoids");
        updateBoidsKernel = _SHADER.FindKernel("UpdateBoids");
    }

    private void InitializeShaderVariables() {
        // Grids specific
        _SHADER.SetInt("numGridCells", _GRID.numGridCells);
        _SHADER.SetFloat("gridCellSize", _GRID.gridCellSize);
        _SHADER.SetFloat("gridScalingX", _GRID.gridScaling[0]);
        _SHADER.SetFloat("gridScalingY", _GRID.gridScaling[1]);
        _SHADER.SetFloat("gridScalingZ", _GRID.gridScaling[2]);
        _SHADER.SetInts("numCellsPerAxis", _GRID.numCellsPerAxis);
        _SHADER.SetInt("numGridBlocks", _numGridBlocks);

        // Boids specific
        _SHADER.SetInt("numBoids", numBoids);
        _SHADER.SetFloat("visualRange", _visualRange);
        _SHADER.SetInt("randomSeed", Random.Range(0,int.MaxValue));

        // Make the call to update the shader variables
        UpdateShaderVariables();
    }

    private void UpdateShaderVariables() {
        float deltaTime = (_dt < 0) ? Time.deltaTime : _dt;
        float innerRange = Mathf.Clamp(_innerRange, 0.001f, _visualRange);
        
        _SHADER.SetFloat("innerRange", innerRange);
        _SHADER.SetFloat("maxSpeed", _maxSpeed);
        _SHADER.SetFloat("minSpeed", _minSpeed);
        _SHADER.SetFloat("cohesionFactor", _cohesionFactor);
        _SHADER.SetFloat("separationFactor", _separationFactor);
        _SHADER.SetFloat("alignmentFactor", _alignmentFactor);
        _SHADER.SetFloat("externalForceFactor", _externalForceFactor);
        _SHADER.SetFloat("turnSpeed", _turnSpeed);
        _SHADER.SetFloat("dt", deltaTime);

        _SHADER.SetInt("useGravity", (_useGravity) ? 1 : 0);
        _SHADER.SetInt("restrictX", (_restrictX) ? 1 : 0);
        _SHADER.SetInt("restrictY", (_restrictY) ? 1 : 0);
        _SHADER.SetInt("restrictZ", (_restrictZ) ? 1 : 0);

        _SHADER.SetFloat("targetBias", (_targetPosition == null) ? 0f : _targetBias);
        _SHADER.SetFloat("targetRange", (_targetPosition == null) ? 0f : _targetRange);
    }

    public ComputeBuffer gridBuffer, innerBoundsBuffer;
    private ComputeBuffer gridOffsetsBuffer, gridSumsBuffer1, gridSumsBuffer2;
    
    public ComputeBuffer boidsBuffer, boidVelocitiesBuffer;
    private ComputeBuffer boidOffsetsBuffer, rearrangedBoidsBuffer;

    private ComputeBuffer boidCurrentDirectionBuffer;
    private ComputeBuffer boidTargetPositionBuffer;

    private ComputeBuffer externalForcesBuffer;
    private bool externalForcesSet = false;

    private void InitializeShaderBuffers() {
        // === CREATING COMPUTE BUFFERS === //
        // Grid-specific buffers
            gridBuffer = new ComputeBuffer(_GRID.numGridCells, sizeof(int));
            innerBoundsBuffer = new ComputeBuffer(6,sizeof(float));
            gridOffsetsBuffer = new ComputeBuffer(_GRID.numGridCells, sizeof(int));
            gridSumsBuffer1 = new ComputeBuffer(_numGridBlocks, sizeof(int));
            gridSumsBuffer2 = new ComputeBuffer(_numGridBlocks, sizeof(int));
        // Boids-specific buffers
            boidsBuffer = new ComputeBuffer(numBoids, sizeof(float)*3 + sizeof(int)*5);
            boidVelocitiesBuffer = new ComputeBuffer(numBoids, sizeof(float)*3);
            _gpuBoids = new OP.Boid[numBoids];
            _gpuBoidVelocities = new float3[numBoids];
            boidOffsetsBuffer = new ComputeBuffer(numBoids, sizeof(int));
            rearrangedBoidsBuffer = new ComputeBuffer(numBoids, sizeof(int));
            boidCurrentDirectionBuffer = new ComputeBuffer(numBoids, sizeof(float)*3);
            _gpuBoidCurrentDirections = new float3[numBoids];
            _gpuBoidDirectionDiffs = new float[numBoids];
            boidTargetPositionBuffer = new ComputeBuffer(1, sizeof(float)*3);
            float3[] tempTargetPositionArray = new float3[1]; 
            tempTargetPositionArray[0] = (_targetPosition == null) 
                ? new(0f,0f,0f)
                : new(_targetPosition.position.x, _targetPosition.position.y, _targetPosition.position.z);
            boidTargetPositionBuffer.SetData(tempTargetPositionArray);
            
        // Obstacle Interaction buffers
            externalForcesBuffer = new ComputeBuffer(numBoids, sizeof(int)*3);

        // POPULATING COMPUTE BUFFERS
        // Grid-specific buffers
            if (_SECTION_INDEX >= 0)    innerBoundsBuffer.SetData(_GRID.sections[_SECTION_INDEX].bounds);
            else                        innerBoundsBuffer.SetData(_GRID.innerBounds);
            //innerBoundsBuffer.SetData(_GRID.innerBounds);

        // SETTING COMPUTE BUFFERS
        // Clear Grid
            _SHADER.SetBuffer(clearGridKernel, "grid", gridBuffer);
        // Generate Boids
            _SHADER.SetBuffer(generateBoidsKernel, "boids", boidsBuffer);
            _SHADER.SetBuffer(generateBoidsKernel, "boidVelocities", boidVelocitiesBuffer);
            _SHADER.SetBuffer(generateBoidsKernel, "innerBounds", innerBoundsBuffer);
        // Update Grid Cell Counts
            _SHADER.SetBuffer(updateGridCellCountsKernel, "boids", boidsBuffer);
            _SHADER.SetBuffer(updateGridCellCountsKernel, "grid", gridBuffer);
            _SHADER.SetBuffer(updateGridCellCountsKernel, "boidOffsets", boidOffsetsBuffer);
        // Prefix Sum
            _SHADER.SetBuffer(prefixSumKernel, "gridOffsetsIn", gridBuffer);
            _SHADER.SetBuffer(prefixSumKernel, "gridOffsets", gridOffsetsBuffer);
            _SHADER.SetBuffer(prefixSumKernel, "gridSumsBuffer", gridSumsBuffer2);
        // Sums 
            // All buffers are declared in the Update loop - no need to declare them here
        // Add Sums
            _SHADER.SetBuffer(addSumsKernel, "gridOffsets", gridOffsetsBuffer);
            // `GridSumsBufferIn` is declared in the upate loop - no need to decalre it here
        // Rearranging Boids
            _SHADER.SetBuffer(rearrangeBoidsKernel, "boids", boidsBuffer);
            _SHADER.SetBuffer(rearrangeBoidsKernel, "rearrangedBoids", rearrangedBoidsBuffer);
            _SHADER.SetBuffer(rearrangeBoidsKernel, "gridOffsets", gridOffsetsBuffer);
            _SHADER.SetBuffer(rearrangeBoidsKernel, "boidOffsets", boidOffsetsBuffer);
        // Update boids
            _SHADER.SetBuffer(updateBoidsKernel, "rearrangedBoids", rearrangedBoidsBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "boids", boidsBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "boidVelocities", boidVelocitiesBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "gridOffsets", gridOffsetsBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "innerBounds", innerBoundsBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "boidCurrentDirections", boidCurrentDirectionBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "externalForces", externalForcesBuffer);
            _SHADER.SetBuffer(updateBoidsKernel, "targetPos", boidTargetPositionBuffer);
    }

    public void SetExternalForcesBuffer(ComputeBuffer newBuffer) {
        externalForcesBuffer.Release();
        externalForcesBuffer = newBuffer;
        _SHADER.SetBuffer(updateBoidsKernel, "externalForces", externalForcesBuffer);
        boidsBuffer.GetData(_gpuBoids);
        for(int i = 0; i < numBoids; i++) _gpuBoids[i].obstacleID = _boids[i].obstacleID;
        boidsBuffer.SetData(_gpuBoids);
        externalForcesSet = true;
    }

    // Update is called once per frame
    void Update() {
        // We can't do anything if `grid` is null or if our compute shader is null
        if (_GRID == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `GRID` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }
        if (_SHADER == null) {
            Debug.LogError("BOIDS - ERROR: Cannot operate if `SHADER` is set to `null`. Please define this reference and restart the simulation.");
            return;
        }

        // Magic happens here
        // 1. Update the shader variables
        UpdateShaderVariables();
        // 2. Clean up `grid` so that all values = 0
        _SHADER.Dispatch(clearGridKernel, _numGridBlocks, 1, 1);
        // 3. We should update each boid with their current indice and projected index, while also declaring offset and performing atomic addition
        _SHADER.Dispatch(updateGridCellCountsKernel, _numBoidBlocks, 1, 1);
        // 4. We perform Prefix Summation
        _SHADER.Dispatch(prefixSumKernel, _numGridBlocks, 1, 1);
        bool swap = false;
        for (int d = 1; d < _numGridBlocks; d *= 2) {
            _SHADER.SetBuffer(sumBlocksKernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
            _SHADER.SetBuffer(sumBlocksKernel, "gridSumsBuffer", swap ? gridSumsBuffer2 : gridSumsBuffer1);
            _SHADER.SetInt("d", d);
            _SHADER.Dispatch(sumBlocksKernel, Mathf.CeilToInt((float)_numGridBlocks / 64f), 1, 1);
            swap = !swap;
        }
        _SHADER.SetBuffer(addSumsKernel, "gridSumsBufferIn", swap ? gridSumsBuffer1 : gridSumsBuffer2);
        _SHADER.Dispatch(addSumsKernel, _numGridBlocks, 1, 1);
        // 5. Rearrange boids based on prefix summation
        _SHADER.Dispatch(rearrangeBoidsKernel, _numBoidBlocks, 1, 1);
        // 6. Update the boids
        _SHADER.Dispatch(updateBoidsKernel, _numBoidBlocks, 1, 1);
        // 7. Update the boid transforms
        UpdateBoidTransforms();
    }

    private void GenerateBoidTransforms() {
        boidsBuffer.GetData(_gpuBoids);
        boidVelocitiesBuffer.GetData(_gpuBoidVelocities);
        for(int i = 0; i < numBoids; i++) {
            float3 v = _gpuBoidVelocities[i];
            MeshObs newBoid = Instantiate(
                _boidPrefab, 
                _gpuBoids[i].position, 
                Quaternion.LookRotation(new Vector3(v[0],v[1],v[2])),
                this.transform
            ) as MeshObs;
            _boids[i].obstacle = newBoid;
            _gpuBoidCurrentDirections[i] = new(
                newBoid.position_transform.forward.x, 
                newBoid.position_transform.forward.y, 
                newBoid.position_transform.forward.z
            );
        }
        boidCurrentDirectionBuffer.SetData(_gpuBoidCurrentDirections);
    }

    private void SetBoidTransforms() {
        boidsBuffer.GetData(_gpuBoids);
        boidVelocitiesBuffer.GetData(_gpuBoidVelocities);
        for(int i = 0; i < numBoids; i++) {
            _gpuBoids[i] = new OP.Boid();
            _gpuBoids[i].obstacleID = _boids[i].obstacleID;
            _gpuBoids[i].position = new(
                _boids[i].obstacle.position_transform.position.x,
                _boids[i].obstacle.position_transform.position.y,
                _boids[i].obstacle.position_transform.position.z
            );
            _gpuBoidVelocities[i] = (float3)new(
                _boids[i].obstacle.position_transform.forward.x, 
                _boids[i].obstacle.position_transform.forward.y, 
                _boids[i].obstacle.position_transform.forward.z
            ) * _maxSpeed;
            _gpuBoidCurrentDirections[i] = new(
                _boids[i].obstacle.position_transform.forward.x, 
                _boids[i].obstacle.position_transform.forward.y, 
                _boids[i].obstacle.position_transform.forward.z
            );
        }
        boidsBuffer.SetData(_gpuBoids);
        boidVelocitiesBuffer.SetData(_gpuBoidVelocities);
        boidCurrentDirectionBuffer.SetData(_gpuBoidCurrentDirections);
    }

    private void UpdateBoidTransforms() {
        boidsBuffer.GetData(_gpuBoids);
        boidVelocitiesBuffer.GetData(_gpuBoidVelocities);
        //float posStep = _meshTranslateSpeed;
        float rotStep = _meshTurnSpeed;
        Vector3 targetPos, targetRotEuler;
        Quaternion targetRot;
        Vector3 localEuler;
        for(int i = 0; i < numBoids; i++) {
            float3 v = _gpuBoidVelocities[i];
            float3 p = _gpuBoids[i].position;

            targetPos = new Vector3(p[0],p[1],p[2]);
            targetRotEuler = new Vector3(v[0],v[1],v[2]);
            targetRot = Quaternion.LookRotation(targetRotEuler);

            _boids[i].obstacle.position_transform.position = targetPos;
            //_boids[i].obstacle.position = Vector3.MoveTowards(_boids[i].obstacle.position, targetPos, posStep);
            _boids[i].obstacle.position_transform.rotation = Quaternion.RotateTowards(_boids[i].obstacle.position_transform.rotation, targetRot, rotStep);
 
            _gpuBoidCurrentDirections[i] = new(
                _boids[i].obstacle.position_transform.forward.x, 
                _boids[i].obstacle.position_transform.forward.y, 
                _boids[i].obstacle.position_transform.forward.z
            );
        }
        boidCurrentDirectionBuffer.SetData(_gpuBoidCurrentDirections);
    }

    void OnDestroy() {
        gridBuffer.Release();
        innerBoundsBuffer.Release();
        gridOffsetsBuffer.Release();
        gridSumsBuffer1.Release();
        gridSumsBuffer2.Release();
        boidsBuffer.Release();
        boidVelocitiesBuffer.Release();
        boidOffsetsBuffer.Release();
        rearrangedBoidsBuffer.Release();
        boidCurrentDirectionBuffer.Release();
        boidTargetPositionBuffer.Release();
        externalForcesBuffer.Release();
    }
}
