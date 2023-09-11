using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;

public class ContinuousFlow : MonoBehaviour
{
    public bool TRANSPOSING = true;

    [Header("=== REFERENCES ===")]
    public BufferManager _BM = null;
    public ParticleController _PC = null;
    public ComputeShader _SHADER;
    public Grid _GRID;
    public int inflowSectionIndex = -1, outflowSectionIndex = -1;
    
    public ComputeBuffer 
        particleBuffer, 
        inflowBoundsBuffer, outflowBoundsBuffer, boundCountsBuffer, inflowPositionsBuffer,
        particlesInOutflowBuffer, particleOffsetsBuffer, particleSumsBuffer1, particleSumsBuffer2,
        rearrangeParticlesBuffer, particlesToTransposeBuffer;
    
    private int calculateInflowPositionsKernel, clearRearrangedParticlesKernel,
        updateCountsKernel, prefixSumKernel, sumBlocksKernel, addSumsKernel, rearrangeParticlesKernel,
        transposeParticlesKernel;
    private int _numParticleBlocks;

    [SerializeField, ReadOnly] private int[] _numParticlesPerAxisInInflow;
    [SerializeField, ReadOnly] private int _numThreshold; 

    [SerializeField] private float[] _timePhases = new float[2]{15f,5f};

    [SerializeField, ReadOnly] private bool _initialized =false;

    [Header("=== DEBUG ===")]
    [SerializeField] private uint2[] particles_in_outflow_array;
    [SerializeField] private uint[] rearranged_particles_array;
    [SerializeField] private int[] bound_counts_array;

    void Awake() {
        if (_BM == null || _PC == null || _SHADER == null || _GRID == null || inflowSectionIndex == -1 || outflowSectionIndex == -1) return;
        // We need to grab our kernels
        calculateInflowPositionsKernel = _SHADER.FindKernel("CalculateInflowPositions");
        clearRearrangedParticlesKernel = _SHADER.FindKernel("ClearRearrangedParticles");
        updateCountsKernel = _SHADER.FindKernel("UpdateCounts");
        //prefixSumKernel = _SHADER.FindKernel("PrefixSum");
        //sumBlocksKernel = _SHADER.FindKernel("SumBlocks");
        //addSumsKernel = _SHADER.FindKernel("AddSums");
        //rearrangeParticlesKernel = _SHADER.FindKernel("RearrangeParticles");
        transposeParticlesKernel = _SHADER.FindKernel("TransposeParticles");
    }

    public void Initialize() {
        if (_BM == null || _PC == null || _SHADER == null || _GRID == null || inflowSectionIndex == -1 || outflowSectionIndex == -1) return;
        // Prep some variables
        _numParticleBlocks = Mathf.CeilToInt((float)_PC.numParticles / 64f);
        
        // Prep the inflow and outflow bounds arrays
        float[] inflowBoundsArray = _GRID.sections[inflowSectionIndex].bounds;
        float[] outflowBoundsArray = _GRID.sections[outflowSectionIndex].bounds;

        // Update our buffers for the boundaries
        inflowBoundsBuffer = new ComputeBuffer(6, sizeof(float));
        outflowBoundsBuffer = new ComputeBuffer(6, sizeof(float));
        inflowBoundsBuffer.SetData(inflowBoundsArray);
        outflowBoundsBuffer.SetData(outflowBoundsArray);

        // Update our buffers for the boundary counts
        boundCountsBuffer = new ComputeBuffer(2, sizeof(int));
        bound_counts_array = new int[2];

        // Update our buffer for the particles in outflow
        particlesInOutflowBuffer = new ComputeBuffer(_PC.numParticles, sizeof(uint)*2);
        particlesToTransposeBuffer = new ComputeBuffer(_PC.numParticles, sizeof(uint)*2);
        particles_in_outflow_array = new uint2[_PC.numParticles];
        for(int i = 0; i < _PC.numParticles; i++) {
            particles_in_outflow_array[i] = new((uint)i, (uint)0);
        }
        particlesInOutflowBuffer.SetData(particles_in_outflow_array);
        particlesToTransposeBuffer.SetData(particles_in_outflow_array);

        // Update the prefix sum-specific buffers
        particleOffsetsBuffer = new ComputeBuffer(_PC.numParticles, sizeof(int));
        particleSumsBuffer1 = new ComputeBuffer(_numParticleBlocks, sizeof(int));
        particleSumsBuffer2 = new ComputeBuffer(_numParticleBlocks, sizeof(int));

        // Update the rearranged particles buffer
        rearrangeParticlesBuffer = new ComputeBuffer(_PC.numParticles, sizeof(uint));
        rearranged_particles_array = new uint[_PC.numParticles];

        // We need to calculate the positions where we can put the particles, as well as the particle # threshold for transposing
        Vector3 bs = new Vector3(
            inflowBoundsArray[3]-inflowBoundsArray[0],
            inflowBoundsArray[4]-inflowBoundsArray[1],
            inflowBoundsArray[5]-inflowBoundsArray[2]
        );
        _numParticlesPerAxisInInflow = new int[3] {
            Mathf.FloorToInt(bs.x / _PC.spawnDistanceBetweenParticles),
            Mathf.FloorToInt(bs.y / _PC.spawnDistanceBetweenParticles),
            Mathf.FloorToInt(bs.z / _PC.spawnDistanceBetweenParticles)
        };
        _numThreshold = Mathf.Clamp(
            _numParticlesPerAxisInInflow[0] * _numParticlesPerAxisInInflow[1] * _numParticlesPerAxisInInflow[2], 
            0, _PC.numParticles
        );
        // Set up the buffers to calculate inflow positions
        inflowPositionsBuffer = new ComputeBuffer(_numThreshold, sizeof(float)*3);

        // Determine direction
        float3 g = new(_PC.g[0], _PC.g[1], _PC.g[2]);
        bool hasFlow = Unity.Mathematics.math.length(g) > 0f;
        float3 direction = (hasFlow) 
            ? Unity.Mathematics.math.normalize(g) 
            : new(1f,0f,0f);
        // Based on this flow direction, we need to calculate the "center" of the endmost possible bound space, starting from the origin of the simulation grid
        Vector3 outflowUpperBoundPos = _GRID.sections[outflowSectionIndex].UPPER_BOUND_REF.position;
        Vector3 dirToUpperBoundPos = outflowUpperBoundPos-_GRID.origin;
        Vector3 directionV3 = new Vector3(direction[0],direction[1],direction[2]);
        Vector3 componentAlongFlowDirToOutflowCenter = directionV3 * Vector3.Dot(dirToUpperBoundPos,directionV3);
        Vector3 outflowCenter = _GRID.origin + componentAlongFlowDirToOutflowCenter;

        // We need to establish some values and buffers
        // 0. Update variables
        _SHADER.SetInt("numParticles",_PC.numParticles);
        _SHADER.SetInt("numParticleBlocks",_numParticleBlocks);
        _SHADER.SetInt("numThreshold", _numThreshold);
        _SHADER.SetInts("numParticlesPerAxis", _numParticlesPerAxisInInflow);
        _SHADER.SetFloat("distanceBetweenParticles",_PC.spawnDistanceBetweenParticles);
        _SHADER.SetFloats("flowDirection", new float[3]{direction[0],direction[1],direction[2]});
        _SHADER.SetFloats("outflowCenter", new float[3]{outflowCenter.x, outflowCenter.y, outflowCenter.z});
        _SHADER.SetInt("defaultLength", (int)(componentAlongFlowDirToOutflowCenter.magnitude * 2f * 1024f));
        // 1. Inflow Positions
        _SHADER.SetBuffer(calculateInflowPositionsKernel, "inflow_positions", inflowPositionsBuffer);
        _SHADER.SetBuffer(calculateInflowPositionsKernel, "inflowBounds", inflowBoundsBuffer);
        // 2. Clear Rearranged Particles
        _SHADER.SetBuffer(clearRearrangedParticlesKernel, "rearranged_particles", rearrangeParticlesBuffer);
        _SHADER.SetBuffer(clearRearrangedParticlesKernel, "particle_in_outflow", particlesInOutflowBuffer);
        // 3. Update Counts
        _SHADER.SetBuffer(updateCountsKernel,"particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(updateCountsKernel,"inflowBounds", inflowBoundsBuffer);
        _SHADER.SetBuffer(updateCountsKernel,"outflowBounds", outflowBoundsBuffer);
        _SHADER.SetBuffer(updateCountsKernel,"boundCounts", boundCountsBuffer);
        _SHADER.SetBuffer(updateCountsKernel,"particle_in_outflow", particlesInOutflowBuffer);
        /*
        // 4. Prefix Sum
        _SHADER.SetBuffer(prefixSumKernel, "particle_in_outflow", particlesInOutflowBuffer);
        _SHADER.SetBuffer(prefixSumKernel, "particle_offsets", particleOffsetsBuffer);
        _SHADER.SetBuffer(prefixSumKernel, "particle_sums_buffer", particleSumsBuffer2);
        // 5. Sum Blocks
            // This is set in the update loop, so we don't do any setting here
        // 6. Add Sums
        _SHADER.SetBuffer(addSumsKernel, "particle_offsets", particleOffsetsBuffer);
        // 7. Rearrange Particles
        _SHADER.SetBuffer(rearrangeParticlesKernel, "rearranged_particles", rearrangeParticlesBuffer);
        _SHADER.SetBuffer(rearrangeParticlesKernel, "particle_offsets", particleOffsetsBuffer);
        */
        // 8. Transpose Particles
        _SHADER.SetBuffer(transposeParticlesKernel, "rearranged_particles", rearrangeParticlesBuffer);
        _SHADER.SetBuffer(transposeParticlesKernel, "particles", _BM.PARTICLES_BUFFER);
        _SHADER.SetBuffer(transposeParticlesKernel, "inflow_positions", inflowPositionsBuffer);
        _SHADER.SetBuffer(transposeParticlesKernel, "particle_velocities", _BM.PARTICLES_VELOCITIES_BUFFER);
        _SHADER.SetBuffer(transposeParticlesKernel, "particles_to_transpose", particlesToTransposeBuffer);

        // We must calculate the inflow positions first before anything
        _SHADER.Dispatch(calculateInflowPositionsKernel, Mathf.CeilToInt((float)_numThreshold/16f), 1, 1);

        // Finally, initialize
        _initialized = true;
    }

    void Update() {
        // We can only run this if we are initialized
        if (!TRANSPOSING || !_initialized) return;

        // First, we clear the counts for our inflow and outflow regions
        boundCountsBuffer.SetData(new int[2]{0,0});
        _SHADER.Dispatch(clearRearrangedParticlesKernel, _numParticleBlocks, 1, 1);
        // Second, we count how many entities are present inside each bound
        _SHADER.Dispatch(updateCountsKernel,_numParticleBlocks,1,1);

        // Grab the particle count results
        particlesInOutflowBuffer.GetData(particles_in_outflow_array);
        boundCountsBuffer.GetData(bound_counts_array);

        // Here, we check our timing. 
        // We need to check if we fit the necessary pre-requisites for transposing
        if (bound_counts_array[0] == 0 && bound_counts_array[1] >= _numThreshold) {
            // We need to reposition the particles now
            // First, rearrange the results from `particles_in_outflow_array` based on the index 1 value of each value
            particles_in_outflow_array = particles_in_outflow_array.OrderBy((d) => d[1]).ToArray();
            // Second, we set the results of the sorting inside our shader
            particlesToTransposeBuffer.SetData(particles_in_outflow_array);
            _SHADER.Dispatch(transposeParticlesKernel, Mathf.CeilToInt((float)_numThreshold / 16f), 1, 1);
        }


        /*
         // Prefix Sum
        _SHADER.Dispatch(prefixSumKernel, _numParticleBlocks, 1, 1);
        // Sum Blocks
        bool swap = false;
        for (int d = 1; d < _numParticleBlocks; d *= 2) {
            _SHADER.SetBuffer(sumBlocksKernel, "particle_sums_buffer_in", swap ? particleSumsBuffer1 : particleSumsBuffer2);
            _SHADER.SetBuffer(sumBlocksKernel, "particle_sums_buffer", swap ? particleSumsBuffer2 : particleSumsBuffer1);
            _SHADER.SetInt("d", d);
            _SHADER.Dispatch(sumBlocksKernel, Mathf.CeilToInt((float)_numParticleBlocks / 64f), 1, 1);
            swap = !swap;
        }
        // Add Sums
        _SHADER.SetBuffer(addSumsKernel, "particle_sums_buffer_in", swap ? particleSumsBuffer1 : particleSumsBuffer2);
        _SHADER.Dispatch(addSumsKernel, _numParticleBlocks, 1, 1);
        // Rearrange Boids
        _SHADER.Dispatch(rearrangeParticlesKernel, _numParticleBlocks, 1, 1);

        // We need to check if we fit the necessary pre-requisites for transposing
        boundCountsBuffer.GetData(bound_counts_array);
        if (bound_counts_array[0] == 0 && bound_counts_array[1] >= _numThreshold) {
            // We need to reposition the particles now
            Debug.Log("Must Transpose Particles");
            _SHADER.Dispatch(transposeParticlesKernel, Mathf.CeilToInt((float)_numThreshold / 16f), 1, 1);
        }

        // For debugging, we grab the rearranged particle indices
        rearrangeParticlesBuffer.GetData(rearranged_particles_array);
        */
        
    }

    void OnDestroy() {
        if (!_initialized) return;
        inflowBoundsBuffer.Release();
        outflowBoundsBuffer.Release();
        boundCountsBuffer.Release();
        particlesInOutflowBuffer.Release();
        particleOffsetsBuffer.Release();
        particleSumsBuffer1.Release();
        particleSumsBuffer2.Release();
        rearrangeParticlesBuffer.Release();
    }
}
