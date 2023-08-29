using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;

public class BufferManager : MonoBehaviour
{
    [Header("=== PARTICLE-related BUFFERS === ")]
    public ComputeBuffer PARTICLES_BUFFER;                      // Stores particles' position data
    public ComputeBuffer PARTICLES_VELOCITIES_BUFFER;           // Stores particles' velocity data
    public ComputeBuffer PARTICLES_DENSITIES_BUFFER;            // Stores particles' density data
    public ComputeBuffer PARTICLES_PRESSURES_BUFFER;            // Stores particles' pressure data (within their kernel)
    public ComputeBuffer PARTICLES_PRESSURE_FORCES_BUFFER;      // Stores the pressure force experienced BY a particle based on its kernel radius
    public ComputeBuffer PARTICLES_VISCOSITY_FORCES_BUFFER;     // Stores the viscosity force experienced BY a particle based on its kernel radius
    public ComputeBuffer PARTICLES_EXTERNAL_FORCES_BUFFER;      // Stores the external forces experienced BY a particle due to forces such as gravity

    [Header("=== GRID-related BUFFERS")]
    public ComputeBuffer PARTICLES_GRID_BUFFER;                 // The grid representation for particle position partitioning and hashing
    public ComputeBuffer PARTICLES_PRESSURE_GRID_BUFFER;        // The grid representation for particle pressure partitioning and hashing

    [Header("=== MESH OBSTACLE-related BUFFERS")]
    public ComputeBuffer MESHOBS_TRANSLATION_FORCES_BUFFER;
    public ComputeBuffer MESHOBS_TORQUE_FORCES_BUFFER;

    [Header("=== DEBUG SETTINGS ===")]
    [SerializeField] private bool _verbose = true;
    [SerializeField] private OP.GridCell[] _particles_grid_array;
    [SerializeField] private float3[] _particles_velocities_array;
    [SerializeField] private float[] _particles_densities_array;
    [SerializeField] private float[] _particles_pressures_array;

    public void InitializeParticleBuffers(int numParticles = 1, int numGridCells = 1) {
        PARTICLES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*6+sizeof(int));
        PARTICLES_GRID_BUFFER = new ComputeBuffer(numGridCells, sizeof(int));

        PARTICLES_VELOCITIES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_DENSITIES_BUFFER = new ComputeBuffer(numParticles, sizeof(float));
        PARTICLES_PRESSURES_BUFFER = new ComputeBuffer(numParticles, sizeof(float));

        PARTICLES_PRESSURE_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_VISCOSITY_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(float)*3);
        PARTICLES_EXTERNAL_FORCES_BUFFER = new ComputeBuffer(numParticles, sizeof(uint) + sizeof(int)*8 + sizeof(float)*27);

        if (_verbose) Debug.Log("[BUFFER MANAGER] Particle buffers initialized!");
    }

    void Update() {
        if (_particles_grid_array.Length > 0) PARTICLES_GRID_BUFFER.GetData(_particles_grid_array);
        if (_particles_velocities_array.Length > 0) PARTICLES_VELOCITIES_BUFFER.GetData(_particles_velocities_array);
        if (_particles_densities_array.Length > 0) PARTICLES_DENSITIES_BUFFER.GetData(_particles_densities_array);
        if (_particles_pressures_array.Length > 0) PARTICLES_PRESSURES_BUFFER.GetData(_particles_pressures_array);
    }

    void OnDestroy() {
        if (PARTICLES_BUFFER != null) PARTICLES_BUFFER.Release();
        if (PARTICLES_GRID_BUFFER != null) PARTICLES_GRID_BUFFER.Release();
        if (PARTICLES_VELOCITIES_BUFFER != null) PARTICLES_VELOCITIES_BUFFER.Release();
        if (PARTICLES_DENSITIES_BUFFER != null) PARTICLES_DENSITIES_BUFFER.Release();
        if (PARTICLES_PRESSURES_BUFFER != null) PARTICLES_PRESSURES_BUFFER.Release();
        if (PARTICLES_PRESSURE_GRID_BUFFER != null) PARTICLES_PRESSURE_GRID_BUFFER.Release();
        if (PARTICLES_PRESSURE_FORCES_BUFFER != null) PARTICLES_PRESSURE_FORCES_BUFFER.Release();
        if (PARTICLES_VISCOSITY_FORCES_BUFFER != null) PARTICLES_VISCOSITY_FORCES_BUFFER.Release();
        if (PARTICLES_EXTERNAL_FORCES_BUFFER != null) PARTICLES_EXTERNAL_FORCES_BUFFER.Release();

        if (MESHOBS_TRANSLATION_FORCES_BUFFER != null) MESHOBS_TRANSLATION_FORCES_BUFFER.Release();
        if (MESHOBS_TORQUE_FORCES_BUFFER != null) MESHOBS_TORQUE_FORCES_BUFFER.Release();
    }
}

