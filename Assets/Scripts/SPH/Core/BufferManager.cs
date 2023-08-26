using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BufferManager : MonoBehaviour
{
    [Header("=== PARTICLE-related BUFFERS === ")]
    public ComputeBuffer PARTICLES_BUFFER;
    public ComputeBuffer PARTICLES_GRID_BUFFER;
    public ComputeBuffer PARTICLES_VELOCITIES_BUFFER;
    public ComputeBuffer PARTICLES_DENSITIES_BUFFER;
    public ComputeBuffer PARTICLES_PRESSURES_BUFFER;
    public ComputeBuffer PARTICLES_PRESSURE_FORCES_BUFFER;
    public ComputeBuffer PARTICLES_VISCOSITY_FORCES_BUFFER;
    public ComputeBuffer PARTICLES_EXTERNAL_FORCES_BUFFER;

    [Header("=== MESH OBSTACLE-related BUFFERS")]
    public ComputeBuffer MESHOBS_TRANSLATION_FORCES_BUFFER;
    public ComputeBuffer MESHOBS_TORQUE_FORCES_BUFFER;

    void OnDestroy() {
        if (PARTICLES_BUFFER != null) PARTICLES_BUFFER.Release();
        if (PARTICLES_GRID_BUFFER != null) PARTICLES_GRID_BUFFER.Release();
        if (PARTICLES_VELOCITIES_BUFFER != null) PARTICLES_VELOCITIES_BUFFER.Release();
        if (PARTICLES_DENSITIES_BUFFER != null) PARTICLES_DENSITIES_BUFFER.Release();
        if (PARTICLES_PRESSURES_BUFFER != null) PARTICLES_PRESSURES_BUFFER.Release();
        if (PARTICLES_PRESSURE_FORCES_BUFFER != null) PARTICLES_PRESSURE_FORCES_BUFFER.Release();
        if (PARTICLES_VISCOSITY_FORCES_BUFFER != null) PARTICLES_VISCOSITY_FORCES_BUFFER.Release();
        if (PARTICLES_EXTERNAL_FORCES_BUFFER != null) PARTICLES_EXTERNAL_FORCES_BUFFER.Release();

        if (MESHOBS_TRANSLATION_FORCES_BUFFER != null) MESHOBS_TRANSLATION_FORCES_BUFFER.Release();
        if (MESHOBS_TORQUE_FORCES_BUFFER != null) MESHOBS_TORQUE_FORCES_BUFFER.Release();
    }
}

