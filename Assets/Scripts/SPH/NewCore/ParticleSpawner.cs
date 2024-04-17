using UnityEngine;
using Unity.Mathematics;

public class ParticleSpawner : MonoBehaviour
{
    public int3 numParticlesPerAxis;
    public float spawnDistanceBetweenParticles;
    private float3 size;
    public float3 initialVel;
    public float jitterStrength;
    public bool showSpawnBounds;
    public Color spawnBoundsColor = Color.yellow;

    [Header("Info")]
    public int debug_numParticles;

    public SpawnData GetSpawnData() {
        int numPoints = numParticlesPerAxis.x * numParticlesPerAxis.y * numParticlesPerAxis.z;
        ParticleStruct[] particles = new ParticleStruct[numPoints];
        float3[] positions = new float3[numPoints];
        float3[] velocities = new float3[numPoints];

        Vector3 center = transform.position;
        int i = 0;

        for (int x = 0; x < numParticlesPerAxis.x; x++) {
            for (int y = 0; y < numParticlesPerAxis.y; y++) {
                for (int z = 0; z < numParticlesPerAxis.z; z++) {
                    float tx = x / (numParticlesPerAxis.x - 1f);
                    float ty = y / (numParticlesPerAxis.y - 1f);
                    float tz = z / (numParticlesPerAxis.z - 1f);

                    float px = (tx - 0.5f) * size.x + center.x;
                    float py = (ty - 0.5f) * size.y + center.y;
                    float pz = (tz - 0.5f) * size.z + center.z;
                    float3 jitter = UnityEngine.Random.insideUnitSphere * jitterStrength;
                    positions[i] = new float3(px, py, pz) + jitter;
                    particles[i] = new ParticleStruct() { position = positions[i], force = new float3(0,0,0), render = 0 };
                    velocities[i] = initialVel;
                    i++;
                }
            }
        }

        return new SpawnData() { particles = particles, positions = positions, velocities = velocities };
    }

    public struct ParticleStruct {
        public float3 position;
        public float3 force;
        public int render;
    }

    public struct SpawnData {
        public ParticleStruct[] particles;
        public float3[] positions;
        public float3[] velocities;
    }

    void OnValidate() {
        // User specifies # of particles OR distance between particles. Let's double-check.
        size.x = (numParticlesPerAxis.x-1) * spawnDistanceBetweenParticles;
        size.y = (numParticlesPerAxis.y-1) * spawnDistanceBetweenParticles;
        size.z = (numParticlesPerAxis.z-1) * spawnDistanceBetweenParticles;
        debug_numParticles = numParticlesPerAxis.x * numParticlesPerAxis.y * numParticlesPerAxis.z;
    }

    void OnDrawGizmos()
    {
        if (showSpawnBounds && !Application.isPlaying) {
            Gizmos.color = spawnBoundsColor;
            Gizmos.DrawWireCube(transform.position, Vector3.one * size);
        }
    }
}
