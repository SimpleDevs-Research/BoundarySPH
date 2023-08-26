using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using OP = ObstaclePrimitives.Structs;

public class CheckProjectionAccuracy : MonoBehaviour
{
    public BufferManager _BM;

    [System.Serializable]
    public class DebugSetup {
        [HideInInspector] public Vector3 particlePosition;
        public Transform targetObstacle;
        public float3 raycastProjection;
        public float3 methodProjection;
        public float displacement;
    }

    [SerializeField] private MeshObsGPU obstacleManager = null;
    public List<DebugSetup> debugSetups = new List<DebugSetup>();

    [SerializeField] private bool showProjections = true;

    void OnDrawGizmos() {
        if (!Application.isPlaying || !showProjections) return;
        for(int i = 0; i < debugSetups.Count; i++) {
            Gizmos.color = Color.black;
            Gizmos.DrawLine(debugSetups[i].particlePosition, debugSetups[i].targetObstacle.position);
            Gizmos.color = new Vector4(0f,0f,1f,0.5f);
            Gizmos.DrawSphere(debugSetups[i].methodProjection,0.2f);
            Gizmos.color = new Vector4(1f,0f,0f,1f);
            Gizmos.DrawSphere(debugSetups[i].raycastProjection,0.1f);
        }
    }

    public ComputeBuffer particleBuffer;
    public ComputeBuffer projectionsBuffer;
    void Update() {
        if (obstacleManager == null) return;
        
        OP.Particle[] particles_array = new OP.Particle[obstacleManager.numParticles];
        OP.Projection[] projections_array = new OP.Projection[obstacleManager.numParticles];

        _BM.PARTICLES_BUFFER.GetData(particles_array);
        _BM.PARTICLES_EXTERNAL_FORCES_BUFFER.GetData(projections_array);
        
        for(int i = 0; i < debugSetups.Count; i++) {
            // Get the projection position. This is the one calculated by our method
            debugSetups[i].methodProjection = new Vector3(projections_array[i].position[0],projections_array[i].position[1],projections_array[i].position[2]);
            // We need to calculate the projection based on SphereCast
            debugSetups[i].particlePosition = new Vector3(particles_array[i].position[0],particles_array[i].position[1],particles_array[i].position[2]);
            Vector3 direction = debugSetups[i].targetObstacle.position - debugSetups[i].particlePosition;
            Vector3 closestPoint = Physics.ClosestPoint(
                debugSetups[i].particlePosition,
                debugSetups[i].targetObstacle.GetComponent<Collider>(),
                debugSetups[i].targetObstacle.position,
                debugSetups[i].targetObstacle.rotation
            );
            // Calculate the displacement
            debugSetups[i].raycastProjection = closestPoint;
            debugSetups[i].displacement = Vector3.Distance(closestPoint, debugSetups[i].methodProjection);
        }
    }
}
