using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;
using OP = ObstaclePrimitives.Structs;

public class PointCloudObstacleManager : MonoBehaviour
{
    [System.Serializable]
    public class PointObstacle {
        public Transform obstacle = null;
        public Color gizmosColor = Color.red;
        public bool isActive = true;
        [HideInInspector, ReadOnlyInsp, Tooltip("World-scale positions of boundary particles")] public List<OP.Particle> boundaryParticles = new List<OP.Particle>();
    }

    public float kernelRadius = 1.5f;
    public List<PointObstacle> obstacles = new List<PointObstacle>();
    [SerializeField] public List<OP.Particle> boundaryParticles;
    [SerializeField] public int numBoundaryParticles = 0;

    [SerializeField] private Color boundaryParticleColor = Color.red;
    //private bool _SHOW_GIZMOS => boundaryParticleColor.a > 0f;

    void OnDrawGizmos() {
        //if (!_SHOW_GIZMOS) return;
        //Gizmos.color = boundaryParticleColor;
        foreach(PointObstacle obstacle in obstacles) {
            if (obstacle.gizmosColor.a == 0f) continue;
            Gizmos.color = obstacle.gizmosColor;
            foreach(OP.Particle particle in obstacle.boundaryParticles) {
                Gizmos.DrawSphere(particle.position,kernelRadius);
            }
        }
        /*
        foreach(ParticleController.Particle particle in boundaryParticles) {
            Gizmos.DrawSphere(particle.position,kernelRadius);
        }
        */
    }

    // We call start becuase ParticleController does its initialization on Awake(). We need to update the Particles compute buffer with our own boundary particles!
    private void Awake() {
        Debug.Log($"Number of Boundary Particles: {numBoundaryParticles}");
    }

    public void ManuallyUpdate() {
        boundaryParticles = new List<OP.Particle>();
        foreach(PointObstacle obs in obstacles) {
            CalculateBoundaryPoints(obs);
            boundaryParticles.AddRange(obs.boundaryParticles);
        }
        numBoundaryParticles = boundaryParticles.Count;
    }

    private void CalculateBoundaryPoints(PointObstacle obs) {
        // Initializations
        obs.boundaryParticles = new List<OP.Particle>();
         if (obs.obstacle == null || kernelRadius <= 0f) return;

        Mesh mesh = obs.obstacle.GetComponent<MeshFilter>().sharedMesh;
        int[] triangles = mesh.triangles;
        Vector3[] vertices = mesh.vertices;

        // world positions for each triangle vertex, as well as centroid
        Vector3 wv1, wv2, wv3, centroid;
        // vector b/w wv1 and wv2, and wv3-closestPointOn_wv1-wv2
        Vector3 wv1wv2, wv3_closest, wv3_wv1wv2;
        // How many points can we fit on each axis?
        int wv1wv2_num, wv3_wv1wv2_num;
        // When we're adding points, we need reference points!
        Vector3 tempPos,posToAdd,normDir;
        // Plane data
        Plane plane;
        // Particle data
        OP.Particle particle;

        for(int t = 0; t < triangles.Length; t += 3) {
            // Calculate world positions for each triangle vertex
            wv1 = obs.obstacle.TransformPoint(mesh.vertices[triangles[t]]);
            wv2 = obs.obstacle.TransformPoint(mesh.vertices[triangles[t+1]]);
            wv3 = obs.obstacle.TransformPoint(mesh.vertices[triangles[t+2]]);
            // Calculate norm of the triangle
            normDir = Vector3.Cross(wv2 - wv1, wv3 - wv1).normalized;
            // Calculate centroid of three points
            centroid = (wv1 + wv2 + wv3)/3f;
            // Calculate world-space vector b/w wv1 and wv2
            wv1wv2 = wv2 - wv1;
            // Calculate closest point to wv3 from `wv1_wv2`
            wv3_closest = FindNearestPointOnLine(wv1,wv2,wv3);
            // Calculate vector from wv3 to wv3_closest
            wv3_wv1wv2 = wv3 - wv3_closest;

            // Calculate how many points we are allowed on each axis
            wv1wv2_num = Mathf.CeilToInt(wv1wv2.magnitude / kernelRadius);
            wv3_wv1wv2_num = Mathf.CeilToInt(wv3_wv1wv2.magnitude / kernelRadius);

            // Get the plane that encompasses wv1,wv2,wv3
            plane = new Plane(wv1, wv2, wv3);

            // Add centroid
            posToAdd = transform.InverseTransformPoint(centroid);
            particle = new OP.Particle();
            particle.position = new(posToAdd.x, posToAdd.y, posToAdd.z);
            obs.boundaryParticles.Add(particle);

            // Add points by iterating across rectangle
            for(int x = 0; x < wv1wv2_num; x++) {
                float wv3_wv1wv2_oddOffset = (x % 2 == 1) ? kernelRadius : 0f;
                for(int y = 0; y < wv3_wv1wv2_num; y++) {
                    /*
                    tempPos = wv1 
                        + (wv1wv2.normalized * (x*(pointCloudResolution*1.5f) + (pointCloudResolution*1.5f))) 
                        + (wv3_wv1wv2.normalized * (y*(pointCloudResolution*1.5f) + (pointCloudResolution*1.5f)));
                    tempPos = plane.ClosestPointOnPlane(tempPos);
                    */
                    tempPos = wv1 
                        + (wv1wv2.normalized * (x*kernelRadius)) 
                        + (wv3_wv1wv2.normalized * (y*kernelRadius));
                    tempPos = plane.ClosestPointOnPlane(tempPos);
                    if(PointInTriangle(tempPos, wv1,wv2,wv3)) {
                        // Convert to local space
                        posToAdd = transform.InverseTransformPoint(tempPos);
                        particle = new OP.Particle();
                        particle.position = new(posToAdd.x, posToAdd.y, posToAdd.z);
                        obs.boundaryParticles.Add(particle);
                        // We also add one layer extra, in case
                        /*
                        particle = new ParticleController.Particle();
                        particle.position = new(posToAdd.x - normDir.x * kernelRadius, posToAdd.y  - normDir.y * kernelRadius, posToAdd.z - normDir.z * kernelRadius);
                        obs.boundaryParticles.Add(particle);
                        */
                    }
                }
            }
        }
    }

    public static bool SameSide(Vector3 p1, Vector3 p2, Vector3 a, Vector3 b) {
        Vector3 cp1 = Vector3.Cross(b-a, p1-a);
        Vector3 cp2 = Vector3.Cross(b-a, p2-a);
        return Vector3.Dot(cp1, cp2) >= -1f;
    }

    public static bool PointInTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c) {
        return SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b);
    }
    public static bool PointInTriangle(Vector3 p, Vector3[] triangle) {
        if (triangle.Length != 3) return false;
        return SameSide(p,triangle[0], triangle[1],triangle[2]) && SameSide(p,triangle[1], triangle[0],triangle[2]) && SameSide(p,triangle[2], triangle[0],triangle[1]);
    }

    public static Vector3 FindNearestPointOnLine(Vector3 origin, Vector3 end, Vector3 point) {
        //Get heading
        Vector3 heading = (end - origin);
        float magnitudeMax = heading.magnitude;
        heading.Normalize();

        //Do projection from the point but clamp it
        Vector3 lhs = point - origin;
        float dotP = Vector3.Dot(lhs, heading);
        dotP = Mathf.Clamp(dotP, 0f, magnitudeMax);
        return origin + heading * dotP;
    }

    /*
    public List<PointCloudObstacle> obstacles;

    [SerializeField]
    private float3[] _points = new float3[0];
    public float3[] points { get => _points; set {} }
    public int numPoints = 0;
    public ComputeBuffer pointsBuffer;

    public bool waitForInitialization = false;

    // Start is called before the first frame update
    void Awake() {
        if (!waitForInitialization) Initialize();
    }

    public void Initialize() {
        GetPoints();
        pointsBuffer = new ComputeBuffer(_points.Length, sizeof(float)*3);
        pointsBuffer.SetData(_points);
    }

    // Update is called once per frame
    void Update() {
        GetPoints();
        pointsBuffer = new ComputeBuffer(_points.Length, sizeof(float)*3);
        pointsBuffer.SetData(_points);
    }

    public void GetPoints() {
        numPoints = 0;
        foreach(PointCloudObstacle obstacle in obstacles) {
            numPoints += obstacle.points.Count;
        }

        if (numPoints == 0) {
            _points = new float3[1] { new(0f,0f,0f) };
        }
        else {
            _points = new float3[numPoints];
            foreach(PointCloudObstacle obstacle in obstacles) {
                _points.Concat(obstacle.pointsF3);
            }
        }
    }
    */
}
