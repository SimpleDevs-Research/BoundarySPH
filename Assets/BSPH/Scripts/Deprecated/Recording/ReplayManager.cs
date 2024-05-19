using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using HelperMethods = Helpers.HelperMethods;
using DrawingMethods = Helpers.DrawingMethods;

[ExecuteInEditMode]
public class ReplayManager : MonoBehaviour
{
    [Header("=== SETTINGS / CONFIGS ===")]
    [SerializeField] private TextAsset _positionsFile, _velocitiesFile;
    [SerializeField] private List<Vector4>[] positions, velocities;

    [Header("=== DRAWING CONFIGS ===")]
    [SerializeField] private bool _drawPositions, _drawVelocities;
    [SerializeField] private Color _positionsColor = Color.blue, _velocitiesColor = Color.red;
    [SerializeField]
    private Vector2 _testBounds = new Vector2(0f,10f);
    public Vector2 testBounds => _testBounds;
    private Vector2 _testVals = new Vector2(1f,5f);
    public Vector2 testVals => _testVals;


    [Header("=== HYPERCUBE PARAMETERS ===")]
    [SerializeField] private Vector4 _hypercubeParams = new Vector4(1f,1f,1f,0.25f);
    
    /*
    [System.Serializable]
    public class HyperCube {
        public List<Vector4> positions;
        public float[] bounds;
        public Vector4 dimensions;
        public HyperCube() {
            positions = new List<Vector3>();
            bounds = new float[8]{0f,0f,0f,0f,0f,0f,0f,0f};
            dimensions = new Vector4(0f,0f,0f,0f);
        }
        public bool CheckPosition(Vector3 pos, float t) {
            if (positions.Count == 0) return AddPosition(pos, t);
            float diffX = Mathf.Max(bounds[4],pos.x) - Mathf.Min(bounds[0],pos.x);
            float diffY = Mathf.Max(bounds[5],pos.y) - Mathf.Min(bounds[1],pos.y);
            float diffZ = Mathf.Max(bounds[6],pos.z) - Mathf.Min(bounds[2],pos.z);
            float diffT = Mathf.Max(bounds[7],t) - Mathf.Min(bounds[3],t);
            if (
                diffX <= _hypercubeParams.x && 
                diffY <= _hypercubeParams.y &&
                diffZ <= _hypercubeParams.z &&
                diffT <= _hypercubeParams.w
            ) return AddPosition(pos, t);
            return false;
        }
        public bool AddPosition() {
            positions.Add(new Vector4(pos.x, pos.y, pos.z, t));
            if (pos.x < bounds[0]) bounds[0] = pos.x;
            if (pos.y < bounds[1]) bounds[1] = pos.y;
            if (pos.z < bounds[2]) bounds[2] = pos.z;
            if (t < bounds[3]) bounds[3] = t;
            if (pos.x > bounds[4]) bounds[4] = pos.x;
            if (pos.y > bounds[5]) bounds[5] = pos.y;
            if (pos.z > bounds[6]) bounds[6] = pos.z;
            if (t > bounds[7]) bounds[7] = t;
            dimensions = new Vector4(
                bounds[4] - bounds[0],
                bounds[5] - bounds[1],
                bounds[6] - bounds[2],
                bounds[7] - bounds[3]
            );
            return true;
        }
    }
    */

    void OnDrawGizmos() {
        if (Application.isPlaying) return;
        if (_drawPositions && positions != null && positions.Length > 0) DrawPositions();
        if (_drawVelocities && positions != null && positions.Length > 0 && velocities != null && velocities.Length > 0) DrawVelocities();
    }

    private void DrawPositions() {
        Gizmos.color = _positionsColor;
        foreach(List<Vector4> current_positions in positions) {
            foreach(Vector4 p in current_positions) {
                Gizmos.DrawSphere(new Vector3(p.x,p.y,p.z), 0.005f);
            }
        }
    }

    private void DrawVelocities() {
        Gizmos.color = _velocitiesColor;
        // Across each particle....
        for(int i = 0; i < velocities.Length; i++) {
            // Across each position/velocity:
            for(int j = 0; j < velocities[i].Count; j++) {
                if (i >= positions.Length || j >= positions[i].Count) continue;
                Vector4 p = positions[i][j];
                Vector4 d = velocities[i][j];
                Vector3 pos = new Vector3(p.x,p.y,p.z);
                Vector3 dir = Vector3.Normalize(new Vector3(d.x,d.y,d.z)) * 0.1f;
                DrawingMethods.ForGizmos(p,d,0.1f);
            }
        }
    }

    public void ReadPositionsFileData() {
        // First, attempt to read the CSV files from both the positions and velocities data
        if (_positionsFile == null) return;
        string[] positions_raw = _positionsFile.text.Split(new string[] {"\n"}, StringSplitOptions.None);
        
        // Second, we need to decode this. `positions_raw` is spilt by line, but each individual item in it still needs to be delimited by ","
        // To put it another way, each row represents a timestamp, and each column represents a particle's position
        // In concept, it would be better to keep a list of particles, with each item in that list being its own list that indicates the position of the i-th particle at a certain timestamp
        // In other words, we need to flip the columns and rows around.
        // Firstly, we can identify the number of particles stored in this data by looking at the header row.
        // The number of particles is (# of items in the header) - 1, as one of the columns represents the timestamps
        int numParticles = positions_raw[0].Split(",").Length - 1;

        // We can now generate the positions array
        positions = new List<Vector4>[numParticles];
        for(int i = 0; i < numParticles; i++) {
            positions[i] = new List<Vector4>();
        }

        // As we iterate across each row (which represents a single timestep)
        for(int i = 1; i < positions_raw.Length; i++) {
            if (positions_raw[i].Length == 0) continue;
            // Get the individual values
            string[] row = positions_raw[i].Split(",");
            // The timestamp is the first item in `row`
            float timestamp = float.Parse(row[0]);
            // The remaining items in `row` represent the particle positions at this timestamp.
            // When iterating across them, we have to create a Vector4 that contains the xyz position + the timestamp, then add them to the
            //   proper ID inisde `positions`
            for(int j = 1; j < row.Length; j++) {
                string[] p_string = row[j].Split("|");
                Vector4 currentPos = new Vector4(
                    float.Parse(p_string[0]),
                    float.Parse(p_string[1]),
                    float.Parse(p_string[2]),
                    timestamp
                );
                positions[j-1].Add(currentPos);
            }
        }
    }

    public void ReadVelocitiesFileData() {
        if (_velocitiesFile == null) return;
        string[] velocities_raw = _velocitiesFile.text.Split(new string[]{"\n"},StringSplitOptions.None);

        int numParticles = velocities_raw[0].Split(",").Length - 1;

        velocities = new List<Vector4>[numParticles];
        for(int i = 0; i < numParticles; i++) {
            velocities[i] = new List<Vector4>();
        }

        for(int i = 1; i < velocities_raw.Length; i++) {
            if (velocities_raw[i].Length == 0) continue;
            string[] row = velocities_raw[i].Split(",");
            float timestamp = float.Parse(row[0]);
            for(int j = 1; j < row.Length; j++) {
                //Vector3 currentPos = HelperMethods.decodeVector3FromInt(int.Parse(row[j]));
                string[] v_string = row[j].Split("|");
                Vector4 currentVel = new Vector4(
                    float.Parse(v_string[0]),
                    float.Parse(v_string[1]),
                    float.Parse(v_string[2]),
                    timestamp
                );
                velocities[j-1].Add(currentVel);
            }
        }
    }

    public void ClearPositionsData() {
        positions = new List<Vector4>[0];
    }

    public void ClearVelocitiesData() {
        velocities = new List<Vector4>[0];
    }

    public void CalculateTrajectories() {
        if (positions == null || positions.Length == 0) {
            Debug.LogError("Cannot create trajectories due to missing Positions data.");
            return;
        }

        // We will iterate across all columns
        for(int i = 0; i < positions.Length; i++) {

        }
    }
}
