using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif

[CustomEditor(typeof(PointCloudObstacleManager))]
public class PointCloudObstacleManagerEditor : Editor
{
    public override void OnInspectorGUI() {
        PointCloudObstacleManager manager = (PointCloudObstacleManager)target;

        DrawDefaultInspector();        
        if (GUILayout.Button("Preprocess Point Clouds")) {
            manager.ManuallyUpdate();
        }
    }
}
