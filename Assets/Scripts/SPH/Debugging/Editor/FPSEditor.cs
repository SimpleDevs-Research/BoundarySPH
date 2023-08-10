using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif

[CustomEditor(typeof(FPS))]
public class FPSEditor : Editor
{
    public override void OnInspectorGUI() {
        FPS fps = (FPS)target;
        DrawDefaultInspector();
        
        /*
        if (GUILayout.Button("Deactivate Debug Manager")) {
            fps.DeactivateDebugManager();
        }
        */
    }
}
