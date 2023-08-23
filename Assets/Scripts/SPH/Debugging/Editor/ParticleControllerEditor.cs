using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif
using System;

[CustomEditor(typeof(ParticleController))]
public class ParticleControllerEditor : Editor
{

    ParticleController _pc;

    public override void OnInspectorGUI() {
        _pc = (ParticleController)target;
        
        DrawDefaultInspector();
        
        if (GUILayout.Button("Start Recording")) _pc.StartRecording();
        
    }
}
