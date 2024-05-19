using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif
using System;

[CustomEditor(typeof(ReplayManager))]
public class ParticleManagerEditor : Editor
{

    ReplayManager _manager;

    public override void OnInspectorGUI() {
        _manager = (ReplayManager)target;
        
        DrawDefaultInspector();
        
        if (GUILayout.Button("Read Positions File Data")) _manager.ReadPositionsFileData();
        if (GUILayout.Button("Clear Positions Data")) _manager.ClearPositionsData();

        if (GUILayout.Button("Read Velocities File Data")) _manager.ReadVelocitiesFileData();
        if (GUILayout.Button("Clear Velocities Data")) _manager.ClearVelocitiesData();
        
    }
}
