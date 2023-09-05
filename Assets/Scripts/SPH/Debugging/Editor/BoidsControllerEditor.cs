using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif
using System;

[CustomEditor(typeof(BoidsController))]
public class BoidsControllerEditor : Editor
{

    BoidsController _bc;

    public override void OnInspectorGUI() {
        _bc = (BoidsController)target;
        
        DrawDefaultInspector();
        
        //if (GUILayout.Button("Generate Boid Data From Transforms")) _bc.GenerateBoidDataFromTransforms();
        //if (GUILayout.Button("Generate Boid Transforms And Data")) _bc.GenerateBoidTransformsAndData();
    }
}
