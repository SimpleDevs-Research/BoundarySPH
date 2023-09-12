using UnityEngine;
#if UNITY_EDITOR
    using UnityEditor;
#endif
using System;

[CustomEditor(typeof(MeshObsGPU))]
public class MeshObsGPUEditor : Editor
{

    MeshObsGPU controller;

    public override void OnInspectorGUI() {
        controller = (MeshObsGPU)target;
        
        DrawDefaultInspector();
        
        if (GUILayout.Button("Generate Boid Transforms")) controller.GenerateBoidTransforms();
        if (GUILayout.Button("Reset Boids to Defaults")) controller.ResetBoidsToDefaults();
        if (GUILayout.Button("Delete Boids")) controller.DeleteBoids();
        //if (GUILayout.Button("Generate Vertices")) _meshObs.GenerateMeshVertices();

        //if (GUILayout.Button("Preprocess Obstacles"))   manager.PreprocessObstacles();
        //if (GUILayout.Button("Update Obstacles"))       manager.UpdateObstacles();

        /*
        serializedObject.Update();
        GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_GRID"));
        GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_SHADER"));
        GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_OBSTACLES"));

        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_VERTICES"), EditorListOption.ListLabel);
        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_TRIANGLES"), EditorListOption.ListLabel);
        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_BOUNDS"), EditorListOption.ListLabel);
        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_VECTOR_FIELD"), EditorListOption.ListLabel);
        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_CELLS_IN_BOUNDS"), EditorListOption.ListLabel);
        //GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_PROJECTION_TEST_RESULTS"), EditorListOption.ListLabel);

        // GPU_ObstacleManagerEditor.RenderVariable(manager, serializedObject.FindProperty("_combined"));
        //serializedObject.ApplyModifiedProperties();
        */
        /*
        if (manager.showEditorControls == GPU_ObstacleManager.ShowHideEnum.Hide) return;
        
        SerializedProperty obstacles = serializedObject.FindProperty("_OBSTACLES");
        for (int i = 0; i < obstacles.arraySize; i++) {
            if (GUILayout.Button($"Mesh Details: Obstacle {i+1}")) manager.PrintObstacleDetails(i);
        }

        if (GUILayout.Button("Manually Initialize"))        manager.Initialize();
        if (GUILayout.Button("Update Obstacles"))           manager.UpdateObstacles();
        if (GUILayout.Button("Manually Check For Updates")) manager.ManualUpdate();
        if (GUILayout.Button("FULL RESET"))                 manager.FullReset();
        */
    }

    public static void RenderVariable(GPU_ObstacleManager manager, SerializedProperty field, EditorListOption options = EditorListOption.All) {
        // Exit early if it's a non-array element
        if(!field.isArray) {
            EditorGUILayout.PropertyField(field);
            return;
        }
        
        // Derive boolean flags
        bool showListLabel = (options & EditorListOption.ListLabel) != 0;
        bool showListSize = (options & EditorListOption.ListSize) != 0;
        bool showListElements = (options & EditorListOption.ListElements) != 0;
        bool showButtons = (options & EditorListOption.Buttons) != 0;

        if (showListLabel) EditorGUILayout.PropertyField(field);
        else if (showListSize) EditorGUILayout.PropertyField(field.FindPropertyRelative("Array.size"));
		
        EditorGUI.indentLevel += 1;

        if (showListElements) {
            for (int i = 0; i < field.arraySize; i++) {
                GPU_ObstacleManagerEditor.RenderVariable(manager, field.GetArrayElementAtIndex(i));
            }	
        }

        EditorGUI.indentLevel -= 1;
    }

}
