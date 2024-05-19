using System;
using UnityEngine;
using Unity;

public class TakeScreenCapture : MonoBehaviour
{
    public string imageName = null;

    void Update() {
        if (Input.GetKeyDown(KeyCode.Space)) {
            DateTime dt = DateTime.Now;
            string saveName = (IsNullOrWhiteSpace(imageName))
                ? dt.ToString("yyyy-MM-dd\\THH:mm:ss\\Z")
                : $"{imageName}.png";
            ScreenCapture.CaptureScreenshot(saveName, 10);
            Debug.Log("Took Screenshot!");
        }
    }

    public static bool IsNullOrWhiteSpace(string value) {
        if (value != null) {
            for (int i = 0; i < value.Length; i++) {
                if (!char.IsWhiteSpace(value[i])) {
                    return false;
                }
            }
        }
        return true;
    }
}
