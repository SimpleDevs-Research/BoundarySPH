using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DebugRotate : MonoBehaviour
{
    public enum Axis {
        X,Y,Z,XY,XZ,YZ,XYZ
    }
    public Axis axis = Axis.Y;
    public Transform centerOfRotation;
    public float deltaTime = -1f;
    private float _deltaTime;
    public float speed = 20f;
    // Update is called once per frame

    void Awake() {
        if (centerOfRotation == null) centerOfRotation = this.transform;
    }
    
    void Update(){
        if (deltaTime < 0f) _deltaTime = Time.deltaTime;
        else _deltaTime = deltaTime;
        Vector3 a;
        switch(axis) {
            case Axis.X:
                a = Vector3.right;
                break;
            case Axis.Y:
                a = Vector3.up;
                break;
            case Axis.XY:
                a = (Vector3.right + Vector3.up).normalized;
                break;
            case Axis.XZ:
                a = (Vector3.right + Vector3.forward).normalized;
                break;
            case Axis.YZ:
                a = (Vector3.up + Vector3.forward).normalized;
                break;
            case Axis.XYZ:
                a = (Vector3.up + Vector3.forward + Vector3.right).normalized;
                break;
            default:
                a = Vector3.forward;
                break;
        }
        transform.RotateAround(centerOfRotation.position, a, speed * _deltaTime);
    }
}
