using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DebugRotate : MonoBehaviour
{
    public enum Axis {
        X,Y,Z,XY,XZ,YZ
    }
    public Axis axis = Axis.Y;
    public Transform centerOfRotation;
    public float speed = 20f;
    // Update is called once per frame

    void Awake() {
        if (centerOfRotation == null) centerOfRotation = this.transform;
    }
    
    void Update(){
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
            default:
                a = Vector3.forward;
                break;
        }
        transform.RotateAround(centerOfRotation.position, a, speed * Time.deltaTime);
    }
}
