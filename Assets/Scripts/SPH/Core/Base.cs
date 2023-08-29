using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Base : MonoBehaviour
{
    public string message = "Hello";
    public virtual void PrintHello() {
        Debug.Log(message);
    }
}
