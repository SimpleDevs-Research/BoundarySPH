using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TestInheritance : MonoBehaviour
{
    public Base _base;

    // Start is called before the first frame update
    void Start() {
        _base.PrintHello();
    }
}
