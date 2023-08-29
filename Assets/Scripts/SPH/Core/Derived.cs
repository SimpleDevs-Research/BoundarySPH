using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Derived : Base
{
    public override void PrintHello() {
        Debug.Log(base.message);
    }
}
