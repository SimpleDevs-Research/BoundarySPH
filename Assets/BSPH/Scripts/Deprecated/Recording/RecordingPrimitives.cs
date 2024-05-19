using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using SerializableTypes;

namespace RecordingPrimitives
{
    [System.Serializable]
    public class RecordPack {
        public float timestamp;
        public int particle_id;
        public int encoded_position;
        public int encoded_velocity;
    }
}
