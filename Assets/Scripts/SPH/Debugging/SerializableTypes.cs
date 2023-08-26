using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace SerializableTypes {

    /// <summary> Serializable version of UnityEngine.Vector2. </summary>
    [System.Serializable]
    public struct SVector2 {
        public float x;
        public float y;
 
        public SVector2(float x, float y) {
            this.x = x;
            this.y = y;
        }
 
        public override string ToString() 
            => $"[{x}, {y}]";
 
        public static implicit operator Vector2(SVector2 s) 
            => new Vector2(s.x, s.y);
 
        public static implicit operator SVector2(Vector2 v)
            => new SVector2(v.x, v.y);
 
 
        public static SVector2 operator +(SVector2 a, SVector2 b) 
            => new SVector2(a.x + b.x, a.y + b.y);
 
        public static SVector2 operator -(SVector2 a, SVector2 b)
            => new SVector2(a.x - b.x, a.y - b.y);
 
        public static SVector2 operator -(SVector2 a)
            => new SVector2(-a.x, -a.y);
 
        public static SVector2 operator *(SVector2 a, float m)
            => new SVector2(a.x * m, a.y * m);
 
        public static SVector2 operator *(float m, SVector2 a)
            => new SVector2(a.x * m, a.y * m);
 
        public static SVector2 operator /(SVector2 a, float d)
            => new SVector2(a.x / d, a.y / d);
    }

    /// <summary> Serializable version of UnityEngine.Vector3. </summary>
    [System.Serializable]
    public struct SVector3 {
        public float x;
        public float y;
        public float z;
 
        public SVector3(float x, float y, float z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
 
        public override string ToString() 
            => $"[{x}, {y}, {z}]";

        public override int GetHashCode()
            => x.GetHashCode() ^ (y.GetHashCode() << 2) ^ (z.GetHashCode() >> 2);

        public override bool Equals(object other) {
            if (!(other is SVector3)) return false;
            return Equals((SVector3)other);
        }

        public bool Equals(SVector3 other)
            => x == other.x && y == other.y && z == other.z;

        public static implicit operator Vector3(SVector3 s) 
            => new Vector3(s.x, s.y, s.z);
 
        public static implicit operator SVector3(Vector3 v)
            => new SVector3(v.x, v.y, v.z);
 
 
        public static SVector3 operator +(SVector3 a, SVector3 b) 
            => new SVector3(a.x + b.x, a.y + b.y, a.z + b.z);
 
        public static SVector3 operator -(SVector3 a, SVector3 b)
            => new SVector3(a.x - b.x, a.y - b.y, a.z - b.z);
 
        public static SVector3 operator -(SVector3 a)
            => new SVector3(-a.x, -a.y, -a.z);
 
        public static SVector3 operator *(SVector3 a, float m)
            => new SVector3(a.x * m, a.y * m, a.z * m);
 
        public static SVector3 operator *(float m, SVector3 a)
            => new SVector3(a.x * m, a.y * m, a.z * m);
 
        public static SVector3 operator /(SVector3 a, float d)
            => new SVector3(a.x / d, a.y / d, a.z / d);
        
        public static bool operator ==(SVector3 a, SVector3 b)
            => (Vector3)a == (Vector3)b;
        public static bool operator ==(SVector3 a, Vector3 b)
            => (Vector3)a == b;
        public static bool operator ==(Vector3 a, SVector3 b)
            => a == (Vector3)b;
        public static bool operator !=(SVector3 a, SVector3 b)
            => (Vector3)a != (Vector3)b;
        public static bool operator !=(SVector3 a, Vector3 b)
            => (Vector3)a != b;
        public static bool operator !=(Vector3 a, SVector3 b)
            => a != (Vector3)b;
        
        public float magnitude {
            get { return (float)Mathf.Sqrt(this.x * this.x + this.y * this.y + this.z * this.z); }
        }
    }

    /// <summary> Serializable version of UnityEngine.Vector4. </summary>
    [System.Serializable]
    public struct SVector4 {
        public float x;
        public float y;
        public float z;
        public float w;
 
        public SVector4(float x, float y, float z, float w) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
 
        public override string ToString() 
            => $"[{x}, {y}, {z}, {w}]";
        
        public override int GetHashCode()
            => x.GetHashCode() ^ (y.GetHashCode() << 2) ^ (z.GetHashCode() >> 2) ^ (w.GetHashCode() >> 1);

        public override bool Equals(object other) {
            if (!(other is SVector4)) return false;
            return Equals((SVector4)other);
        }
 
        public static implicit operator Vector4(SVector4 s) 
            => new Vector4(s.x, s.y, s.z, s.w);
 
        public static implicit operator SVector4(Vector4 v)
            => new SVector4(v.x, v.y, v.z, v.w);
 
 
        public static SVector4 operator +(SVector4 a, SVector4 b) 
            => new SVector4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
 
        public static SVector4 operator -(SVector4 a, SVector4 b)
            => new SVector4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
 
        public static SVector4 operator -(SVector4 a)
            => new SVector4(-a.x, -a.y, -a.z, -a.w);
 
        public static SVector4 operator *(SVector4 a, float m)
            => new SVector4(a.x * m, a.y * m, a.z * m, a.w * m);
 
        public static SVector4 operator *(float m, SVector4 a)
            => new SVector4(a.x * m, a.y * m, a.z * m, a.w * m);
 
        public static SVector4 operator /(SVector4 a, float d)
            => new SVector4(a.x / d, a.y / d, a.z / d, a.w / d);

        public static bool operator ==(SVector4 a, SVector4 b)
            => (Vector4)a == (Vector4)b;
        public static bool operator ==(SVector4 a, Vector4 b)
            => (Vector4)a == b;
        public static bool operator ==(Vector4 a, SVector4 b)
            => a == (Vector4)b;
        public static bool operator !=(SVector4 a, SVector4 b)
            => (Vector4)a != (Vector4)b;
        public static bool operator !=(SVector4 a, Vector4 b)
            => (Vector4)a != b;
        public static bool operator !=(Vector4 a, SVector4 b)
            => a != (Vector4)b;

        public float magnitude {
            get { return (float)Mathf.Sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w); }
        }
    }
 
    /// <summary> Serializable version of UnityEngine.Quaternion. </summary>
    [System.Serializable]
    public struct SQuaternion {
        public float x;
        public float y;
        public float z;
        public float w;
 
        public SQuaternion(float x, float y, float z, float w) {
            this.x = x;
            this.y = y;
            this.z = z;
            this.w = w;
        }
 
        public override string ToString()
            => $"[{x}, {y}, {z}, {w}]";
        
        public override int GetHashCode()
            => x.GetHashCode() ^ (y.GetHashCode() << 2) ^ (z.GetHashCode() >> 2) ^ (w.GetHashCode() >> 1);
 
        public override bool Equals(object other) {
            if (!(other is SQuaternion)) return false;
            return Equals((SQuaternion)other);
        }

        public bool Equals(SQuaternion other)
            => x.Equals(other.x) && y.Equals(other.y) && z.Equals(other.z) && w.Equals(other.w);

        public static implicit operator Quaternion(SQuaternion s)
            => new Quaternion(s.x, s.y, s.z, s.w);
 
        public static implicit operator SQuaternion(Quaternion q)
            => new SQuaternion(q.x, q.y, q.z, q.w);
        
        public static bool operator ==(SQuaternion a, SQuaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) < 0.01f;
        public static bool operator ==(SQuaternion a, Quaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) < 0.01f;
        public static bool operator ==(Quaternion a, SQuaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) < 0.01f;
        public static bool operator !=(SQuaternion a, SQuaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) >= 0.01f;
        public static bool operator !=(SQuaternion a, Quaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) >= 0.01f;
        public static bool operator !=(Quaternion a, SQuaternion b)
            => 1f - Mathf.Abs(Quaternion.Dot(a, b)) >= 0.01f;
    }
 
    /// <summary> Serializable version of UnityEngine.Color32 without transparency. </summary>
    [System.Serializable]
    public struct SColor32 {
        public byte r;
        public byte g;
        public byte b;
 
        public SColor32(byte r, byte g, byte b) {
            this.r = r;
            this.g = g;
            this.b = b;
        }
 
        public SColor32(Color32 c) {
            r = c.r;
            g = c.g;
            b = c.b;
        }
 
        public override string ToString()
            => $"[{r}, {g}, {b}]";
 
        public static implicit operator Color32(SColor32 rValue)
            => new Color32(rValue.r, rValue.g, rValue.b, a: byte.MaxValue);
 
        public static implicit operator SColor32(Color32 rValue)
            => new SColor32(rValue.r, rValue.g, rValue.b);
    }

    [System.Serializable]
    public struct SRaycastHit {
        // The ones we generally want to care about
        public float distance;
        public SVector3 normal;
        public SVector3 point;
        // The things we generally don't care about
        public SVector3 barycentricCoordinate;
        public int colliderInstanceID;
        public SVector2 lightmapCoord;
        public SVector2 textureCoord;
        public SVector2 textureCoord2;
        public int triangleIndex;

        public SRaycastHit(float distance, SVector3 normal, SVector3 point) {
            this.distance = distance;
            this.normal = normal;
            this.point = point;
            this.barycentricCoordinate = default(SVector3);
            this.colliderInstanceID = default(int);
            this.lightmapCoord = default(SVector2);
            this.textureCoord = default(SVector2);
            this.textureCoord2 = default(SVector2);
            this.triangleIndex = default(int);
        }
        public SRaycastHit(
            float distance, SVector3 normal, SVector3 point, 
            SVector3 barycentricCoordinate, int colliderInstanceID,
            SVector2 lightmapCoord, 
            SVector2 textureCoord, SVector2 textureCoord2,
            int triangleIndex
        ) {
            this.distance = distance;   this.normal = normal;   this.point = point;
            this.barycentricCoordinate = barycentricCoordinate;
            this.colliderInstanceID = colliderInstanceID;
            this.lightmapCoord = lightmapCoord;
            this.textureCoord = textureCoord;
            this.textureCoord2 = textureCoord2;
            this.triangleIndex = triangleIndex;
        }

        public override string ToString()
            => $"Distance: {this.distance}\nNormal: {this.normal}\nPoint: {this.point}]";
 
        public static implicit operator SRaycastHit(RaycastHit hit)
            => new SRaycastHit(
                hit.distance, hit.normal, hit.point,
                hit.barycentricCoordinate, hit.colliderInstanceID,
                hit.lightmapCoord,
                hit.textureCoord, hit.textureCoord2,
                hit.triangleIndex
            );
    }

}
