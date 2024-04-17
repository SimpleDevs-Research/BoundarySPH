using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;

public class MeshObs : MonoBehaviour
{
    [SerializeField] private Transform _position_transform = null;
    public Transform position_transform {
        get { return _position_transform; }
        set {_position_transform = value;}
    }

    [SerializeField] private Transform _mesh_transform = null;
    public Transform mesh_transform {
        get { return _mesh_transform; }
        set { _mesh_transform = value; }
    }

    [SerializeField] private MeshFilter _meshFilter = null;
    public MeshFilter mf {
        get { return _meshFilter; }
        set { _meshFilter = value; }
    }
    public bool hasMeshFilter => _meshFilter != null;

    [SerializeField] private SkinnedMeshRenderer _skinnedMeshRenderer = null;
    public SkinnedMeshRenderer smr {
        get { return _skinnedMeshRenderer; }
        set { _skinnedMeshRenderer = value; }
    }
    public bool hasSkinnedMeshFilter => _skinnedMeshRenderer != null;
    [SerializeField] private Vector3 _skinnedMeshRendererRotation = new Vector3(0f,0f,0f);

    [SerializeField] private Rigidbody _rigidbody = null;
    public Rigidbody rb => _rigidbody;
    public bool hasRigidbody => _rigidbody != null;

    public enum InitializationSettings {
        OnAwake,
        Manual
    }
    [SerializeField] private InitializationSettings _initializationSetting = InitializationSettings.OnAwake;
    private bool initialized = false;

    private Mesh _mesh;
    public Mesh mesh => _mesh;
    [SerializeField] private Vector3[] _vertices;
    [SerializeField] private int[] _triangles;
    public Vector3[] vertices => _vertices;
    public int[] triangles => _triangles;
    private Dictionary<Vector3,int> condensed_vertices = new Dictionary<Vector3,int>();
    [HideInInspector] public List<float3> fixed_vs;
    [HideInInspector] public uint[] vs_map;

    private void Awake() {
        if (_initializationSetting == InitializationSettings.OnAwake || !initialized) Initialize();
    }

    private void Start() {
        //if (transform.GetComponent<Animation>() != null) transform.GetComponent<Animation>().Play();
    }

    public void ManualInitialization() {
        Initialize();
        _initializationSetting = InitializationSettings.Manual;
    }

    public void Initialize() {
        if (_position_transform == null) _position_transform = this.transform;
        
        if (_mesh_transform == null) _mesh_transform = _position_transform;

        if (_meshFilter == null) {
            MeshFilter m = _mesh_transform.GetComponent<MeshFilter>();
            if (m != null) _meshFilter = m;
        }
        if (_skinnedMeshRenderer == null) {
            SkinnedMeshRenderer s = _mesh_transform.GetComponent<SkinnedMeshRenderer>();
            if (s != null) _skinnedMeshRenderer = s;
        }
        if (_rigidbody == null) {
            Rigidbody r = _position_transform.GetComponent<Rigidbody>();
            if (r != null) _rigidbody = r;
        }

        _mesh = GetMesh();
        if (_mesh != null) {
            if (_skinnedMeshRenderer != null) {
                _vertices = new Vector3[_mesh.vertices.Length];
                for(int i = 0; i < _mesh.vertices.Length; i++) {
                     _vertices[i] = RotMultVec3(Quaternion.Euler(_skinnedMeshRendererRotation.x, _skinnedMeshRendererRotation.y, _skinnedMeshRendererRotation.z), _mesh.vertices[i]);
                }
            } 
            else {
                // Should be local space meshes...
                _vertices = _mesh.vertices;
            }
            Debug.Log(_vertices.Length);
            _triangles = _mesh.triangles;
        }

        initialized = true;
    }

    public Mesh GetMesh() {
        if (hasSkinnedMeshFilter) return _skinnedMeshRenderer.sharedMesh;
        if (hasMeshFilter) return _meshFilter.sharedMesh;
        return null;
    }

    public float3[] GetWorldVertices() {
        if (!hasSkinnedMeshFilter) return null;

        Mesh temp = new Mesh();
        _skinnedMeshRenderer.BakeMesh(temp);
        
        List<uint> added_vertex_indices = new List<uint>();
        List<float3> toReturn = new List<float3>();
        for(int i = 0; i < temp.vertices.Length; i++) {
            if (!added_vertex_indices.Contains(vs_map[i])) {
                added_vertex_indices.Add(vs_map[i]);
                Vector3 toRotate = (_mesh_transform.rotation * temp.vertices[i]) + _mesh_transform.position;
                toReturn.Add(new(toRotate.x, toRotate.y, toRotate.z));
            }
        }
        return toReturn.ToArray();
    }

    public static Vector3 RotMultVec3(Quaternion quat, Vector3 vec){
        float num = quat.x * 2f;
        float num2 = quat.y * 2f;
        float num3 = quat.z * 2f;
        float num4 = quat.x * num;
        float num5 = quat.y * num2;
        float num6 = quat.z * num3;
        float num7 = quat.x * num2;
        float num8 = quat.x * num3;
        float num9 = quat.y * num3;
        float num10 = quat.w * num;
        float num11 = quat.w * num2;
        float num12 = quat.w * num3;

        float x = (1f - (num5 + num6)) * vec.x + (num7 - num12) * vec.y + (num8 + num11) * vec.z;
        if (Mathf.Abs(x) < 0.000000000001f) x = 0f;
        float y = (num7 + num12) * vec.x + (1f - (num4 + num6)) * vec.y + (num9 - num10) * vec.z;
        if (Mathf.Abs(y) < 0.000000000001f) y = 0f;
        float z = (num8 - num11) * vec.x + (num9 + num10) * vec.y + (1f - (num4 + num5)) * vec.z;
        if (Mathf.Abs(z) < 0.000000000001f) z = 0f;
        return new Vector3(x,y,z);
    }

}
