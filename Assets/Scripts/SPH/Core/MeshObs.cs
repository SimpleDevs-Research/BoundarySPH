using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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
    public bool hasSkinnedMeshRenderer => _skinnedMeshRenderer != null;


    [SerializeField] private Rigidbody _rigidbody = null;
    public Rigidbody rb => _rigidbody;
    public bool hasRigidbody => _rigidbody != null;

    private void Awake() {
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
    }

    public Mesh GetMesh() {
        if (hasSkinnedMeshRenderer) return _skinnedMeshRenderer.sharedMesh;
        if (hasMeshFilter) return _meshFilter.sharedMesh;
        return null;
    }

}
