using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;

public class CPU_Obstacle : SPH_Obstacle
{
    [Header("== MESH VARIABLES ==")]
    // The boundaries of the obstacle in world-scale. Extracted from the mesh renderer, not from any vertex calculation
    [SerializeField, ReadOnlyInsp] private float[] _bounds;
     // The world-scale vertices of the mesh.
    [SerializeField, ReadOnlyInsp] private float3[] _vertices;
    // List to store the triangles of the simulation.
    [SerializeField, ReadOnlyInsp] private ParticleTriangle[] _triangles;

    [Header("== DEBUGGING ==")]
    // Counters for the intersection calculation for each debug particle. If the counter > 0, then the partice is inside the mesh
    [SerializeField, ReadOnlyInsp] private int[] _counters;
    public int[] counters => _counters;

    // DEBUG ELEMENT: stores all projections of the current debug particle so that we can render them via Gizmos
    private List<Vector3> _projections = new List<Vector3>();
    // DEBUG ELEMENT: stores the closest projection and is used for the raycast direction.
    private Vector3[] _closestPoints = new Vector3[0];
    private bool[] _isIntersecting = new bool[0];

    void OnDrawGizmos() {
        if (_bounds.Length != 6) return;
        Vector3 b = new Vector3(_bounds[3]-_bounds[0], _bounds[4]-_bounds[1], _bounds[5]-_bounds[2]);
        Vector3 c = new Vector3(_bounds[3]+_bounds[0], _bounds[4]+_bounds[1], _bounds[5]+_bounds[2]) / 2f;
        Gizmos.color = Color.yellow;
        Gizmos.DrawWireCube(c,b);

        if (_debugParticles.Count > 0) {
            for(int i = 0; i < _debugParticles.Count; i++) {
                Vector3 pos = _debugParticles[i].position;
                Gizmos.color = Color.blue;
                Gizmos.DrawRay(pos, _closestPoints[i] - pos);

            }

            if (_projections.Count > 0) {
                Gizmos.color = Color.yellow;
                foreach(Vector3 p in _projections) {
                    Gizmos.DrawSphere(p,0.05f);
                }
            }
        }
    }

    void Start() {
        UpdateWorlScaleVariables();
        
        _closestPoints = new Vector3[_debugParticles.Count];
        for(int i = 0; i < _debugParticles.Count; i++) {
            _closestPoints[i] = _debugParticles[i].position;
        }

        transform.hasChanged = false;
    }

    private void UpdateWorlScaleVariables() {
        // Need to update the bounds from the mesh
        base.UpdateBounds(out _bounds);

        // Need to update the triangles and vertices from the mesh
        base.CalculateTrianglesAndVertices(out _vertices, out _triangles);
    }

    void Update() {
        _projections = new List<Vector3>();
        _closestPoints = new Vector3[_debugParticles.Count];
        _counters = new int[_debugParticles.Count];
        _isIntersecting = new bool[_debugParticles.Count];

        Vector3 query;
        int closeTriangle;
        
        for(int i = 0; i < _debugParticles.Count; i++) {
            query = _debugParticles[i].position;
            if (CheckIfInBounds(query)) {
                _closestPoints[i] = FindClosestPoint(query, out closeTriangle);
                _isIntersecting[i] = CheckIfIntersecting(query, _closestPoints[i], out _counters[i]);
            } else {
                _closestPoints[i] = query;
                _isIntersecting[i] = false;
            }
        }
    }

    private bool CheckIfInBounds(Vector3 query) {
        return (
            _bounds[0] <= query.x 
            && _bounds[1] <= query.y
            && _bounds[2] <= query.z
            && _bounds[3] >= query.x
            && _bounds[4] >= query.y
            && _bounds[5] >= query.z
        );
    }

    private Vector3 FindClosestPoint(Vector3 query, out int closestTriangleIndex) {
        Vector3 v1, v2, v3, c, n;
        Vector3 projection;

        Vector3 closestPoint = query;
        closestTriangleIndex = 0;
        float dist, closestDistance = Mathf.Infinity;

        for(int ti = 0; ti < _triangles.Length; ti++) {
            v1 = _vertices[_triangles[ti].vertexIndices[0]];
            v2 = _vertices[_triangles[ti].vertexIndices[1]];
            v3 = _vertices[_triangles[ti].vertexIndices[2]];
            c = new Vector3(_triangles[ti].c[0], _triangles[ti].c[1], _triangles[ti].c[2]);
            n = new Vector3(_triangles[ti].n[0], _triangles[ti].n[1], _triangles[ti].n[2]).normalized;

            GetRayProjectionOntoPlane(
                query, 
                Mathf.Sign(Vector3.Dot(n, c - query)) * n, 
                n, 
                _vertices[_triangles[ti].vertexIndices[0]], 
                out projection
            );
            if (CheckIfPointInTriangle(projection, v1, v2, v3, n) == 1) {
                dist = Vector3.Distance(query, projection);
                if (dist < closestDistance) {
                    closestDistance = dist;
                    closestPoint = projection;
                    closestTriangleIndex = ti;
                }
            }
        }

        return closestPoint;
    }

    private bool CheckIfIntersecting(Vector3 query, Vector3 closest, out int count) {        
        int isValid;
        Vector3 v1, v2, v3, c, n;
        Vector3 projection;

        int counter = 0;
        Vector4 rot = new Vector4(transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w);

        Vector3 raycast = closest - query;
        
        foreach(ParticleTriangle t in _triangles) {
            v1 = _vertices[t.vertexIndices[0]];
            v2 = _vertices[t.vertexIndices[1]];
            v3 = _vertices[t.vertexIndices[2]];
            c = new Vector3(t.c[0], t.c[1], t.c[2]);
            n = new Vector3(t.n[0], t.n[1], t.n[2]).normalized;

            isValid = GetRayProjectionOntoPlane(query, raycast, n, v1, out projection);
            if (isValid == 0) continue;
        
            if (CheckIfPointInTriangle(projection, v1, v2, v3, n) == 1) {
                counter += (Vector3.Dot(raycast,n) <= 0) ? -1 : 1;
                _projections.Add(projection);
            }
        }
        count = counter;
        return counter > 0;
    }


}
