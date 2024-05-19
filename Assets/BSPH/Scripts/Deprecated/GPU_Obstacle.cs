using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using System.Linq;
using OP = ObstaclePrimitives.Structs;

public class GPU_Obstacle : SPH_Obstacle
{
    public BufferManager _BM;

    // This is for GPU Threading!
    private Vector3Int _NUM_THREADS = new Vector3Int(256,4,1);
    private int _NUM_BLOCKS_PARTICLES, _NUM_BLOCKS_TRIANGLES;

    [Header("== GPU VARIABLES ==")]
    [SerializeField] private ParticleController _particleController = null;
    // We need a compute shader for the GPU calculations. Make sure to set one up and reference it here.
    [SerializeField, Tooltip("We need a compute shader for the GPU calculations. Make sure to set one up and reference it here.")] 
    private ComputeShader _obstacleShader = null;
    
    private int _numParticles = -1, _numTriangles = -1, _numVertices = -1;
    [SerializeField] private int[] _counters;
    [SerializeField] private OP.Projection[] _projections;
    //[SerializeField] private int[] _projection_counts;

    [Header("== GIZMOS SETTINGS ==")]
    [SerializeField] private bool _showBounds = true;
    [SerializeField] private bool _showProjections = false;

    void OnDrawGizmos() {
        if (!Application.isPlaying || !this.enabled) return;

        if (_showBounds) {
            float[] bs = new float[6];
            _BOUNDS_BUFFER.GetData(bs);
            Vector3 bsf = new Vector3(bs[3] - bs[0], bs[4]-bs[1], bs[5]-bs[2]);
            Vector3 cen = new Vector3(bs[3] + bs[0], bs[4] + bs[1], bs[5] + bs[2]) / 2f;
            Gizmos.color = Color.yellow;
            Gizmos.DrawWireCube(cen, bsf);
        }
        
        if (!_showProjections) return;
        
        OP.Particle[] ps = new OP.Particle[_numParticles];
        _BM.PARTICLES_BUFFER.GetData(ps);
        float3[] vs = new float3[_numVertices];
        _VERTICES_BUFFER.GetData(vs);
        ParticleTriangle[] ts = new ParticleTriangle[_numTriangles];
        _TRIANGLES_BUFFER.GetData(ts);
        
        for(int i = 0; i < _numParticles; i++) {
            //if (i >= _projections.Length || _projections[i].triangleIndex == -1) continue;
            Vector3 pos = ps[i].position;
            Vector3 proj = new Vector3(_projections[i].position[0], _projections[i].position[1], _projections[i].position[2]);
            int triIndex = (int)_projections[i].triangleID;
            Vector3 ray = proj - pos;

            Gizmos.color = Color.blue;
            Gizmos.DrawRay(pos, ray);

            Gizmos.color = Color.green;
            Gizmos.DrawSphere(proj, 0.025f);

            if (triIndex != -1) {
                Gizmos.color = Color.red;
                ParticleTriangle tri = ts[triIndex];
                float3 v1 = vs[tri.vertexIndices[0]];
                Gizmos.DrawSphere(v1,0.1f);
                float3 v2 = vs[tri.vertexIndices[1]];
                Gizmos.DrawSphere(v2,0.1f);
                float3 v3 = vs[tri.vertexIndices[2]];
                Gizmos.DrawSphere(v3,0.1f);
            }
        }
        /*
        Gizmos.color = Color.red;
        for(int v = 0; v < _numVertices; v++) {
            Gizmos.DrawSphere(vs[v], 0.15f);
        }
        */

    }
    
    void Start() {
        // We call on update to the shader
        InitializeShader();
    }

    private void InitializeShader() {
        // We instantiate a local copy of the compute shader
        _obstacleShader = Instantiate(_obstacleShader);
        // We initialize the kernels first and foremost!
        InitializeKernels();
        // We update the buffer to store particles inside.
        // This function also updates any buffers that are tied to particle count
        UpdateParticles();
        // Initialize the boundaries formed by the mesh. We also simultaneously initialize `cellLimits`
        InitializeBounds();
        // Update the vertices and triangles
        InitializeVertsAndTriangles();
        // Finally, we update this obstacle's world-scale position, rotation, and scale
        UpdateTransform();
    }

    private int _CHECK_IN_BOUNDS_KERNEL;
    private int _FIND_CLOSEST_POINT_KERNEL;
    private int _CLEAR_COUNTERS_KERNEL;
    //private int _CHECK_INTERSECTIONS_1_KERNEL;
    private int _CHECK_INTERSECTIONS_2_KERNEL;
    private int _PUSH_PARTICLES_KERNEL;
    private void InitializeKernels() {
        _CHECK_IN_BOUNDS_KERNEL = _obstacleShader.FindKernel("CheckInBounds");
        _FIND_CLOSEST_POINT_KERNEL = _obstacleShader.FindKernel("FindClosestPoint");
        _CLEAR_COUNTERS_KERNEL = _obstacleShader.FindKernel("ClearCounters");
        //_CHECK_INTERSECTIONS_1_KERNEL = _obstacleShader.FindKernel("CheckIntersections1");
        _CHECK_INTERSECTIONS_2_KERNEL = _obstacleShader.FindKernel("CheckIntersections2");
        _PUSH_PARTICLES_KERNEL = _obstacleShader.FindKernel("PushParticles");
    }
    
    // THIS IS PURELY FOR DEBUGGING PURPOSES
    private void UpdateParticles() {
        // We need to check if we have any debug particles or not.
        // If we do, that means we'll be using that array as the basis for our particles.
        // Otherwise, we'll be using the particle buffer offered by ParticleController, if it happens to exist.
        if (_debugParticles.Count > 0) {
            if (_numParticles != _debugParticles.Count) {
                // Uh oh, we have a new count for then umber of particles! We need to adjust course!
                // Set the number of particles as a value
                _numParticles = _debugParticles.Count;
                // Initialize the _BM.PARTICLES_BUFFER buffer if it doesn't exist yet
                if (_BM.PARTICLES_BUFFER != null) _BM.PARTICLES_BUFFER.Release();
                if (_BM.PARTICLES_VELOCITIES_BUFFER != null) _BM.PARTICLES_VELOCITIES_BUFFER.Release();
                _BM.PARTICLES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
                _BM.PARTICLES_VELOCITIES_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*3);
                _obstacleShader.SetBuffer(_CLEAR_COUNTERS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_PARTICLE_VELOCITIES", _BM.PARTICLES_VELOCITIES_BUFFER);
                _obstacleShader.SetInt("numParticles",_numParticles);
                _obstacleShader.SetFloat("particleRadius",_particle_radius);
                // Update the thread groups count
                _NUM_BLOCKS_PARTICLES = Mathf.CeilToInt((float)_numParticles / (float)_NUM_THREADS.x);            
                // Update the relevant buffers on their counts as well
                InitializeResults();
            }

            // Update the buffer's data with the most recent version of our particles
            OP.Particle[] ps = new OP.Particle[_numParticles];
            for(int i = 0; i < _numParticles; i++) {
                OP.Particle p = new OP.Particle();
                Vector3 pos  =_debugParticles[i].position;
                p.position = new(pos.x, pos.y, pos.z);
                ps[i] = p;
            }
            _BM.PARTICLES_BUFFER.SetData(ps);
        } else if (_particleController != null) {
            if (_numParticles != _particleController.numParticles) {
                _numParticles = _particleController.numParticles;
                _particle_radius = _particleController.particleRenderRadius;
                _obstacleShader.SetBuffer(_CLEAR_COUNTERS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_PARTICLES", _BM.PARTICLES_BUFFER);
                _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_PARTICLE_VELOCITIES", _BM.PARTICLES_VELOCITIES_BUFFER);
                _obstacleShader.SetInt("numParticles",_numParticles);
                _obstacleShader.SetFloat("particleRadius", _particle_radius);
                _NUM_BLOCKS_PARTICLES = Mathf.CeilToInt((float)_numParticles / (float)_NUM_THREADS.x);
                InitializeResults();
            }
        }
    }

    // This method initializes any results-based buffers. Primarily for debugging.
    // The following buffers in this case are:
    //  1. _IN_BOUNDS_BUFFER: to take into account the current number of particles
    //  2. _PROJECTIONS_BUFFER: to take into account the direction of the raycast
    //  3. _CLOSEST_TRIANGLE_INDICES_BUFFER: to take into account the closest triangle of a particle to this mesh
    //  4. _COUNTERS_BUFFER: to finally trakc how many intersects a raycast has for each particle towards the direction of the closest triangle mesh
    // This one is updated in the GPU. So if you call this, this is effecitvely re-setting the IN_BOUNDS_BUFFER buffer
    private ComputeBuffer _IN_BOUNDS_BUFFER = null;
    private ComputeBuffer _PROJECTIONS_BUFFER = null;
    //private ComputeBuffer _CLOSEST_TRIANGLE_INDICES_BUFFER = null;
    private ComputeBuffer _COUNTERS_BUFFER = null;

    private void InitializeResults() {
        // If the in_bounds buffer is already initialized, then we have to release it first
        if (_IN_BOUNDS_BUFFER != null) _IN_BOUNDS_BUFFER.Release();
        _IN_BOUNDS_BUFFER = new ComputeBuffer(_numParticles, sizeof(int));
        _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "_IN_BOUNDS", _IN_BOUNDS_BUFFER);
        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_IN_BOUNDS",_IN_BOUNDS_BUFFER);
        _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_IN_BOUNDS", _IN_BOUNDS_BUFFER);

        // If the _PROJECTIONS buffer is already initialized, then we have to release it first
        if (_PROJECTIONS_BUFFER != null) _PROJECTIONS_BUFFER.Release();
        _PROJECTIONS_BUFFER = new ComputeBuffer(_numParticles, sizeof(float)*7 + sizeof(int)*2);
        OP.Projection[] projections = new OP.Projection[_numParticles];
        for(int i = 0; i < _numParticles; i++) {
            OP.Projection p = new OP.Projection();
            //p.triangleID = -1;
            p.position = new(0f,0f,0f);
            //p.d = 0f;
            p.counter = 0;
            //p.check = new(0f,0f,0f);
            projections[i] = p;
        }
        _PROJECTIONS_BUFFER.SetData(projections);
        _obstacleShader.SetBuffer(_CLEAR_COUNTERS_KERNEL, "_PROJECTIONS", _PROJECTIONS_BUFFER);
        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_PROJECTIONS", _PROJECTIONS_BUFFER);
        //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_PROJECTIONS", _PROJECTIONS_BUFFER);
        _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_PROJECTIONS", _PROJECTIONS_BUFFER);
        _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_PROJECTIONS", _PROJECTIONS_BUFFER);

        // If the _COUNTERS buffer is already initialized, then we have to release it first
        if (_COUNTERS_BUFFER != null) _COUNTERS_BUFFER.Release();
        _COUNTERS_BUFFER = new ComputeBuffer(_numParticles, sizeof(int));
        _obstacleShader.SetBuffer(_CLEAR_COUNTERS_KERNEL, "_COUNTERS", _COUNTERS_BUFFER);
        //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_COUNTERS", _COUNTERS_BUFFER);
        _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_COUNTERS", _COUNTERS_BUFFER);
        _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_COUNTERS", _COUNTERS_BUFFER);
        
        _projections = new OP.Projection[_numParticles];
        _counters = new int[_numParticles];
    }


    private ComputeBuffer _BOUNDS_BUFFER = null;
    private void InitializeBounds() {
        // If the bounds buffer is not initialized, do it.
        _BOUNDS_BUFFER = new ComputeBuffer(6, sizeof(float));
        _obstacleShader.SetBuffer(_CHECK_IN_BOUNDS_KERNEL, "_BOUNDS", _BOUNDS_BUFFER);
        // Since this is mere initialization, we also push for updating the bounds
        UpdateBounds();
    }

    // We need a method to update the bounds with the latest info
    // We don't need to worry about size of bounds because the bounds will alweays be 6 floats in length
    private void UpdateBounds() {
        float[] bs = new float[6];
        base.UpdateBounds(out bs);
        _BOUNDS_BUFFER.SetData(bs);
    }

    private ComputeBuffer _VERTICES_BUFFER = null, _TRIANGLES_BUFFER = null;
    private void InitializeVertsAndTriangles() {
        // We use the base class function to get vertices and triangles
        float3[] vertices;
        ParticleTriangle[] triangles;
        base.CalculateTrianglesAndVertices(out vertices, out triangles);
        _numTriangles = triangles.Length;
        _numVertices = vertices.Length;
        // With the vertices and triangles calculated, we need to pass them to the buffers
        // Firstly, release them if they're set already
        if (_VERTICES_BUFFER != null) _VERTICES_BUFFER.Release();
        if (_TRIANGLES_BUFFER != null) _TRIANGLES_BUFFER.Release();
        // Secondly, re-spec the buffers again
        _VERTICES_BUFFER = new ComputeBuffer(_numVertices, sizeof(float)*3);
        _VERTICES_BUFFER.SetData(vertices);
        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_VERTICES", _VERTICES_BUFFER);
        //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_VERTICES", _VERTICES_BUFFER);
        _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_VERTICES", _VERTICES_BUFFER);

        _TRIANGLES_BUFFER = new ComputeBuffer(_numTriangles, sizeof(int)*3 + sizeof(float)*7);
        _TRIANGLES_BUFFER.SetData(triangles);
        _obstacleShader.SetBuffer(_FIND_CLOSEST_POINT_KERNEL, "_TRIANGLES", _TRIANGLES_BUFFER);
        //_obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_1_KERNEL, "_TRIANGLES", _TRIANGLES_BUFFER);
        _obstacleShader.SetBuffer(_CHECK_INTERSECTIONS_2_KERNEL, "_TRIANGLES", _TRIANGLES_BUFFER);
        _obstacleShader.SetBuffer(_PUSH_PARTICLES_KERNEL, "_TRIANGLES", _TRIANGLES_BUFFER);
        // Thirdly, update `numTriangles` and `numVertices` in the shader
        _obstacleShader.SetInt("numVertices", _numVertices);
        _obstacleShader.SetInt("numTriangles", _numTriangles);
        
        // Update _NUM_BLOCKS_TRIANGLES
        _NUM_BLOCKS_TRIANGLES = Mathf.CeilToInt((float)_numTriangles / (float)_NUM_THREADS.y);

    }

    private void UpdateTransform() {
        _obstacleShader.SetFloats("position", new float[3]{transform.position.x, transform.position.y, transform.position.z});
        _obstacleShader.SetFloats("rotation", new float[4]{transform.rotation.x, transform.rotation.y, transform.rotation.z, transform.rotation.w});
        _obstacleShader.SetFloats("scale", new float[3]{transform.lossyScale.x, transform.lossyScale.y, transform.lossyScale.z});
    }

    void Update() {
        // At this poitn, we must have a compute buffer aside for particles. If this doesn't work, then... welp.
        if (_BM.PARTICLES_BUFFER == null) {
            Debug.Log("PARTICLE BUFFER NULL!");
            return;
        }
        // This is immediately done, no questions asked
        UpdateParticles();
        // If the transform has changed, we need to update 1) bounds, 2) vertices, and 3) triangles
        if (_obstacleType == ObstacleType.Dynamic && transform.hasChanged) {
            UpdateBounds();
            InitializeVertsAndTriangles();
            UpdateTransform();
            transform.hasChanged = false;
        }
        // Now perform the necessary dispatches
        // Firstly, clear the counters and projections
        _obstacleShader.Dispatch(_CLEAR_COUNTERS_KERNEL, _NUM_BLOCKS_PARTICLES,1,1);
        // Secondly, check which particles are currently in bounds
        
        _obstacleShader.Dispatch(_CHECK_IN_BOUNDS_KERNEL, _NUM_BLOCKS_PARTICLES,1,1);
        // Thirdly, Find the projections of each particle onto the mesh
        _obstacleShader.Dispatch(_FIND_CLOSEST_POINT_KERNEL, _NUM_BLOCKS_PARTICLES,1,1);
        // Fourthly, Get the number of intersections
        // _obstacleShader.Dispatch(_CHECK_INTERSECTIONS_1_KERNEL, _NUM_BLOCKS_PARTICLES,1,1);
        _obstacleShader.Dispatch(_CHECK_INTERSECTIONS_2_KERNEL, _NUM_BLOCKS_PARTICLES,_NUM_BLOCKS_TRIANGLES,1);

        // Finally, alter the velocities and positions of particles if they are intersecting. We push objects to their projection point and set their velocities to -0.5 of their current velocity
        _obstacleShader.Dispatch(_PUSH_PARTICLES_KERNEL, _NUM_BLOCKS_PARTICLES,1,1);

        // For debugging, we'll update the _counters array we set up as a global variable
        OP.Projection[] projs = new OP.Projection[_numParticles];
        _PROJECTIONS_BUFFER.GetData(projs);
        // _projections = projs.Where(p => p.count == 1).ToArray();
        _projections = projs;
        int[] counts = new int[_numParticles];
        _COUNTERS_BUFFER.GetData(counts);
        _counters = counts.Where(c => c == 0).ToArray();
        
    }
    
    void OnDestroy() {
        if (_BOUNDS_BUFFER != null) {
            _BOUNDS_BUFFER.Release();
            _BOUNDS_BUFFER = null;
        }
        //_CLOSEST_TRIANGLE_INDICES_BUFFER.Release();
        if (_COUNTERS_BUFFER != null) {
            _COUNTERS_BUFFER.Release();
            _COUNTERS_BUFFER = null;
        }
        if (_IN_BOUNDS_BUFFER != null) {
            _IN_BOUNDS_BUFFER.Release();
            _IN_BOUNDS_BUFFER = null;
        }
        if (_PROJECTIONS_BUFFER != null) {
            _PROJECTIONS_BUFFER.Release();
            _PROJECTIONS_BUFFER = null;
        }
        if (_TRIANGLES_BUFFER != null) {
            _TRIANGLES_BUFFER.Release();
            _TRIANGLES_BUFFER = null;
        }
        if (_VERTICES_BUFFER != null) {
            _VERTICES_BUFFER.Release();
            _VERTICES_BUFFER = null;
        }
        
    }

    public void Initialize() {
        float3[] vertices;
        ParticleTriangle[] triangles;
        base.CalculateTrianglesAndVertices(out vertices, out triangles);
        _numTriangles = triangles.Length;
        _numVertices = vertices.Length;
    }
}
