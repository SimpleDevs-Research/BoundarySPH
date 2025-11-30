using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;
using UnityEditor;
using OP = ObstaclePrimitives.Classes;

[ExecuteInEditMode]
public class Grid : MonoBehaviour
{

    [Header("== GRID CONFIGURATIONS ==")]
        #if UNITY_EDITOR
        [Help("These configurations will adjust the grid conditions for the simulation.", UnityEditor.MessageType.None)]
        #endif

        [SerializeField, Tooltip("The transform reference for the lower bound of the grid")]
        private Transform _LOWER_BOUND_TRANSFORM = null;
        [SerializeField, Tooltip("The transform reference for the upper bound of the grid")]
        private Transform _UPPER_BOUND_TRANSFORM = null;
        [SerializeField, Tooltip("How big are our grid cells? For Prefix Summation (no longer used), recommended to match `smoothingRadius`")] 
        private float _gridCellSize = 2f;
        public float gridCellSize => _gridCellSize;
        [SerializeField, Tooltip("In our grid space, how many buffer cells are added to each axis? Limited to 2 at a minimum, to add one cell on each end of the axis.")]
        private int _bufferCellsPerAxis = 2;
        [ReadOnlyInsp, SerializeField, Tooltip("Where in world space are we centering the simulation around?")]
        private Vector3 _origin = new Vector3(0f,0f,0f);
        public Vector3 origin => _origin;
        public float[] originF => new float[3]{origin.x, origin.y, origin.z};
        [ReadOnlyInsp, SerializeField, Tooltip("What are the world space length (per axis) is the simulation, limited by buffer cells?")]
        private float[] _innerBounds = new float[6]{0f, 0f, 0f, 0f, 0f, 0f};
        public float[] innerBounds => _innerBounds;
        public Vector3 innerBoundsV3 => new Vector3(_innerBounds[3]-_innerBounds[0], _innerBounds[4]-_innerBounds[1], _innerBounds[5]-_innerBounds[2]);
        private float[] innerBoundsF => new float[3]{_innerBounds[3]-_innerBounds[0], _innerBounds[4]-_innerBounds[1], _innerBounds[5]-_innerBounds[2]};
        [ReadOnlyInsp, SerializeField, Tooltip("What are the total world space length (per axis) is the simulation?")]
        private float[] _outerBounds = new float[6]{0f, 0f, 0f, 0f, 0f, 0f};
        public float[] outerBounds => _outerBounds;
        private Vector3 outerBoundsV3 => new Vector3(_outerBounds[3]-_outerBounds[0], _outerBounds[4]-_outerBounds[1], _outerBounds[5]-_outerBounds[2]);
        private float[] outerBoundsF => new float[3]{_outerBounds[3]-_outerBounds[0], _outerBounds[4]-_outerBounds[1], _outerBounds[5]-_outerBounds[2]};
        [ReadOnlyInsp, SerializeField, Tooltip("Given `_outerBounds`, how many grid cells are along each axis?")]
        private int[] _numCellsPerAxis = new int[3]{0,0,0};
        public int[] numCellsPerAxis => _numCellsPerAxis;
        public Vector3Int numCellsPerAxisV3I => new Vector3Int(_numCellsPerAxis[0], _numCellsPerAxis[1], _numCellsPerAxis[2]);
        public float[] numCellsPerAxisF => new float[3]{(float)_numCellsPerAxis[0], (float)_numCellsPerAxis[1], (float)_numCellsPerAxis[2]};
        // How many grid cells do we have in total?
        public int numGridCells => _numCellsPerAxis[0] * _numCellsPerAxis[1] * _numCellsPerAxis[2];
        [SerializeField, ReadOnlyInsp] private float[] _gridScaling = new float[3]{1f,1f,1f};
        public float[] gridScaling => _gridScaling;

    
    // Helper function called while NOT IN PLAY MODE (aka during editing mode).
    // This automatically calculates anything grid-related, such as what the boundaries for the simulation are, how many grid cells we have, etc.
    // There are some key concepts to understand first:
    //  1) WE don't set the origin - the SYSTEM does. The system calculates the origin based on the placement of two GameObjects in the unity scene.
    //  2) The space defined by the two GameObjects is actually the outermost possible boundary. The simulation and all its grid cells are placed inside this space
    //  3) We have an `inner` and an `outer` bounds. The `inner` bounds is essentially where all particles will remain. The `outer` bounds is buffer cell space so that the particles don't accidentally look at a nonexistent cell when looking in neighbor cells
    private void UpdateBounds() {
        if (_LOWER_BOUND_TRANSFORM == null || _UPPER_BOUND_TRANSFORM == null) return;
        // First step: Determine our intended origin point from `_LOWER_BOUND_TRANSFORM` and `_UPPER_BOUND_TRANSFORM`
        _origin = (_LOWER_BOUND_TRANSFORM.position + _UPPER_BOUND_TRANSFORM.position) / 2f;
        // Now, we need to calculate the bounds out from origin, based on gridcellsize
        // To do that, we calculate the 3D bounds by doing the following:
        Vector3 space = _UPPER_BOUND_TRANSFORM.position - _LOWER_BOUND_TRANSFORM.position;
        // Then, for each dimension (x, y, and z), we calculate how many cells we can ideally fit inside that.
        // AKA: Calculate how many cells will fit within the provided dimensions
        _numCellsPerAxis = new int[3] {
            Mathf.FloorToInt(space.x / _gridCellSize),
            Mathf.FloorToInt(space.y / _gridCellSize),
            Mathf.FloorToInt(space.z / _gridCellSize)
        };
        // The new bounds are defined, from the origin, how many cells on the X, Y, and Z axes we can fit
        // We define two different versions: a `bounds` that doesn't include the buffer cells, and and `outerBounds` that does consider buffer cells.
        // Note that the space encased by `_UPPER/_LOWER_BOUND_TRANSFORM` contains `outerBounds
        //if (_bufferCellsPerAxis < 2) _bufferCellsPerAxis = 2;
        //else if (_bufferCellsPerAxis % 2 == 1) _bufferCellsPerAxis = Mathf.Max(2,Mathf.RoundToInt(_bufferCellsPerAxis/2));
        if (_bufferCellsPerAxis % 2 == 1) _bufferCellsPerAxis = Mathf.Max(2,Mathf.RoundToInt(_bufferCellsPerAxis/2));
        Vector3 innerBoundsHalf = new Vector3(
            (_numCellsPerAxis[0] - _bufferCellsPerAxis) * _gridCellSize,
            (_numCellsPerAxis[1] - _bufferCellsPerAxis) * _gridCellSize,
            (_numCellsPerAxis[2] - _bufferCellsPerAxis) * _gridCellSize
        ) / 2f;
        _innerBounds = new float[6] {
            _origin.x - innerBoundsHalf.x,
            _origin.y - innerBoundsHalf.y,
            _origin.z - innerBoundsHalf.z,
            _origin.x + innerBoundsHalf.x,
            _origin.y + innerBoundsHalf.y,
            _origin.z + innerBoundsHalf.z
        };
        Vector3 outerBoundsHalf = new Vector3(
            _numCellsPerAxis[0] * _gridCellSize, 
            _numCellsPerAxis[1] * _gridCellSize, 
            _numCellsPerAxis[2] * _gridCellSize
        ) / 2f;
        _outerBounds = new float[6] {
            _origin.x - outerBoundsHalf.x,
            _origin.y - outerBoundsHalf.y,
            _origin.z - outerBoundsHalf.z,
            _origin.x + outerBoundsHalf.x,
            _origin.y + outerBoundsHalf.y,
            _origin.z + outerBoundsHalf.z
        };
        _gridScaling = new float[3]{
            _origin.x - (_numCellsPerAxis[0] * _gridCellSize)/2f,
            _origin.y - (_numCellsPerAxis[1] * _gridCellSize)/2f,
            _origin.z - (_numCellsPerAxis[2] * _gridCellSize)/2f
        };
    }

    [Tooltip("Create subsections of this grid for other particle systems. Subsections adhere to the inner bounds of the grid and will adopt its grid cell size and other properties.")]
    public List<OP.ParticleGridSection> sections = new List<OP.ParticleGridSection>();

    [Header("== PARTICLE CONFIGURATIONS ==")]
        #if UNITY_EDITOR
        [Help("These configurations adjust the particles themselves, such as how many particles are present and how big they are in world scale.", UnityEditor.MessageType.None)]
        #endif
        
        [Tooltip("The mesh used to render each particle in the simulation. Usually just the default `Sphere` mesh from Unity.")]
        public Mesh particle_mesh;
        public Material particle_material;

        [Tooltip("The mesh and shader used to render each grid cell in the simulation. Usually the `Cube` mesh from Unity.")]
        public Mesh grid_cell_mesh;
        public Material grid_cell_material;


    [SerializeField]
    private bool showGizmos = true;
    [SerializeField]
    private List<ParticleController> _controllers = new List<ParticleController>();

     void OnDrawGizmos() {
        if (!showGizmos) return;

        if (_LOWER_BOUND_TRANSFORM == null || _UPPER_BOUND_TRANSFORM == null) return;
        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(origin, 0.1f);
        Gizmos.DrawWireCube(origin, innerBoundsV3);
        Gizmos.color = Color.white;
        
        if (Application.isPlaying) return;
        
        Gizmos.DrawWireCube(origin, outerBoundsV3);
        Gizmos.color = new Vector4(1f,1f,1f,0.5f);
        Gizmos.DrawWireCube(origin, _UPPER_BOUND_TRANSFORM.position - _LOWER_BOUND_TRANSFORM.position);

        if (sections.Count > 0) {
            foreach(OP.ParticleGridSection section in sections) {
                Gizmos.color = section.gizmosColor;
                Gizmos.DrawCube(section.origin, section.dimensionsV3);
            }
        }

        /*
        Gizmos.color = new Vector4(0f,0f,1f,0.75f);
        float waterHeight = Mathf.FloorToInt((_WATER_LEVEL_TRANSFORM.position.y-_innerBounds[1])/gridCellSize) * gridCellSize;
        Vector3 waterDisplayCenter = new Vector3(origin.x, _innerBounds[1] +  waterHeight * 0.5f, origin.z);
        Gizmos.DrawCube(waterDisplayCenter, new Vector3(innerBounds.x, waterHeight, innerBounds.z));

        Vector3 handlePos = _WATER_LEVEL_TRANSFORM.position + Vector3.left * 20f;
        Handles.Label(handlePos, $"Grid Cell Dimensions: ({_numCellsPerAxis[0]},{_numCellsPerAxis[1]},{_numCellsPerAxis[2]})\n# Particles: {numParticles}");
        */
    }

    [SerializeField, ReadOnlyInsp]
    private bool started = false;

    private void Start() {
        if (Application.isPlaying && !started) {
            // Have to deactivate all bound transforms
            if (_LOWER_BOUND_TRANSFORM != null) _LOWER_BOUND_TRANSFORM.gameObject.SetActive(false);
            if (_UPPER_BOUND_TRANSFORM != null) _UPPER_BOUND_TRANSFORM.gameObject.SetActive(false);
            if (sections.Count > 0) {
                foreach(OP.ParticleGridSection section in sections) {
                    if (section.LOWER_BOUND_REF != null) section.LOWER_BOUND_REF.gameObject.SetActive(false);
                    if (section.UPPER_BOUND_REF != null) section.UPPER_BOUND_REF.gameObject.SetActive(false);
                }
            }
            started = true; 
        } else if (started) {
            // Have to activate all bound transforms
            if (_LOWER_BOUND_TRANSFORM != null) _LOWER_BOUND_TRANSFORM.gameObject.SetActive(true);
            if (_UPPER_BOUND_TRANSFORM != null) _UPPER_BOUND_TRANSFORM.gameObject.SetActive(true);
            if (sections.Count > 0) {
                foreach(OP.ParticleGridSection section in sections) {
                    if (section.LOWER_BOUND_REF != null) section.LOWER_BOUND_REF.gameObject.SetActive(true);
                    if (section.UPPER_BOUND_REF != null) section.UPPER_BOUND_REF.gameObject.SetActive(true);
                }
            }
            started = false;
        }
    }

    // Update() is only called when something in the scene changes
    private void Update() {
        if (Application.isPlaying) return;
        UpdateBounds();
        DetectControllers();
        if (sections.Count > 0) {
            foreach(OP.ParticleGridSection section in sections) {
                section.UpdateSegment(_innerBounds, _numCellsPerAxis, _gridCellSize);
            }
        }
    }

    private void DetectControllers() {
        Component[] cs = GetComponentsInChildren<ParticleController>();
        _controllers = new List<ParticleController>();
        foreach(ParticleController c in cs) {
            c._GRID = this;
            c.CalculateParticles();
            _controllers.Add(c);
        }
    }

    public int3 GetGridXYZIndices(Vector3 position) {
        return new(
            Mathf.FloorToInt((position.x - (_origin.x - (_numCellsPerAxis[0] * _gridCellSize)/2f))/_gridCellSize),
            Mathf.FloorToInt((position.y - (_origin.y - (_numCellsPerAxis[1] * _gridCellSize)/2f))/_gridCellSize),
            Mathf.FloorToInt((position.z - (_origin.z - (_numCellsPerAxis[2] * _gridCellSize)/2f))/_gridCellSize)
        );
    }
    public int3 GetGridXYZIndices(float3 position) {
        return new(
            Mathf.FloorToInt((position[0] - (_origin.x - (_numCellsPerAxis[0] * _gridCellSize)/2f))/_gridCellSize),
            Mathf.FloorToInt((position[1] - (_origin.y - (_numCellsPerAxis[1] * _gridCellSize)/2f))/_gridCellSize),
            Mathf.FloorToInt((position[2] - (_origin.z - (_numCellsPerAxis[2] * _gridCellSize)/2f))/_gridCellSize)
        );
    }

    public int GetProjectedGridIndexFromXYZ(int x, int y, int z) {
        return x + (_numCellsPerAxis[0] * y) + (_numCellsPerAxis[0] * _numCellsPerAxis[1] * z);
    }

    public float3 GetWorldPositionFromXYZ(int x, int y, int z) {
        return new(
            _outerBounds[0] + (_gridCellSize * 0.5f) + (_gridCellSize * x),
            _outerBounds[1] + (_gridCellSize * 0.5f) + (_gridCellSize * y),
            _outerBounds[2] + (_gridCellSize * 0.5f) + (_gridCellSize * z)
        );
    }
}
