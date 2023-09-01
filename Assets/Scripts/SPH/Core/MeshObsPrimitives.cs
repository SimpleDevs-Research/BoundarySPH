using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Mathematics;


namespace ObstaclePrimitives {

    namespace Classes {

        // =================================== //
        // ===== OBSTACLE REPRESENTATION ===== //
        // =================================== //

        /// <summary>
        /// Represents the STATIC properties of an obstacle. This includes the vertex, triangle, and edge indices.
        /// </summary>
        /// <param name="index">Ref. index inside the obstacles arrays</param>
        /// <param name="vs">[0] = starting index in the vertices arrays, [1] = the number of vertices</param>
        /// <param name="ts">[0] = starting index in the triangles arrays, [1] = the number of triangles</param>
        /// <param name="es">[0] = starting index in the edges arrays, [1] = the number of edges</param>
        /// <param name="mass">The mass of the obstacle. Considered a "static" variable because mass doesn't change.</param>
        [System.Serializable]
        public class ObstacleStatic {
            public int index;
            public int2 vs;
            public int2 ts;
            public int2 es;
            public float mass;
            public uint has_rb;
            public uint has_smr;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of an obstacle that change during the update lop.
        /// </summary>
        /// <param name="index">Ref. index inside the obstacles arrays.</param>
        /// <param name="position">The world-space position of the transform.</param>
        /// <param name="rotation">The world-space rotation of the transform.</param>
        /// <param name="scale">The world-space scale of the transform. Typically via `Transform.LossyScale`.</param>
        /// <param name="lowerBound">The lower boundary position of this obstacle's bounding box.</param>
        /// <param name="upperBound">The upper boundary position of this obstacle's bounding box.</param>
        /// <param name="frictionCoefficient">The amount of friction the object is meant to exhibit on particles.</param>
        /// <param name="checkObstacleBounds">Determines if we should check if the particle is within the bounds of an obstacle or not</param>
        /// <param name="checkTriangleBounds">Determines if we should check if the particle is within the bounds of a triangle or not</param>
        /// <param name="hasChanged">A boolean indicator to determine if this obstacle had transformed in some way between the previous and current frame</param>
        [System.Serializable]
        public class ObstacleDynamic {
            public uint index;
            public uint2 vs;
            public uint2 ts;
            public uint2 es;
            public float3 position;
            public float4 rotation;
            public float3 scale;
            public float density;
            public float3 lowerBound;
            public float3 upperBound;
            public uint isStatic;
            public float3 centerOfMass;
            public float frictionCoefficient;
            public uint checkObstacleBounds;
            public uint checkTriangleBounds;
            public uint hasChanged;
        }


        // ================================= //
        // ===== VERTEX REPRESENTATION ===== //
        // ================================= //

        /// <summary>
        /// Represents the STATIC properties of each vertex. These values do not change at any point after the pre-processing stage.
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="localPosition">The position of the vertex relative to the local space of the mesh. This will never change.</param>
        /// <param name="localNormalEnd">The local endpoint of the 3D normal vector if starting from `localPosition`</param>
        [System.Serializable]
        public class VertexStatic {
            public int obstacleIndex;
            public float3 localPosition;
            public float3 localNormalEnd;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of each vertex. These values have to be updated if the parent obstacle has been transformed between the previous and current frames.
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="position">The current world-space position of this vertex.</param>
        /// <param name="normal">The current world-space, 3D normal vector of this vertex. Direction only.</param>
        /// <param name="force">The calculated force based on the current velocity and previous velocity, which is F=m*(v/t)</param>
        [System.Serializable]
        public class VertexDynamic {
            public int obstacleIndex;
            public float3 position;
            public float3 normal;
            public float3 force;
        }


        // =================================== //
        // ===== TRIANGLE REPRESENTATION ===== //
        // =================================== //

        /// <summary>
        /// Represents the STATIC properties of each triangle. These values never change after the pre-processing stage
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="vertices">The indices of the three vertices of this triangle, referencing values in the vertices arrays. Organized in order of v1-v2-v3.</param>
        /// <param name="angles">The angles of each vertex. Used as weights in the calculation of the 3D normal vectors for each vertex. Organized in order of v1-v2-v3.</param>
        /// <param name="edges">The indices of the three edges of this triangle, referencing values in the edges array. Organized in order of v1v2-v1v3-v2v3.</param>
        /// <param name="localCenter">The local-space centroid of the triangle.</param>
        /// <param name="localNormalEnd">The local endpoint of the 3D normal vector if starting from `localCenter`</param>
        [System.Serializable]
        public class TriangleStatic {
            public int obstacleIndex;
            public int3 vertices;   // v1,v2,v3
            public float3 angles;   // v1,v2,v3
            public int3 edges;      // v1v2,v1v3,v2v3
            public float3 localCenter;
            public float3 localNormalEnd;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of each triangle. These values are updated whenever the parent object is transformed between the previous and current frames.
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="center">The world-space centroid of the triangle, defined by calculating the average of the world-space positions of each vertex.</param>
        /// <param name="normal">The world-space face normal vector of the triangle. Direction only.</param>
        /// <param name="d">The signed distance of the plane to the origin (0,0,0). Needed to calculate the projection points onto the plane defined by this triangle.</param>
        /// <param name="v1v2n">The world-space 2D normal vector of the v1v2 edge.</param>
        /// <param name="v1v3n">The world-space 2D normal vector of the v1v3 edge.</param>
        /// <param name="v2v3n">The world-space 2D normal vector of the v2v3 edge.</param>
        /// <param name="lowerBound">The lower boundary position of this triangle's bounding box.</param>
        /// <param name="upperBound">The upper boundary position of this triangle's bounding box.</param>
        [System.Serializable]
        public class TriangleDynamic {
            public int obstacleIndex;
            public float3 center;
            public float3 normal;
            public float d;
            public uint3 vertices;
            public uint3 edges;
            public float3 v1v2n;
            public float3 v1v3n;
            public float3 v2v3n;
            public float3 lowerBound;
            public float3 upperBound;
        }


        // =============================== //
        // ===== EDGE REPRESENTATION ===== //
        // =============================== //
        
        /// <summary>
        /// Represents the STATIC properties of each edge. These values will never change after the pre-processing stage.
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="vertices">The indices of the two vertices of this edge, referencing values in the vertices arrays. Order is typically ascending in nature.</param>
        /// <param name="triangles">The indices of the two triangles of this edge, referencing values in the triangles arrays. Order is not important.</param>
        /// <param name="midpoint">The middle point on the edge between the two vertices.</param>
        /// <param name="localNormalEnd">The local-space endpoint of the 3D normal vector of the edge, if starting from the local midpoint of the edge line.</param>
        [System.Serializable]
        public class EdgeStatic {
            public int obstacleIndex;
            public int2 vertices;
            public int2 triangles;
            public float3 midpoint;
            public float3 localNormalEnd;
        }


        // =============================== //
        // ===== GRID REPRESENTATION ===== //
        // =============================== //

        [System.Serializable]
        public class ParticleGridSection {

            [SerializeField, Tooltip("The transform reference for the lower bound of this subsection")]
            public Transform LOWER_BOUND_REF = null;
            [SerializeField, Tooltip("The transform reference for the upper bound of this subsection")]
            public Transform UPPER_BOUND_REF = null;
            [SerializeField, ReadOnly, Tooltip("The origin (center) of the section formed from the average of the bound refs")]
            private Vector3 _origin;
            public Vector3 origin => _origin;
            public float[] originF => new float[3]{origin.x, origin.y, origin.z};
            // Store the bounds formed by this section's upper and lower bounds
            [SerializeField, ReadOnly]
            private float[] _bounds = new float[6];
            public float[] bounds => _bounds;
            public Vector3 dimensionsV3 => new Vector3(_bounds[3]-_bounds[0], _bounds[4]-_bounds[1],_bounds[5]-_bounds[2]);
            [SerializeField, Tooltip("Color the bounds of this subsection")]
            private Color _gizmosColor = new Vector4(1f,1f,1f,0.5f);
            public Color gizmosColor => _gizmosColor;

            public void UpdateSegment(float[] restrictionBounds,  int[] numCellsPerAxis, float gridCellSize) {
                if (LOWER_BOUND_REF == null || UPPER_BOUND_REF == null) return;
                // First, make sure that the transforms for the bounds are within the inner bounds of the parent grid
                LOWER_BOUND_REF.position = new Vector3(
                    Mathf.Clamp(
                        LOWER_BOUND_REF.position.x, 
                        restrictionBounds[0], 
                        Mathf.Min(UPPER_BOUND_REF.position.x, restrictionBounds[3])
                    ),
                    Mathf.Clamp(
                        LOWER_BOUND_REF.position.y, 
                    restrictionBounds[1], 
                        Mathf.Min(UPPER_BOUND_REF.position.y, restrictionBounds[4])
                    ),
                    Mathf.Clamp(
                        LOWER_BOUND_REF.position.z, 
                        restrictionBounds[2], 
                        Mathf.Min(UPPER_BOUND_REF.position.z, restrictionBounds[5])
                    )
                );
                UPPER_BOUND_REF.position = new Vector3(
                    Mathf.Clamp(
                        UPPER_BOUND_REF.position.x, 
                        Mathf.Max(LOWER_BOUND_REF.position.x, restrictionBounds[0]), 
                        restrictionBounds[3]
                    ),
                    Mathf.Clamp(
                        UPPER_BOUND_REF.position.y, 
                        Mathf.Max(LOWER_BOUND_REF.position.y, restrictionBounds[1]), 
                        restrictionBounds[4]
                    ),
                    Mathf.Clamp(
                        UPPER_BOUND_REF.position.z, 
                        Mathf.Max(
                            LOWER_BOUND_REF.position.z, 
                            restrictionBounds[2]), 
                            restrictionBounds[5]
                        )
                );
                _origin = (LOWER_BOUND_REF.position + UPPER_BOUND_REF.position) / 2f;
                _bounds = new float[6] {
                    LOWER_BOUND_REF.position.x, LOWER_BOUND_REF.position.y, LOWER_BOUND_REF.position.z,
                    UPPER_BOUND_REF.position.x, UPPER_BOUND_REF.position.y, UPPER_BOUND_REF.position.z
                };
            }
        }


        // ==================================== //
        // ===== PARTICLES REPRESENTATION ===== //
        // ==================================== //

        /// <summary>
        /// Represents the particles of the SPH simulation
        /// Size: sizeof(float)*3
        /// </summary>
        /// <param name="position">The position of the particle in 3D space</param>
        [System.Serializable]
        public class Particle {
            public float3 position;
            public float3 force;
            public int render;
        };

        /// <summary>
        /// Stores some information about the grid cell in the 3D grid. Contains info such as the XYZ index of the cell, the world position, and the lower/upper limits
        /// Size: sizeof(int)*9 + sizeof(float)*3
        /// </summary>
        /// <param name="id">Reference to the XYZ grid indices of this current cell.</param>
        /// <param name="position">Reference to the 3D world position of this grid cell.</param>
        /// <param name="lowerLimits">The lower bounds position of this cell.</param>
        /// <param name="upperLimits">The upper bounds position of this cell.</param>
        [System.Serializable]
        public class CellLimits {
            public int3 id;
            public float3 position;
            public int3 lowerLimits;
            public int3 upperLimits;
        }


        // ================================ //
        // ===== BOIDS REPRESENTATION ===== //
        // ================================ //

        [System.Serializable]
        public class Boid {
            public int obstacleID;
            public float3 position;
            public int3 gridIndices;
            public int projectedGridIndex;
        }

    }

    namespace Structs {

        // =================================== //
        // ===== OBSTACLE REPRESENTATION ===== //
        // =================================== //

        /// <summary>
        /// Represents the STATIC properties of an obstacle. This includes the vertex, triangle, and edge indices.
        /// Size: sizeof(uint)*7 + sizeof(float)*6
        /// </summary>
        /// <param name="index">Ref. index inside the obstacles arrays</param>
        /// <param name="vs">[0] = starting index in the vertices arrays, [1] = the number of vertices</param>
        /// <param name="ts">[0] = starting index in the triangles arrays, [1] = the number of triangles</param>
        /// <param name="es">[0] = starting index in the edges arrays, [1] = the number of edges</param>
        /// <param name="localLowerBound">The local-space lower bound of this obstacle</param>
        /// <param name="localUpperBound">The local-space upper bound of this obstacle</param>
        /// <param name="mass">The mass of the obstacle. Considered a "static" variable because mass doesn't change.</param>
        [System.Serializable]
        public struct ObstacleStatic {
            public uint index;
            public uint2 vs;
            public uint2 ts;
            public uint2 es;
            public float mass;
            public uint has_rb;
            public uint has_smr;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of an obstacle that change during the update lop.
        /// Size: sizeof(uint)*2 + sizeof(float)*17
        /// </summary>
        /// <param name="index">Ref. index inside the obstacles arrays.</param>
        /// <param name="position">The world-space position of the transform.</param>
        /// <param name="rotation">The world-space rotation of the transform.</param>
        /// <param name="scale">The world-space scale of the transform. Typically via `Transform.LossyScale`.</param>
        /// <param name="lowerBound">The lower boundary position of this obstacle's bounding box.</param>
        /// <param name="upperBound">The upper boundary position of this obstacle's bounding box.</param>
        /// <param name="frictionCoefficient">The amount of friction the object is meant to exhibit on particles.</param>
        /// <param name="checkObstacleBounds">Determines if we should check if the particle is within the bounds of an obstacle or not</param>
        /// <param name="checkTriangleBounds">Determines if we should check if the particle is within the bounds of a triangle or not</param>
        /// <param name="hasChanged">A boolean indicator to determine if this obstacle had transformed in some way between the previous and current frame</param>
        [System.Serializable]
        public struct ObstacleDynamic {
            public uint index;
            public uint2 vs;
            public uint2 ts;
            public uint2 es;
            public float3 position;
            public float4 rotation;
            public float3 scale;
            public float density;
            public float3 lowerBound;
            public float3 upperBound;
            public uint isStatic;
            public float3 centerOfMass;
            public float frictionCoefficient;
            public uint checkObstacleBounds;
            public uint checkTriangleBounds;
            public uint hasChanged;
        }


        // ================================= //
        // ===== VERTEX REPRESENTATION ===== //
        // ================================= //

        /// <summary>
        /// Represents the STATIC properties of each vertex. These values do not change at any point after the pre-processing stage.
        /// Size: sizeof(uint) + sizeof(float)*6
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="localPosition">The position of the vertex relative to the local space of the mesh. This will never change.</param>
        /// <param name="localNormal">The local endpoint of the 3D normal vector if starting from `localPosition`</param>
        [System.Serializable]
        public struct VertexStatic {
            public uint obstacleIndex;
            public float3 localPosition;
            public float3 localNormal;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of each vertex. These values have to be updated if the parent obstacle has been transformed between the previous and current frames.
        /// Size: sizeof(uint) + sizeof(float)*6
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="position">The current world-space position of this vertex.</param>
        /// <param name="normal">The current world-space, 3D normal vector of this vertex. Direction only.</param>
        /// <param name="force">The calculated force based on the current velocity and previous velocity, which is F=m*(v/t)</param>
        [System.Serializable]
        public struct VertexDynamic {
            public uint obstacleIndex;
            public float3 position;
            public float3 normal;
            public float3 force;
        }


        // =================================== //
        // ===== TRIANGLE REPRESENTATION ===== //
        // =================================== //

        /// <summary>
        /// Represents the STATIC properties of each triangle. These values never change after the pre-processing stage
        /// Size: sizeof(uint)*7 + sizeof(float)*9
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="vertices">The indices of the three vertices of this triangle, referencing values in the vertices arrays. Organized in order of v1-v2-v3.</param>
        /// <param name="angles">The angles of each vertex. Used as weights in the calculation of the 3D normal vectors for each vertex. Organized in order of v1-v2-v3.</param>
        /// <param name="edges">The indices of the three edges of this triangle, referencing values in the edges array. Organized in order of v1v2-v1v3-v2v3.</param>
        /// <param name="localCenter">The local-space centroid of the triangle.</param>
        /// <param name="localCenter">The local-space centroid of the triangle.</param>
        /// <param name="localNormal">The local endpoint of the 3D normal vector if starting from `localCenter`</param>
        [System.Serializable]
        public struct TriangleStatic {
            public uint obstacleIndex;
            public uint3 vertices;   // v1,v2,v3
            public float3 angles;   // v1,v2,v3
            public uint3 edges;      // v1v2,v1v3,v2v3
            public float3 localCenter;
            public float3 localNormal;
        }

        /// <summary>
        /// Represents the DYNAMIC properties of each triangle. These values are updated whenever the parent object is transformed between the previous and current frames.
        /// Size: sizeof(uint) + sizeof(float)*22
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="center">The world-space centroid of the triangle, defined by calculating the average of the world-space positions of each vertex.</param>
        /// <param name="normal">The world-space face normal vector of the triangle. Direction only.</param>
        /// <param name="d">The signed distance of the plane to the origin (0,0,0). Needed to calculate the projection points onto the plane defined by this triangle.</param>
        /// <param name="v1v2n">The world-space 2D normal vector of the v1v2 edge.</param>
        /// <param name="v1v3n">The world-space 2D normal vector of the v1v3 edge.</param>
        /// <param name="v2v3n">The world-space 2D normal vector of the v2v3 edge.</param>
        /// <param name="lowerBound">The lower boundary position of this triangle's bounding box.</param>
        /// <param name="upperBound">The upper boundary position of this triangle's bounding box.</param>
        [System.Serializable]
        public struct TriangleDynamic {
            public uint obstacleIndex;
            public float3 center;
            public float3 normal;
            public float d;
            public uint3 vertices;
            public uint3 edges;
            public float3 v1v2n;
            public float3 v1v3n;
            public float3 v2v3n;
            public float3 lowerBound;
            public float3 upperBound;
        }


        // =============================== //
        // ===== EDGE REPRESENTATION ===== //
        // =============================== //
        
        /// <summary>
        /// Represents the STATIC properties of each edge. These values will never change after the pre-processing stage.
        /// Size: sizeof(uint)*5 + sizeof(float)*6
        /// </summary>
        /// <param name="obstacleIndex">Ref. index to its parent obstacle in the obstacles arrays.</param>
        /// <param name="vertices">The indices of the two vertices of this edge, referencing values in the vertices arrays. Order is typically ascending in nature.</param>
        /// <param name="triangles">The indices of the two triangles of this edge, referencing values in the triangles arrays. Order is not important.</param>
        /// <param name="midpoint">The middle point on the edge between the two vertices.</param>
        /// <param name="localNormal">The local-space endpoint of the 3D normal vector of the edge, if starting from the local midpoint of the edge line.</param>
        [System.Serializable]
        public struct EdgeStatic {
            public uint obstacleIndex;
            public uint2 vertices;
            public uint2 triangles;
            public float3 midpoint;
            public float3 localNormal;
        }


        // =============================== //
        // ===== GRID REPRESENTATION ===== //
        // =============================== //

        // PARTICLE GRID SECTION (in Classes) DOES NOT HAVE A STRUCT REPRESENTATION
        
        [System.Serializable]
        public struct GridCell {
            public int n;
            public float3 worldPosition;
            public int cellPressure;
        }

        // ==================================== //
        // ===== PARTICLES REPRESENTATION ===== //
        // ==================================== //

        /// <summary>
        /// Represents the particles of the SPH simulation
        /// Size: sizeof(float)*3
        /// </summary>
        /// <param name="position">The position of the particle in 3D space</param>
        [System.Serializable]
        public struct Particle {
            public float3 position;
            public float3 force;
            public int render;
        };

        /// <summary>
        /// Stores some information about the grid cell in the 3D grid. Contains info such as the XYZ index of the cell, the world position, and the lower/upper limits
        /// Size: sizeof(int)*9 + sizeof(float)*3
        /// </summary>
        /// <param name="id">Reference to the XYZ grid indices of this current cell.</param>
        /// <param name="position">Reference to the 3D world position of this grid cell.</param>
        /// <param name="lowerLimits">The lower bounds position of this cell.</param>
        /// <param name="upperLimits">The upper bounds position of this cell.</param>
        [System.Serializable]
        public struct CellLimits {
            public int3 id;
            public float3 position;
            public int3 lowerLimits;
            public int3 upperLimits;
        }
        
        [System.Serializable]
        public struct Projection {
            public uint triangleID;
            public float3 projection;
            public float3 position;
            public float3 normal;
            public float3 particle_force;
            public float3 external_force;
            public int counter;
            public float frictionCoefficient;

            public float3 e1;
            public float3 e2;
            public float3 e1_2DN;
            public float3 e2_2DN;
            public float3 e1_3DN;
            public float3 e2_3DN;
        }


        // ================================ //
        // ===== BOIDS REPRESENTATION ===== //
        // ================================ //
        
        [System.Serializable]
        public struct Boid {
            public int obstacleID;
            public float3 position;
            public int3 gridIndices;
            public int projectedGridIndex;
        }


    }
}