Shader "SPH/GridCell"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "Queue" = "Transparent" "RenderType" = "Transparent"}
        LOD 200

        CGPROGRAM
		#include "../../Noise/noise.cginc"
		#pragma surface surf Standard addshadow fullforwardshadows alpha:fade
		#pragma multi_compile_instancing
		#pragma instancing_options procedural:setup

        sampler2D _MainTex;
		float size;

        struct Input {
			float2 uv_MainTex;
		};

        struct GridCell {
            int n;
            float3 worldPosition;
            int cellPressure;
        };

        #ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
			StructuredBuffer<GridCell> grid_cell_buffer;
            StructuredBuffer<float> pressures_buffer;
            StructuredBuffer<float> render_limits_buffer;
		#endif

        void setup() {
            #ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
                float3 pos = grid_cell_buffer[unity_InstanceID].worldPosition;
                unity_ObjectToWorld._11_21_31_41 = float4(size, 0, 0, 0);
                unity_ObjectToWorld._12_22_32_42 = float4(0, size, 0, 0);
                unity_ObjectToWorld._13_23_33_43 = float4(0, 0, size, 0);
                unity_ObjectToWorld._14_24_34_44 = float4(pos.xyz, 1);
                unity_WorldToObject = unity_ObjectToWorld;
                unity_WorldToObject._14_24_34 *= -1;
                unity_WorldToObject._11_22_33 = 1.0f / unity_WorldToObject._11_22_33;
            #endif
        }

	    void surf(Input IN, inout SurfaceOutputStandard o) {
		    // blue
		    float4 lowPressureColor = float4(0.25, 0.5, 1.0, 1.0);
		    // yellow
		    float4 highPressureColor = float4(1.0, 0.0, 0.0, 1.0);
		    float4 finalColor = lowPressureColor;

    		#ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
                float p = max(0.0, min(1.0, pressures_buffer[unity_InstanceID]/(0.080214*150)));
                //Gizmos.color = Color.red * p + Color.blue * (1f-p);
			    //float a = clamp(length(grid_cell_buffer[unity_InstanceID].cellPressure)/25.0, 0.0, 1.0);
			    finalColor = lowPressureColor * (1.0-p) + highPressureColor * p; 

                // Also do not render if we have a limit set
                float3 pos = grid_cell_buffer[unity_InstanceID].worldPosition;
                if (
                    pos[0] < render_limits_buffer[0] || pos[1] < render_limits_buffer[1] || pos[2] < render_limits_buffer[2] 
                    || pos[0] > render_limits_buffer[3] || pos[1] > render_limits_buffer[4] || pos[2] > render_limits_buffer[5]
                ) finalColor.a = 0.0;

		    #endif

		    float4 c = finalColor;
		    /* #ifdef UNITY_PROCEDURAL_INSTANCING_ENABLED
		    c = particle_buffer[unity_InstanceID].color;
		    #endif */
		    c *= tex2D(_MainTex, IN.uv_MainTex);
		    o.Albedo = c.rgb + noise(c.rgb);
		    o.Alpha = c.a;
	    }

    ENDCG
    }

    FallBack "Diffuse"

}
