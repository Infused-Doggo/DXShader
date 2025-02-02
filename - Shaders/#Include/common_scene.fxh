//==============================//
//           GENERAL : 
//==============================//

static float4 g_material = float4(1, 1, 1, saturate(MaterialDiffuse.w) * 1 );
float2 gl_FragCoord;

static float4 g_npr_cloth_spec_color = float4(LightSpecular.xyz * 2, 0.20);
#ifdef PMX_Color
float4 PMX_Specular = Specular;
	static float4 g_material_state_specular = float4(MaterialSpecular.xyz, length(LightAmbient.xyz)* Specular.w);
#else
	static float4 g_material_state_specular = float4(Specular.xyz, length(LightAmbient.xyz) * Specular.w);
#endif

//////////////////////////////////////////////////////////////////////////////////////////

#ifdef obj_ID
static float4 g_texcoord_transforms[2] = {
#ifdef Has_TX1
float4((TX1_RPT + 1 == 1 ? 1 : TX1_RPT).x, (TX1_RPT + 1 == 1 ? 1 : TX1_RPT).y, -TX1_TRF.xy),
float4((TX0_RPT + 1 == 1 ? 1 : TX0_RPT).x, (TX0_RPT + 1 == 1 ? 1 : TX0_RPT).y, -TX0_TRF.xy)};
#else
float4((TX0_RPT + 1 == 1 ? 1 : TX0_RPT).x, (TX0_RPT + 1 == 1 ? 1 : TX0_RPT).y, -TX0_TRF.xy),
float4((TX1_RPT + 1 == 1 ? 1 : TX1_RPT).x, (TX1_RPT + 1 == 1 ? 1 : TX1_RPT).y, -TX1_TRF.xy)};
#endif
#else
float4 g_texcoord_transforms[2] = {
float4(1.00, 1.00, 1.00, 1.00),
float4(1.00, 1.00, 1.00, 1.00)};
#endif
	

static float4 g_light_env_chara_diffuse = float4(Override ? 1 : LightAmbient.xyz * 1.6, 1.0);



static float4 g_view_position = float4(CameraPosition.xyz, 1.0);

//#Batch
static float4 g_material_state_diffuse = PMX_Color ? MaterialDiffuse : Diffuse;
static float4 g_material_state_ambient = PMX_Color ? MaterialAmbient : Ambient;
static float4 g_material_state_emission = PMX_Color ? Emission : Emission;
