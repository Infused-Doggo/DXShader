
float R_OffsetA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="R_Offset +";>;
float R_OffsetB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="R_Offset -";>;
float G_OffsetA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="G_Offset +";>;
float G_OffsetB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="G_Offset -";>;
float B_OffsetA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="B_Offset +";>;
float B_OffsetB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="B_Offset -";>;
float R_ScaleA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="R_Scale +";>;
float R_ScaleB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="R_Scale -";>;
float G_ScaleA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="G_Scale +";>;
float G_ScaleB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="G_Scale -";>;
float B_ScaleA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="B_Scale +";>;
float B_ScaleB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="B_Scale -";>;

float SaturationA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Saturation +";>;
float SaturationB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Saturation -";>;
float ExposureA   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Exposure +";>;
float ExposureB   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Exposure -";>;
float GammaA      : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma +";>;
float GammaB      : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma -";>;

float tone_type : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Tonemap_Type";>;
float auto_exposure  : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Auto_Exposure";>;
float Saturation_Pow : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Saturation_Pow";>;
float Fade_alpha     : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Fade_Alpha";>;
float Fade_A         : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Fade_R +";>;
float Fade_B         : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Fade_G +";>;
float Fade_C         : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Fade_B +";>;
float Override_TM    : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Override";>;

float set(float A, float B) {
    return lerp(1 + (A * 1.5) * 1 - B, A, (int)Override_TM);
}

float set2(float A, float B) {
    return lerp((A * 1.5) - (B * 1.5), A, (int)Override_TM);
}

#define EPS 0.001
#define PI  3.1415926535

static float exposure = set(ExposureA, ExposureB);
static float exposure_rate = 1.5;
static float4 g_exposure    = float4(exposure * exposure_rate, 0.0625f, exposure * exposure_rate * 0.5f, auto_exposure ? 1.0f : 0.0f);
static float4 g_fade_color  = float4(Fade_A, Fade_B, Fade_C, Fade_alpha);
static float4 g_tone_scale  = float4(set(R_ScaleA, R_ScaleB), set(G_ScaleA, G_ScaleB), set(B_ScaleA, B_ScaleB), 0.00);
static float4 g_tone_offset = float4(set2(R_OffsetA, R_OffsetB), set2(G_OffsetA, G_OffsetB), set2(B_OffsetA, B_OffsetB), 0.66667);

texture2D RampTex;
sampler2D g_ramp_s = sampler_state {
	texture = <RampTex>;
    FILTER = ANISOTROPIC;
    ADDRESSU = CLAMP;
    ADDRESSV = CLAMP;
};

shared texture2D g_tonemap : RENDERCOLORTARGET;
sampler2D g_tone_map_s = sampler_state {
    texture = <g_tonemap>;
    MinFilter = LINEAR;
    MagFilter = LINEAR;
    MipFilter = LINEAR;
    AddressU  = WRAP;
    AddressV  = WRAP;
};

bool TM_x : CONTROLOBJECT < string name = "ToneMap.x"; >;
bool EX_x : CONTROLOBJECT < string name = "Bloom.x"; >;
const float4 to_ybr = float4(0.3, 0.59, 0.11, 1.0);
const float4 to_rgb = float4(-0.508475, 1.0, -0.186441, 1.0);

float3 apply_tonemap(float3 color) {

	g_tone_scale.xyz = g_tone_scale.xyz * 1.1;
	g_tone_offset.xyz = g_tone_offset.xyz * lerp(1.1, -1.1, (int)Override_TM);
	float2 frg_exposure = float2(g_exposure.x, g_exposure.y * g_exposure.x);
	
	float4 col = 0;
    float3 sum = color;  
	
	float3 res;
	if (tone_type >= 0.66) { // #if TONE_MAP_2_DEF
        res = min(sum.rgb * 0.25 * frg_exposure.x, 0.80);
    }
    else if (tone_type >= 0.33) { // #elif TONE_MAP_1_DEF
        res = min(sum.rgb * 0.48 * frg_exposure.x, 0.96);
    }
    else { // #else
        float3 ybr;
        ybr.y = dot(sum.rgb, to_ybr.xyz);
        ybr.xz = sum.xz - ybr.y;

			col = tex2D(g_tone_map_s, float2(ybr.y * frg_exposure.y, 0.0)).xxxy;

        col.xz = col.w * frg_exposure.x * ybr.xz;
        res.rb = col.xz + col.y;
        res.g = dot(col.rgb, to_rgb.xyz);
    } // #endif
	
    res = clamp(res * g_tone_scale.rgb + g_tone_offset.rgb, (0.0), (1.0));

        const float blend = g_tone_scale.w;
        bool3 cc = (blend) == float3(0.0, 1.0, 2.0);
        res = lerp(res, lerp(res, g_fade_color.rgb, g_fade_color.a), float(cc.x));
        res = lerp(res, res * g_fade_color.rgb, float(cc.y));
        res = lerp(res, res + g_fade_color.rgb, float(cc.z));
  
	if (TM_x || EX_x) {
		res = color;
	}
	return res;
}

float3 apply_fog_color(float3 color, float4 fog_color) {
    return apply_tonemap(lerp(color, fog_color.rgb, fog_color.w));
}

float3 get_ibl_diffuse(samplerCUBE tex, float3 ray, float lc) {
    float3 col0 = texCUBElod(tex, float4(ray, 0.0)).rgb;
    float3 col1 = texCUBElod(tex, float4(ray, 1.0)).rgb;
    return lerp(col1, col0, lc);
}

// - - - ° Addition ° - - - //

float3 inv(float3 x) {
    return x * float3(1, 1, -1);
}

float4x4 CTF(float3 frg_position, float4 frg_normal, float4 frg_texcoord) {
	float4x4 Out;
	float3 p_dx = ddx(frg_position.xyz);
	float3 p_dy = ddy(frg_position.xyz);
	float2 tc_dx = ddx(frg_texcoord.xy);
	float2 tc_dy = ddy(frg_texcoord.xy);
	float direction = tc_dx.x * tc_dy.y - tc_dx.y * tc_dy.x > 0.0f ? 1.0f : -1.0f;
	float3 t = normalize(tc_dy.y * p_dx - tc_dx.y * p_dy) * direction;
	float3 b = normalize( (tc_dy.x * p_dx - tc_dx.x * p_dy) * direction );
	float3 n = normalize(frg_normal.xyz);
	float3 x = cross(n, t);
	t = cross(x, n);
	t = normalize(t);
	x = cross(b, n);
	b = cross(n, x);
	b = normalize(b);
	
	Out[0].xyz = t;
	Out[1].xyz = b;
	Out[2] = frg_normal;
	return Out;
}

float3 Light_Position(float3 pos_dir)
{
		pos_dir *= float3(1, 1, -1);
		float flength = length(pos_dir);
        if (flength <= 0.000001)
            pos_dir = float3(0.0f, 1.0f, 0.0f);
        else
            pos_dir *= 1.0f / flength;
        return pos_dir;
}


float DistributionGGX(float NdotH, float roughness){
    float a = roughness * roughness;
	float NdotH2 = NdotH * NdotH;
    float denom = NdotH2 * (a - 1.0) + 1.0;
    return 1 / (PI * denom * denom);
}

float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float denom = NdotV * (1.0 - k) + k;
    return (1.0 / denom) * NdotV;
}

float GeometrySmith(float NdotL, float NdotV, float roughness)
{
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    return ggx1 * ggx2;
}

float FresnelSchlick(float F0, float cosTheta)
{
    return min(1.0, F0 + (1.0 - F0) * pow(1.0 - saturate(cosTheta), 5.0));
}

float3 PBR(float3 N, float3 L, float3 V, float roughness, float metalness)
{
	L = normalize( L );
	V = normalize( V );
	float3 H = normalize(L + V);
	
    float NdotL = max(dot(N, L), 0.0);
    float NdotV = max(dot(N, V), 0.0);
    float NdotH = max(dot(N, H), 0.0);
	
	float Exp = EPS + roughness;
	
	float NDF = DistributionGGX(NdotH, Exp);
	float G = GeometrySmith(NdotL, NdotV, roughness);
	float F = FresnelSchlick(metalness, dot(H, V));
	
	float nominator    = NDF * G;
	float denominator  = 1.0 / (4.0 * NdotV * NdotL);
	float BRDF = nominator * denominator * (Exp * Exp);
	return float3(BRDF, NdotL, F);
}

float3 sRGB2Lin(float3 color)
{
    return exp2(2.2 * log2(abs(color)));
	//return color;
}

float4 IBLTransform(float3 target, float3 viewVec, float3 refl, float idkkboth)
{
    float fixed = dot(target, viewVec);
	float dir = dot(target, target);
	
	float4 output;
	float _1113 = idkkboth - dir;
	float _1222 = sqrt(fixed * fixed + _1113) + fixed;
	float3 rotateVec = normalize(_1222 * refl + target);
	rotateVec.xy = (rotateVec.xy * (1 / abs(1)) + 1.5) - 1.5;
	
	output.xyz = rotateVec;
	output.w = sqrt(idkkboth) + (1 / -rsqrt(dir));
	return output;
}