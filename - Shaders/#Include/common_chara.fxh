
float R_GammaA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_R +";>;
float R_GammaB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_R -";>;
float G_GammaA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_G +";>;
float G_GammaB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_G -";>;
float B_GammaA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_B +";>;
float B_GammaB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Gamma_B -";>;

float HueA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Hue +";>;
float HueB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Hue -";>;
float SaturationA : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Saturation +";>;
float SaturationB : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Saturation -";>;
float LightnessA  : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Lightness +";>;
float LightnessB  : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Lightness -";>;
float ExposureA   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Exposure +";>;
float ExposureB   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Exposure -";>;
float ContrastA   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Contrast +";>;
float ContrastB   : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Contrast -";>;

float Override_TM    : CONTROLOBJECT <string name="#ToneMap_Controller.pmx"; string item="Override";>;

float set(float A, float B, float C) {
    return lerp(A, (A + C) * (1 - B), 1-(int)Override_TM);
}

float set2(float A, float B) {
    return lerp(A, 1 + (A * 1) - (B * 1), 1-(int)Override_TM); }

#define EPS 0.001
#define PI  3.1415926535

bool TM_x : CONTROLOBJECT < string name = "ToneMap.x"; >;
bool EX_x : CONTROLOBJECT < string name = "Bloom.x"; >;
const float4 to_ybr = float4(0.3, 0.59, 0.11, 0.59);
const float4 to_rgb = float4(-0.508475, 1.0, -0.186441, 1.0);

float3 Apply_Tonemap(float3 color) {
	
Hue = 1 - (Hue + HueA - HueB);
Saturation = set(SaturationA, SaturationB, Saturation);
Lightness = set(LightnessA, LightnessB, Lightness);
Exposure = set(ExposureA, ExposureB, Exposure);
Contrast = set(ContrastA, ContrastB, Contrast);

static float3 Gamma_o = float3(set2(R_GammaA, R_GammaB), set2(G_GammaA, G_GammaB), set2(B_GammaA, B_GammaB))* Gamma;

	float3 Bloom_Color = float3(1, 1, 1);
		
	float3 Bloom_Bright = float3(Saturation, Lightness, Hue * Saturation)*Exposure;
	
	Gamma_o = 0.9 / Gamma_o;
	
	#define M_LOG2E 1.4426950
	

	float3 result;
	
	float _240 = 1.3804;
	precise float _241 = _240 * 0.06667;
    precise float _243 = (-_240) * 0.06667;
	
	
    float4 sum = 1;
	sum.rgb = color;
	
	float3 col;
    col = 0;
	col *= Bloom_Color;
	sum.rgb += col;
	
	sum.rgb *= 2.2; // ?? lightness, exposure??
	
	float stillidk = 0.15;
	
	col = _240 * (1.0 / (sum.rgb * (sum.rgb * stillidk + 0.5) + (0.20 * 0.3)));
	col *= (sum.rgb * stillidk + 0.05098) * sum.rgb + 0.00392;
	col.xz -= _241;
	
	sum.rgb = col;
	sum.a = _243;
	
	float3 ybr;
	ybr.y = dot(sum.rgba, to_ybr.rgba);
    ybr.xz = sum.xz - ybr.y;
	ybr.xyz *= Bloom_Bright;
	
	float3 res;
	res.xz = ybr.xz + ybr.y;
	res.y = dot(ybr.rgb, to_rgb.xyz);
	
	res = Gamma_o * log2(abs(res));
	res = exp2(res) * (-PI * 2) * Contrast;
	res *= M_LOG2E;
	
	float idk = 12.3451;

	res = 1.0 / (exp2(res) * idk + 1);
	
	result.rgb = 1.09000003337860107421875 * res + -0.045000016689300537109375;
	return result;
}

float3 apply_fog_color(float3 color, float4 fog_color) {
    return Apply_Tonemap(lerp(color, fog_color.rgb, fog_color.w));
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