/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 */

export const ssao_frag = `
precision highp float;
precision highp int;
precision highp sampler2D;

#include common

uniform sampler2D tDepth;
uniform vec2 uTexSize;
uniform vec4 uBounds;

uniform float uNear;
uniform float uFar;

#if dLightCount != 0
    uniform vec3 uLightDirection[dLightCount];
    uniform vec3 uLightColor[dLightCount];
#endif

uniform vec3 uSamples[dNSamples];

uniform mat4 uProjection;
uniform mat4 uView;
uniform mat4 uInvProjection;

uniform float uRadius;
uniform float uBias;

// shadow uniform
uniform float uSDistance;
uniform float uSTolerance;
uniform float uSBias;
uniform int uShadow;

//ssao-pro uniform
uniform int uCloseAO;
uniform float uCloseBias;
uniform float uCloseDistance;
uniform float uCDistanceCutoff;
uniform float uCCutoffFalloff;
uniform float uCIntensity;
uniform float uCDistance;
    
//ssao-old-blender uniform
uniform int uSoftAO;
uniform float uAorange;
uniform float uDepthTolerance;
uniform float uAoMultiplier;
uniform float uAoCap;
uniform float uAScale;
uniform int uARings;
uniform int uASamples;

#define PI 3.14159265
#define SAMPLES_HIGH 1
#define SAMPLES_ULTRA 0
#define SAMPLE_NOISE 1


float smootherstep(float edge0, float edge1, float x) {
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

float noise(const in vec2 coords) {
    float a = 12.9898;
    float b = 78.233;
    float c = 43758.5453;
    float dt = dot(coords, vec2(a,b));
    float sn = mod(dt, PI);
    return abs(fract(sin(sn) * c)); // is abs necessary?
}

vec2 getNoiseVec2(const in vec2 coords) {
    return vec2(noise(coords), noise(coords + vec2(PI, 2.71828)));
}

bool isBackground(const in float depth) {
    return depth == 1.0;
}

bool outsideBounds(const in vec2 p) {
    return p.x < uBounds.x || p.y < uBounds.y || p.x > uBounds.z || p.y > uBounds.w;
}

float getViewZ(in float depth) {
    #if dOrthographic == 1
        return orthographicDepthToViewZ(depth, uNear, uFar);
    #else
        return perspectiveDepthToViewZ(depth, uNear, uFar);
    #endif
}

float getDepth(const in vec2 coords) {
    if (outsideBounds(coords)) {
        return 1.0;
    } else {
        #ifdef depthTextureSupport
            return texture2D(tDepth, coords).r;
        #else
            return unpackRGBAToDepth(texture2D(tDepth, coords));
        #endif
    }
}

vec3 normalFromDepth(const in float depth, const in float depth1, const in float depth2, vec2 offset1, vec2 offset2) {
    vec3 p1 = vec3(offset1, depth1 - depth);
    vec3 p2 = vec3(offset2, depth2 - depth);

    vec3 normal = cross(p1, p2);
    normal.z = -normal.z;

    return normalize(normal);
}

float readDepth( in vec2 coord ) {
	return (2.0 * uNear) / (uFar + uNear - getDepth(coord) * (uFar - uNear));	
}

float compareDepths( in float d1, in float d2 ){
    float near = uNear;
    float far = uFar;
    //float aorange = 160.0; //uniform
    //float depthTolerance = 0.0;//uniform
    //float aoMultiplier = 100.0;//uniform
    //float aoCap = 1.0;//uniform
    //go linear ?
    float depth1 = d1;//getViewZ(d1);
    float depth2 = d2;//getViewZ(d2);
    //float diff = sqrt(clamp(1.0-(depth1-depth2) / (uAorange),0.0,1.0));
    float diff = sqrt(clamp(1.0-(depth1-depth2) / (uAorange/(far-near)),0.0,1.0));
    float ao = min(uAoCap,max(0.0,depth1-depth2-uDepthTolerance) * uAoMultiplier) * diff;
    return ao;
}

float computeAO(in vec2 scrCoord){
    float depth = readDepth(scrCoord);
    vec2 invTexSize = 1.0 / uTexSize;
    int do_noise = 0;
    
    float scale = uAScale; //uniform
    float aspect = uTexSize.x/uTexSize.y;
    int rings = uARings;//min(6,int(uRadius)); //uniform
    int samples = uASamples;//min(6,int(dNSamples)); //uniform
    //vec3 randomVec = normalize(vec3(getNoiseVec2(scrCoord), 0.0));
    vec2 noise = getNoiseVec2(scrCoord);//getRandom(srcCoord);//
    float w;
	float h;
	if (do_noise == 1) {
       w = invTexSize.x/clamp(depth,0.05,1.0)+(noise.x*(1.0-noise.x))*scale;
       h = invTexSize.y/clamp(depth,0.05,1.0)+(noise.y*(1.0-noise.y))*scale;
	   }
    else {
       w = invTexSize.x/clamp(depth,0.05,1.0)+0.001*scale;//+(noise.x*(1.0-noise.x));
       h = invTexSize.y/clamp(depth,0.05,1.0)+0.001*scale;//+(noise.y*(1.0-noise.y));
	}
    float pw;
    float ph;

    float ao;
    float s;

    int ringsamples;
    for (int i = 1; i <= rings; i += 1){
        ringsamples = i * samples;
        for (int j = 0 ; j < ringsamples ; j += 1)   {
            float step = PI*2.0 / float(ringsamples);
            pw = (cos(float(j)*step)*float(i));
            ph = (sin(float(j)*step)*float(i))*aspect;
            float v = readDepth( vec2(scrCoord.s+pw*w,scrCoord.t+ph*h) );
            ao += compareDepths(depth, v);
            s += 1.0;
            }
        }
    ao /= s;
    // ao = 1.0-ao;
    return ao;
}


float computeOcclusion(in float aradius, in mat3 TBN, in vec3 selfViewPos ){
    float occlusion = 0.0;
    for(int i = 0; i < dNSamples; i++){
        vec3 sampleViewPos = TBN * uSamples[i];
        sampleViewPos = selfViewPos + sampleViewPos * aradius;

        vec4 offset = vec4(sampleViewPos, 1.0);
        offset = uProjection * offset;
        offset.xyz = (offset.xyz / offset.w) * 0.5 + 0.5;

        float sampleViewZ = screenSpaceToViewSpace(vec3(offset.xy, getDepth(offset.xy)), uInvProjection).z;

        occlusion += step(sampleViewPos.z + 0.025, sampleViewZ) * smootherstep(0.0, 1.0, aradius / abs(selfViewPos.z - sampleViewZ));
    }
    return occlusion;
}

float calcAO(in vec2 tcoord, in vec2 uv, in vec3 p, in vec3 cnorm)
{
    float _Bias = uCloseBias;
    float _Intensity = uCIntensity;
    float _Distance = uCDistance;
    vec2 t = tcoord + uv;
    float depth = getDepth(t);
    vec3 diff = screenSpaceToViewSpace(vec3(t, depth), uInvProjection) - p;
    vec3 v = normalize(diff);
    float d = length(diff) * _Distance;
    // cnorm = normalize(_WorldSpaceCameraPos - p);
    return max(0.0, dot(cnorm, v) - _Bias) * (1.0 / (1.0 + d)) * _Intensity;
}

float invlerp(float from, float to, float value)
{
    return (value - from) / (to - from);
}

// Gold Noise function
float PHI = 1.61803398874989484820459 * 00000.1; // Golden Ratio   
float PIT  = 3.14159265358979323846264 * 00000.1; // PI
float SRT = 1.41421356237309504880169 * 10000.0; // Square Root of Two

float random_0t1(in vec2 coordinate, in float seed)
{
    return fract(sin(dot(coordinate*seed, vec2(PHI, PIT)))*SRT);
}

float ssao(in vec2 uv, in vec3 normal)
{
    float _SampleRadius = 5.0;
    float _DistanceCutoff = uCDistanceCutoff;//100.0;
    float _CutoffFalloff = uCCutoffFalloff;//25.0;

    vec2 CROSS[4] = vec2[4]( vec2(1.0, 0.0), vec2(-1.0, 0.0), vec2(0.0, 1.0), vec2(0.0, -1.0) );
    float depth = getDepth(uv);
    float eyeDepth = getViewZ(depth);
    vec3 position = screenSpaceToViewSpace(vec3(uv, depth), uInvProjection);
    float radius =  uCloseDistance; // original was max(_SampleRadius / eyeDepth, 0.005);
    // clip(_DistanceCutoff - eyeDepth); // Skip out of range pixels
    if (_DistanceCutoff - abs(eyeDepth) < 0.0) return 1.0;
    #if defined(SAMPLE_NOISE)
        float a = random_0t1(uv,depth);
        float b = random_0t1(uv,eyeDepth);
        vec2 random = normalize(vec2(a,b));
        // original used a texture for noise
        // vec2 random = normalize(tex2D(_NoiseTex, _ScreenParams.xy * uv / _NoiseSize).rg * 2.0 - 1.0);
    #endif    
    float ao = 0.0;
    // Sampling
    for (int j = 0; j < 4; j++)
    {
        vec2 coord1;

        #if defined(SAMPLE_NOISE)
            coord1 = reflect(CROSS[j], random) * radius;
        #else
            coord1 = CROSS[j] * radius;
        #endif

        // #if !SAMPLES_VERY_LOW
        vec2 coord2 = coord1 * 0.707;
        coord2 = vec2(coord2.x - coord2.y, coord2.x + coord2.y);
        // #endif 

        #if defined(SAMPLES_ULTRA)    // 20
        ao += calcAO(uv, coord1 * 0.20, position, normal);
        ao += calcAO(uv, coord2 * 0.40, position, normal);
        ao += calcAO(uv, coord1 * 0.60, position, normal);
        ao += calcAO(uv, coord2 * 0.80, position, normal);
        ao += calcAO(uv, coord1, position, normal);
        #elif defined(SAMPLES_HIGH)          // 16
        ao += calcAO(uv, coord1 * 0.25, position, normal);
        ao += calcAO(uv, coord2 * 0.50, position, normal);
        ao += calcAO(uv, coord1 * 0.75, position, normal);
        ao += calcAO(uv, coord2, position, normal);
        #elif defined(SAMPLES_MEDIUM)        // 12
        ao += calcAO(uv, coord1 * 0.30, position, normal);
        ao += calcAO(uv, coord2 * 0.60, position, normal);
        ao += calcAO(uv, coord1 * 0.90, position, normal);
        #elif defined(SAMPLES_LOW )          // 8
        ao += calcAO(uv, coord1 * 0.30, position, normal);
        ao += calcAO(uv, coord2 * 0.80, position, normal);
        #else   // 4
        ao += calcAO(uv, coord1 * 0.50, position, normal);
        #endif
    }
    
    #if SAMPLES_ULTRA
    ao /= 20.0;
    #elif SAMPLES_HIGH
    ao /= 16.0;
    #elif SAMPLES_MEDIUM
    ao /= 12.0;
    #elif SAMPLES_LOW
    ao /= 8.0;
    #else
    ao /= 4.0;
    #endif

    // Distance cutoff
    ao = mix(1.0 - ao, 1.0, saturate(invlerp(_DistanceCutoff - _CutoffFalloff, _DistanceCutoff, eyeDepth)));

    return ao;
}

float ScreenSpaceShadows(in vec2 uv, in vec3 position, in vec3 light_direction)
{
    // Settings
    int  g_sss_steps            = dSSample; // Quality/performancedNSamples
    float g_sss_ray_max_distance = uSDistance; // Max shadow length
    float g_sss_tolerance        = uSTolerance; // Error in favor of reducing gaps
    float g_sss_step_length      = g_sss_ray_max_distance / float(g_sss_steps);
    float uvdepth = getDepth(uv);
    
    float eyeDepth = getViewZ(uvdepth);
    // Compute ray position and direction (in view-space)
    vec3 ray_pos = position;
    vec3 ray_dir = -light_direction; //light direction in View space
	vec2 uv_pos = uv;
    // Compute ray step
    vec3 ray_step = ray_dir * g_sss_step_length;
	vec2 uv_step = ray_dir.xy * g_sss_step_length;
    // Ray march towards the light
    float occlusion = 0.0;
    
    vec4 ray_uv   = vec4(0.0,0.0,0.0,0.0);
    for (int i = 0; i < g_sss_steps; i++)
    {
        // Step the ray
        uv_pos += uv_step;
        ray_pos += ray_step;
        // Compute the difference between the ray's and the camera's depth
        ray_uv  = (uProjection * vec4(ray_pos,1.0));
        ray_uv.xyz = (ray_uv.xyz / ray_uv.w) * 0.5 + 0.5;
        float depth = getDepth(ray_uv.xy);       
        float depth_z = getViewZ(depth);
        float depth_delta = ray_pos.z - depth_z;
        if (depth_delta < g_sss_tolerance){
        // original test : if (abs(g_sss_tolerance - depth_delta) < g_sss_tolerance){
            occlusion = 1.0;
            vec2 fade = max(12.0 * abs(ray_uv.xy - 0.5) - 5.0, vec2(0.0,0.0));
            occlusion *= saturate(1.0 - dot(fade, fade));             
            break;
        }
    }
    // Fade out as we approach the edges of the screen
    // occlusion *= screen_fade(ray_uv);return 1.0 - occlusion;
    occlusion = 1.0 - (uSBias * occlusion);
    return occlusion;
}

// StarCraft II Ambient Occlusion by [Filion and McNaughton 2008]
void main(void) {
    vec2 invTexSize = 1.0 / uTexSize;
    vec2 selfCoords = gl_FragCoord.xy * invTexSize;

    float selfDepth = getDepth(selfCoords);
    vec2 selfPackedDepth = packUnitIntervalToRG(selfDepth);

    if (isBackground(selfDepth)) {
        gl_FragColor = vec4(packUnitIntervalToRG(0.0), selfPackedDepth);
        return;
    }

    vec2 offset1 = vec2(0.0, invTexSize.y);
    vec2 offset2 = vec2(invTexSize.x, 0.0);

    float selfDepth1 = getDepth(selfCoords + offset1);
    float selfDepth2 = getDepth(selfCoords + offset2);

    vec3 selfViewNormal = normalFromDepth(selfDepth, selfDepth1, selfDepth2, offset1, offset2);
    vec3 selfViewPos = screenSpaceToViewSpace(vec3(selfCoords, selfDepth), uInvProjection);

    vec3 randomVec = normalize(vec3(getNoiseVec2(selfCoords) * 2.0 - 1.0, 0.0));

    vec3 tangent = normalize(randomVec - selfViewNormal * dot(randomVec, selfViewNormal));
    vec3 bitangent = cross(selfViewNormal, tangent);
    mat3 TBN = mat3(tangent, bitangent, selfViewNormal);

    float occlusion = computeOcclusion(uRadius, TBN, selfViewPos);
    occlusion = 1.0 - (uBias * occlusion / float(dNSamples));
    float ao1=0.0;
    // alternative ao algo
    if (uSoftAO == 1)
    {
        ao1 = computeAO(selfCoords);
        ao1 = clamp(ao1,0.0,1.0);
        if ( ao1 > 1.0 ) {ao1 = 1.0 ;}
        if ( ao1 < 0.0 ) {ao1 = 0.0 ;}
        if (selfDepth > 1.0 ) {ao1 = 1.0 ;}
        if (selfDepth < 0.0 ) {ao1 = 0.0 ;}
        ao1 = 1.0 - (ao1);
    }
    
    bool isClose = true;
    if (abs(selfViewPos.z) > 1200.0) isClose = false;

    float ao = 1.0;
    if (uCloseAO == 1){
        ao = saturate(ssao(selfCoords, selfViewNormal));
    }
    float o = 9999.9;
    if (uShadow == 1) {
        #if dLightCount != 0
        float sh[dLightCount];
        #pragma unroll_loop_start
        for (int i = 0; i < dLightCount; ++i) {
            sh[i] = ScreenSpaceShadows(selfCoords, selfViewPos, uLightDirection[i]);
            o = min(o,min(min(sh[i],ao),occlusion));
        }
        #pragma unroll_loop_end
        #endif
    }
    else{ 
        o = min(ao,occlusion);
    }
    if (uSoftAO==1){
        o = min(ao1,o);
    }
    vec2 packedOcclusion = packUnitIntervalToRG(o);
    gl_FragColor = vec4(packedOcclusion, selfPackedDepth);
}
`;