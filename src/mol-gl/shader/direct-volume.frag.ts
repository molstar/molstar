/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export default `
precision highp float;
precision highp int;

#include common
#include light_frag_params

#include read_from_texture
#include texture3d_from_2d_nearest
#include texture3d_from_2d_linear

uniform mat4 uProjection, uTransform, uModelView, uView;
uniform vec3 uCameraDir;

uniform sampler2D tDepth;
uniform vec4 uViewport;
uniform float uNear;
uniform float uFar;

varying vec3 unitCoord;
varying vec3 origPos;
varying float instance;

uniform mat4 uInvView;
uniform vec2 uIsoValue;
uniform vec3 uGridDim;
uniform sampler2D tTransferTex;
uniform float uStepFactor;

uniform int uObjectId;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
uniform vec2 uMarkerTexDim;
uniform sampler2D tMarker;

uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform int uTransparentBackground;

uniform float uInteriorDarkening;
uniform int uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;

uniform float uIsOrtho;

uniform vec3 uCellDim;
uniform vec3 uCameraPosition;
uniform mat4 uCartnToUnit;
uniform mat4 uUnitToCartn;

#if __VERSION__ == 300
    // for webgl1 this is given as a 'define'
    uniform int uMaxSteps;
#endif

#if defined(dGridTexType_2d)
    precision highp sampler2D;
    uniform sampler2D tGridTex;
    uniform vec3 uGridTexDim;
#elif defined(dGridTexType_3d)
    precision highp sampler3D;
    uniform sampler3D tGridTex;
#endif

#if defined(dRenderVariant_color)
    #if defined(dColorType_uniform)
        uniform vec3 uColor;
    #elif defined(dColorType_texture)
        uniform vec2 uColorTexDim;
        uniform sampler2D tColor;
    #endif

    #ifdef dOverpaint
        varying vec4 vOverpaint;
        uniform vec2 uOverpaintTexDim;
        uniform sampler2D tOverpaint;
    #endif
#endif

#if defined(dGridTexType_2d)
    vec4 textureVal(vec3 pos) {
        return texture3dFrom2dLinear(tGridTex, pos + (vec3(0.5, 0.5, 0.0) / uGridDim), uGridDim, uGridTexDim.xy);
    }
    vec4 textureGroup(vec3 pos) {
        return texture3dFrom2dNearest(tGridTex, pos + (vec3(0.5, 0.5, 0.0) / uGridDim), uGridDim, uGridTexDim.xy);
    }
#elif defined(dGridTexType_3d)
    vec4 textureVal(vec3 pos) {
        return texture(tGridTex, pos + (vec3(0.5) / uGridDim));
    }
    vec4 textureGroup(vec3 pos) {
        return texelFetch(tGridTex, ivec3(pos * uGridDim), 0);
    }
#endif

vec4 transferFunction(float value) {
    return texture2D(tTransferTex, vec2(value, 0.0));
}

// Calculate depth based on the given camera position.
float calcDepth(const in vec3 cameraPos){
    vec2 clipZW = cameraPos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}

float getDepth(const in vec2 coords) {
    #ifdef dPackedDepth
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #else
        return texture2D(tDepth, coords).r;
    #endif
}

const float gradOffset = 0.5;

vec3 toUnit(vec3 p) {
    return (uCartnToUnit * vec4(p, 1.0)).xyz;
}

vec4 raymarch(vec3 startLoc, vec3 step) {
    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec3 pos = startLoc;
    vec4 cell;
    vec4 prevCell = vec4(-1);
    float prevValue = -1.0;
    float value = 0.0;
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);
    bool hit = false;

    vec3 posMin = vec3(0.0);
    vec3 posMax = vec3(1.0) - vec3(1.0) / uGridDim;

    vec3 unitPos;
    vec3 isoPos;

    vec3 color = vec3(0.45, 0.55, 0.8);
    vec3 gradient = vec3(1.0);
    vec3 dx = vec3(gradOffset * scaleVol.x, 0.0, 0.0);
    vec3 dy = vec3(0.0, gradOffset * scaleVol.y, 0.0);
    vec3 dz = vec3(0.0, 0.0, gradOffset * scaleVol.z);

    for(int i = 0; i < uMaxSteps; ++i){
        unitPos = toUnit(pos);
        if(unitPos.x > posMax.x || unitPos.y > posMax.y || unitPos.z > posMax.z || unitPos.x < posMin.x || unitPos.y < posMin.y || unitPos.z < posMin.z) {
            if (hit) break;
            prevValue = value;
            prevCell = cell;
            pos += step;
            continue;
        }

        cell = textureVal(unitPos);
        value = cell.a; // current voxel value

        #if defined(dRenderMode_isosurface)
            if(prevValue > 0.0 && ( // there was a prev Value
                (prevValue < uIsoValue.x && value > uIsoValue.x) || // entering isosurface
                (prevValue > uIsoValue.x && value < uIsoValue.x) // leaving isosurface
            )) {
                isoPos = toUnit(mix(pos - step, pos, ((prevValue - uIsoValue.x) / ((prevValue - uIsoValue.x) - (value - uIsoValue.x)))));

                vec4 mvPosition = uModelView * uTransform * vec4(isoPos * uGridDim, 1.0);
                float depth = calcDepth(mvPosition.xyz);
                if (depth > getDepth(gl_FragCoord.xy / uViewport.zw))
                    break;

                #ifdef enabledFragDepth
                    if (!hit) {
                        gl_FragDepthEXT = depth;
                        hit = true;
                    }
                #endif

                #if defined(dRenderVariant_pickObject)
                    return vec4(encodeFloatRGB(float(uObjectId)), 1.0);
                #elif defined(dRenderVariant_pickInstance)
                    return vec4(encodeFloatRGB(instance), 1.0);
                #elif defined(dRenderVariant_pickGroup)
                    #ifdef dPackedGroup
                        return vec4(textureGroup(floor(isoPos * uGridDim + 0.5) / uGridDim).rgb, 1.0);
                    #else
                        vec3 g = floor(isoPos * uGridDim + 0.5);
                        return vec4(encodeFloatRGB(g.z + g.y * uGridDim.z + g.x * uGridDim.z * uGridDim.y), 1.0);
                    #endif
                #elif defined(dRenderVariant_depth)
                    #ifdef enabledFragDepth
                        return packDepthToRGBA(gl_FragDepthEXT);
                    #else
                        return packDepthToRGBA(gl_FragCoord.z);
                    #endif
                #elif defined(dRenderVariant_color)
                    #ifdef dPackedGroup
                        float group = decodeFloatRGB(textureGroup(floor(isoPos * uGridDim + 0.5) / uGridDim).rgb);
                    #else
                        vec3 g = floor(isoPos * uGridDim + 0.5);
                        float group = g.z + g.y * uGridDim.z + g.x * uGridDim.z * uGridDim.y;
                    #endif

                    #if defined(dColorType_instance)
                        color = readFromTexture(tColor, instance, uColorTexDim).rgb;
                    #elif defined(dColorType_group)
                        color = readFromTexture(tColor, group, uColorTexDim).rgb;
                    #elif defined(dColorType_groupInstance)
                        color = readFromTexture(tColor, instance * float(uGroupCount) + group, uColorTexDim).rgb;
                    #elif defined(dColorType_uniform)
                        color = uColor;
                    #endif

                    // handle flipping and negative isosurfaces
                    #ifdef dFlipSided
                        bool flipped = value < uIsoValue.y; // flipped
                    #else
                        bool flipped = value > uIsoValue.y;
                    #endif
                    interior = value < uIsoValue.x && flipped;
                    #ifndef dDoubleSided
                        if (interior) {
                            prevValue = value;
                            prevCell = cell;
                            pos += step;
                            continue;
                        }
                    #endif
                    vec3 vViewPosition = mvPosition.xyz;
                    vec4 material = vec4(color, uAlpha);

                    #ifdef dIgnoreLight
                        gl_FragColor = material;
                    #else
                        #if defined(dFlatShaded)
                            // nearest grid point
                            isoPos = floor(isoPos * uGridDim + 0.5) / uGridDim;
                        #endif
                        #ifdef dPackedGroup
                            // compute gradient by central differences
                            gradient.x = textureVal(isoPos - dx).a - textureVal(isoPos + dx).a;
                            gradient.y = textureVal(isoPos - dy).a - textureVal(isoPos + dy).a;
                            gradient.z = textureVal(isoPos - dz).a - textureVal(isoPos + dz).a;
                        #else
                            gradient = textureVal(isoPos).xyz * 2.0 - 1.0;
                        #endif
                        mat3 normalMatrix = transpose3(inverse3(mat3(uModelView)));
                        vec3 normal = -normalize(normalMatrix * normalize(gradient));
                        normal = normal * (float(flipped) * 2.0 - 1.0);
                        normal = normal * -(float(interior) * 2.0 - 1.0);
                        #include apply_light_color
                    #endif

                    float vMarker = readFromTexture(tMarker, instance * float(uGroupCount) + group, uMarkerTexDim).a;
                    #include apply_interior_color
                    #include apply_marker_color
                    #include apply_fog

                    src.rgb = gl_FragColor.rgb;
                    src.a =  gl_FragColor.a;

                    src.rgb *= src.a;
                    dst = (1.0 - dst.a) * src + dst; // standard blending
                #endif

                #ifdef dSingleLayer
                    break;
                #endif
            }
            prevValue = value;
        #endif

        #if defined(dRenderMode_volume)
            isoPos = toUnit(pos);
            vec4 mvPosition = uModelView * uTransform * vec4(isoPos * uGridDim, 1.0);
            if (calcDepth(mvPosition.xyz) > getDepth(gl_FragCoord.xy / uViewport.zw))
                break;

            // bool flipped = value > uIsoValue.y; // negative isosurfaces
            // interior = value < uIsoValue.x && flipped;
            vec3 vViewPosition = mvPosition.xyz;
            vec4 material = transferFunction(value);
            src.a = material.a * uAlpha;

            if (material.a >= 0.01) {
                #ifdef dPackedGroup
                    // compute gradient by central differences
                    gradient.x = textureVal(isoPos - dx).a - textureVal(isoPos + dx).a;
                    gradient.y = textureVal(isoPos - dy).a - textureVal(isoPos + dy).a;
                    gradient.z = textureVal(isoPos - dz).a - textureVal(isoPos + dz).a;
                #else
                    gradient = cell.xyz * 2.0 - 1.0;
                #endif
                mat3 normalMatrix = transpose3(inverse3(mat3(uModelView * uUnitToCartn)));
                vec3 normal = -normalize(normalMatrix * normalize(gradient));
                // normal = normal * (float(flipped) * 2.0 - 1.0);
                // normal = normal * -(float(interior) * 2.0 - 1.0);
                #include apply_light_color
                src.rgb = gl_FragColor.rgb;
            } else {
                src.rgb = material.rgb;
            }

            src.rgb *= src.a;
            dst = (1.0 - dst.a) * src + dst; // standard blending
        #endif

        // break if the color is opaque enough
        if (dst.a > 0.95)
            break;

        pos += step;
    }
    return dst;
}

// TODO calculate normalMatrix on CPU
// TODO fix near/far clipping
// TODO support clip objects
// TODO support float texture for higher precision values

void main () {
    // TODO handle on CPU in renderloop
    #if defined(dRenderVariant_pick)
        #if defined(dRenderMode_volume)
            discard;
        #elif defined(dRenderMode_isosurface)
            if (uAlpha < uPickingAlphaThreshold)
                discard; // ignore so the element below can be picked
        #endif
    #endif

    vec3 cameraPos = uInvView[3].xyz / uInvView[3].w;
    vec3 rayDir = mix(normalize(origPos - uCameraPosition), uCameraDir, uIsOrtho);;

    // TODO: set the scale as uniform?
    float stepScale = min(uCellDim.x, min(uCellDim.y, uCellDim.z)) * uStepFactor;
    vec3 step = rayDir * stepScale;

    gl_FragColor = raymarch(origPos, step);
    if (length(gl_FragColor.rgb) < 0.00001)
        discard;
}
`;