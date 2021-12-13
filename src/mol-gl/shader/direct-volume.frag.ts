/**
 * Copyright (c) 2017-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

export const directVolume_frag = `
precision highp float;
precision highp int;

#include common
#include light_frag_params

#if dClipObjectCount != 0
    uniform int uClipObjectType[dClipObjectCount];
    uniform bool uClipObjectInvert[dClipObjectCount];
    uniform vec3 uClipObjectPosition[dClipObjectCount];
    uniform vec4 uClipObjectRotation[dClipObjectCount];
    uniform vec3 uClipObjectScale[dClipObjectCount];
#endif
#include common_clip

#include read_from_texture
#include texture3d_from_1d_trilinear
#include texture3d_from_2d_nearest
#include texture3d_from_2d_linear

uniform mat4 uProjection, uTransform, uModelView, uModel, uView;
uniform vec3 uCameraDir;

uniform sampler2D tDepth;
uniform vec2 uDrawingBufferSize;

varying vec3 vOrigPos;
varying float vInstance;
varying vec4 vBoundingSphere;
varying mat4 vTransform;

uniform mat4 uInvView;
uniform vec2 uIsoValue;
uniform vec3 uGridDim;
uniform vec3 uBboxSize;
uniform sampler2D tTransferTex;
uniform float uTransferScale;
uniform float uStepScale;
uniform float uJumpLength;

uniform int uObjectId;
uniform int uVertexCount;
uniform int uInstanceCount;
uniform int uGroupCount;

uniform vec3 uHighlightColor;
uniform vec3 uSelectColor;
uniform float uHighlightStrength;
uniform float uSelectStrength;
uniform int uMarkerPriority;

#if defined(dMarkerType_uniform)
    uniform float uMarker;
#elif defined(dMarkerType_groupInstance)
    uniform vec2 uMarkerTexDim;
    uniform sampler2D tMarker;
#endif

uniform float uMetalness;
uniform float uRoughness;

uniform float uFogNear;
uniform float uFogFar;
uniform vec3 uFogColor;

uniform float uAlpha;
uniform float uPickingAlphaThreshold;
uniform bool uTransparentBackground;

uniform float uInteriorDarkening;
uniform bool uInteriorColorFlag;
uniform vec3 uInteriorColor;
bool interior;

uniform bool uRenderWboit;

uniform float uNear;
uniform float uFar;
uniform float uIsOrtho;

uniform vec3 uCellDim;
uniform vec3 uCameraPosition;
uniform mat4 uCartnToUnit;

#if __VERSION__ != 100
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
        #if defined(dOverpaintType_groupInstance) || defined(dOverpaintType_vertexInstance)
            uniform vec2 uOverpaintTexDim;
            uniform sampler2D tOverpaint;
        #endif
    #endif

    #ifdef dSubstance
        #if defined(dSubstanceType_groupInstance) || defined(dSubstanceType_vertexInstance)
            uniform vec2 uSubstanceTexDim;
            uniform sampler2D tSubstance;
        #endif
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

float calcDepth(const in vec3 pos) {
    vec2 clipZW = pos.z * uProjection[2].zw + uProjection[3].zw;
    return 0.5 + 0.5 * clipZW.x / clipZW.y;
}

vec4 transferFunction(float value) {
    return texture2D(tTransferTex, vec2(value, 0.0));
}

float getDepth(const in vec2 coords) {
    #ifdef depthTextureSupport
        if (!uRenderWboit) {
            // in case of opaque volumes (and depth texture support)
            return texture2D(tDepth, coords).r;
        } else {
            return unpackRGBAToDepth(texture2D(tDepth, coords));
        }
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

const float gradOffset = 0.5;

vec3 v3m4(vec3 p, mat4 m) {
    return (m * vec4(p, 1.0)).xyz;
}

float preFogAlphaBlended = 0.0;

vec4 raymarch(vec3 startLoc, vec3 step, vec3 rayDir) {
    #if defined(dRenderVariant_color) && !defined(dIgnoreLight)
        mat3 normalMatrix = transpose3(inverse3(mat3(uModelView * vTransform)));
    #endif
    mat4 cartnToUnit = uCartnToUnit * inverse4(vTransform);
    #if defined(dClipVariant_pixel) && dClipObjectCount != 0
        mat4 modelTransform = uModel * vTransform * uTransform;
    #endif
    mat4 modelViewTransform = uModelView * vTransform * uTransform;

    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec3 pos = startLoc;
    vec4 cell;
    float prevValue = -1.0;
    float value = 0.0;
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);
    bool hit = false;
    float fragmentDepth;

    vec3 posMin = vec3(0.0);
    vec3 posMax = vec3(1.0) - vec3(1.0) / uGridDim;

    vec3 unitPos;
    vec3 isoPos;

    vec3 nextPos;
    float nextValue;

    vec3 color = vec3(0.45, 0.55, 0.8);
    vec4 overpaint = vec4(0.0);
    vec3 substance = vec3(0.0);
    float metalness = uMetalness;
    float roughness = uRoughness;

    vec3 gradient = vec3(1.0);
    vec3 dx = vec3(gradOffset * scaleVol.x, 0.0, 0.0);
    vec3 dy = vec3(0.0, gradOffset * scaleVol.y, 0.0);
    vec3 dz = vec3(0.0, 0.0, gradOffset * scaleVol.z);

    float maxDist = min(vBoundingSphere.w * 2.0, uFar - uNear);
    float maxDistSq = maxDist * maxDist;

    for (int i = 0; i < uMaxSteps; ++i) {
        // break when beyond bounding-sphere or far-plane
        vec3 distVec = startLoc - pos;
        if (dot(distVec, distVec) > maxDistSq) break;

        unitPos = v3m4(pos, cartnToUnit);

        // continue when outside of grid
        if (unitPos.x > posMax.x || unitPos.y > posMax.y || unitPos.z > posMax.z ||
            unitPos.x < posMin.x || unitPos.y < posMin.y || unitPos.z < posMin.z
        ) {
            if (hit) break;
            prevValue = value;
            pos += step;
            continue;
        }

        cell = textureVal(unitPos);
        value = cell.a; // current voxel value

        if (uJumpLength > 0.0 && value < 0.01) {
            nextPos = pos + rayDir * uJumpLength;
            nextValue = textureVal(v3m4(nextPos, cartnToUnit)).a;
            if (nextValue < 0.01) {
                prevValue = nextValue;
                pos = nextPos;
                continue;
            }
        }

        #if defined(dRenderMode_isosurface)
            if (prevValue > 0.0 && ( // there was a prev Value
                (prevValue < uIsoValue.x && value > uIsoValue.x) || // entering isosurface
                (prevValue > uIsoValue.x && value < uIsoValue.x) // leaving isosurface
            )) {
                isoPos = v3m4(mix(pos - step, pos, ((prevValue - uIsoValue.x) / ((prevValue - uIsoValue.x) - (value - uIsoValue.x)))), cartnToUnit);

                vec4 mvPosition = modelViewTransform * vec4(isoPos * uGridDim, 1.0);

                #if defined(dClipVariant_pixel) && dClipObjectCount != 0
                    vec3 vModelPosition = v3m4(isoPos * uGridDim, modelTransform);
                    if (clipTest(vec4(vModelPosition, 0.0), 0)) {
                        prevValue = value;
                        pos += step;
                        continue;
                    }
                #endif

                float depth = calcDepth(mvPosition.xyz);
                if (depth > getDepth(gl_FragCoord.xy / uDrawingBufferSize))
                    break;

                #ifdef enabledFragDepth
                    if (!hit) {
                        gl_FragDepthEXT = depth;
                    }
                #endif

                #if defined(dRenderVariant_pickObject)
                    return vec4(encodeFloatRGB(float(uObjectId)), 1.0);
                #elif defined(dRenderVariant_pickInstance)
                    return vec4(encodeFloatRGB(vInstance), 1.0);
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
                        return packDepthToRGBA(depth);
                    #endif
                #elif defined(dRenderVariant_color)
                    #ifdef dPackedGroup
                        float group = decodeFloatRGB(textureGroup(floor(isoPos * uGridDim + 0.5) / uGridDim).rgb);
                    #else
                        vec3 g = floor(isoPos * uGridDim + 0.5);
                        float group = g.z + g.y * uGridDim.z + g.x * uGridDim.z * uGridDim.y;
                    #endif

                    #if defined(dColorType_uniform)
                        color = uColor;
                    #elif defined(dColorType_instance)
                        color = readFromTexture(tColor, vInstance, uColorTexDim).rgb;
                    #elif defined(dColorType_group)
                        color = readFromTexture(tColor, group, uColorTexDim).rgb;
                    #elif defined(dColorType_groupInstance)
                        color = readFromTexture(tColor, vInstance * float(uGroupCount) + group, uColorTexDim).rgb;
                    #elif defined(dColorType_vertex)
                        color = texture3dFrom1dTrilinear(tColor, isoPos, uGridDim, uColorTexDim, 0.0).rgb;
                    #elif defined(dColorType_vertexInstance)
                        color = texture3dFrom1dTrilinear(tColor, isoPos, uGridDim, uColorTexDim, vInstance * float(uVertexCount)).rgb;
                    #endif

                    #ifdef dOverpaint
                        #if defined(dOverpaintType_groupInstance)
                            overpaint = readFromTexture(tOverpaint, vInstance * float(uGroupCount) + group, uOverpaintTexDim);
                        #elif defined(dOverpaintType_vertexInstance)
                            overpaint = texture3dFrom1dTrilinear(tOverpaint, isoPos, uGridDim, uOverpaintTexDim, vInstance * float(uVertexCount));
                        #endif

                        color = mix(color, overpaint.rgb, overpaint.a);
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
                        vec3 normal = -normalize(normalMatrix * normalize(gradient));
                        normal = normal * (float(flipped) * 2.0 - 1.0);
                        normal = normal * -(float(interior) * 2.0 - 1.0);
                        #ifdef dSubstance
                            #if defined(dSubstanceType_groupInstance)
                                substance = readFromTexture(tSubstance, vInstance * float(uGroupCount) + group, uSubstanceTexDim).rgb;
                            #elif defined(dSubstanceType_vertexInstance)
                                substance = texture3dFrom1dTrilinear(tSubstance, isoPos, uGridDim, uSubstanceTexDim, vInstance * float(uVertexCount)).rgb;
                            #endif
                            metalness = mix(metalness, substance.r, substance.b);
                            roughness = mix(roughness, substance.g, substance.b);
                        #endif
                        #include apply_light_color
                    #endif

                    #if defined(dMarkerType_uniform)
                        float marker = uMarker;
                    #elif defined(dMarkerType_groupInstance)
                        float marker = readFromTexture(tMarker, vInstance * float(uGroupCount) + group, uMarkerTexDim).a;
                        marker = floor(marker * 255.0 + 0.5); // rounding required to work on some cards on win
                    #endif
                    #include apply_interior_color
                    #include apply_marker_color

                    preFogAlphaBlended = (1.0 - preFogAlphaBlended) * gl_FragColor.a + preFogAlphaBlended;
                    fragmentDepth = depth;
                    #include apply_fog

                    src = gl_FragColor;

                    if (!uTransparentBackground) {
                        // done in 'apply_fog' otherwise
                        src.rgb *= src.a;
                    }
                    dst = (1.0 - dst.a) * src + dst; // standard blending
                #endif

                #ifdef dSingleLayer
                    break;
                #endif

                hit = true;
            }
            prevValue = value;
        #elif defined(dRenderMode_volume)
            vec4 mvPosition = modelViewTransform * vec4(unitPos * uGridDim, 1.0);
            if (calcDepth(mvPosition.xyz) > getDepth(gl_FragCoord.xy / uDrawingBufferSize))
                break;

            #if defined(dClipVariant_pixel) && dClipObjectCount != 0
                vec3 vModelPosition = v3m4(unitPos * uGridDim, modelTransform);
                if (clipTest(vec4(vModelPosition, 0.0), 0)) {
                    prevValue = value;
                    pos += step;
                    continue;
                }
            #endif

            #if defined(dRenderVariant_color)
                vec3 vViewPosition = mvPosition.xyz;
                vec4 material = transferFunction(value);

                #ifdef dIgnoreLight
                    gl_FragColor.rgb = material.rgb;
                #else
                    if (material.a >= 0.01) {
                        #ifdef dPackedGroup
                            // compute gradient by central differences
                            gradient.x = textureVal(unitPos - dx).a - textureVal(unitPos + dx).a;
                            gradient.y = textureVal(unitPos - dy).a - textureVal(unitPos + dy).a;
                            gradient.z = textureVal(unitPos - dz).a - textureVal(unitPos + dz).a;
                        #else
                            gradient = cell.xyz * 2.0 - 1.0;
                        #endif
                        vec3 normal = -normalize(normalMatrix * normalize(gradient));
                        #include apply_light_color
                    } else {
                        gl_FragColor.rgb = material.rgb;
                    }
                #endif

                gl_FragColor.a = material.a * uAlpha * uTransferScale;

                #if defined(dMarkerType_uniform)
                    float marker = uMarker;
                #elif defined(dMarkerType_groupInstance)
                    #ifdef dPackedGroup
                        float group = decodeFloatRGB(textureGroup(floor(unitPos * uGridDim + 0.5) / uGridDim).rgb);
                    #else
                        vec3 g = floor(unitPos * uGridDim + 0.5);
                        float group = g.z + g.y * uGridDim.z + g.x * uGridDim.z * uGridDim.y;
                    #endif
                    float marker = readFromTexture(tMarker, vInstance * float(uGroupCount) + group, uMarkerTexDim).a;
                    marker = floor(marker * 255.0 + 0.5); // rounding required to work on some cards on win
                #endif
                #include apply_marker_color

                preFogAlphaBlended = (1.0 - preFogAlphaBlended) * gl_FragColor.a + preFogAlphaBlended;
                fragmentDepth = calcDepth(mvPosition.xyz);
                #include apply_fog

                src = gl_FragColor;

                if (!uTransparentBackground) {
                    // done in 'apply_fog' otherwise
                    src.rgb *= src.a;
                }
                dst = (1.0 - dst.a) * src + dst; // standard blending
            #endif
        #endif

        // break if the color is opaque enough
        if (dst.a > 0.95)
            break;

        pos += step;
    }

    #if defined(dRenderMode_isosurface) && defined(enabledFragDepth)
        // ensure depth is written everywhere
        if (!hit)
            gl_FragDepthEXT = 1.0;
    #endif

    return dst;
}

// TODO: support float texture for higher precision values???
// TODO: support clipping exclusion texture support

void main() {
    if (gl_FrontFacing)
        discard;

    #ifdef dRenderVariant_marking
        // not supported
        discard;
    #endif

    #if defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
        #if defined(dRenderMode_volume)
            // always ignore pick & depth for volume
            discard;
        #elif defined(dRenderMode_isosurface)
            if (uAlpha < uPickingAlphaThreshold)
                discard; // ignore so the element below can be picked
        #endif
    #endif

    vec3 rayDir = mix(normalize(vOrigPos - uCameraPosition), uCameraDir, uIsOrtho);
    vec3 step = rayDir * uStepScale;

    float boundingSphereNear = distance(vBoundingSphere.xyz, uCameraPosition) - vBoundingSphere.w;
    float d = max(uNear, boundingSphereNear) - mix(0.0, distance(vOrigPos, uCameraPosition), uIsOrtho);
    vec3 start = mix(uCameraPosition, vOrigPos, uIsOrtho) + (d * rayDir);
    gl_FragColor = raymarch(start, step, rayDir);

    #if defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
        // discard when nothing was hit
        if (gl_FragColor == vec4(0.0))
            discard;
    #endif

    #if defined(dRenderVariant_color)
        #if defined(dRenderMode_isosurface) && defined(enabledFragDepth)
            float fragmentDepth = gl_FragDepthEXT;
        #else
            float fragmentDepth = calcDepth((uModelView * vec4(start, 1.0)).xyz);
        #endif
        float preFogAlpha = clamp(preFogAlphaBlended, 0.0, 1.0);
        interior = false;
        #include wboit_write
    #endif
}
`;