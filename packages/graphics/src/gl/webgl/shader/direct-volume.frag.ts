/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { common } from './chunks/common.glsl.js';
import { light_frag_params } from './chunks/light-frag-params.glsl.js';
import { common_clip } from './chunks/common-clip.glsl.js';
import { read_from_texture } from './chunks/read-from-texture.glsl.js';
import { texture3d_from_1d_trilinear } from './chunks/texture3d-from-1d-trilinear.glsl.js';
import { texture3d_from_2d_nearest } from './chunks/texture3d-from-2d-nearest.glsl.js';
import { texture3d_from_2d_linear } from './chunks/texture3d-from-2d-linear.glsl.js';
import { apply_light_color } from './chunks/apply-light-color.glsl.js';
import { apply_marker_color } from './chunks/apply-marker-color.glsl.js';
import { apply_fog } from './chunks/apply-fog.glsl.js';
import { wboit_write } from './chunks/wboit-write.glsl.js';
import { apply_marking } from './chunks/apply-marking.glsl.js';

export const directVolume_frag = `
precision highp float;
precision highp int;

${common}
${light_frag_params}

// per-group clip flags come from a varying, which direct volumes don't have
#undef dClipping
${common_clip}

${read_from_texture}
${texture3d_from_1d_trilinear}
${texture3d_from_2d_nearest}
${texture3d_from_2d_linear}

varying vec3 vModelPosition;
varying float vInstance;
varying vec4 vBoundingSphere;
varying mat4 vTransform;
varying mat3 vNormalMatrix;

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
    return clipZToDepth(clipZW.x / clipZW.y);
}

float transferFunction(float value) {
    return texture2D(tTransferTex, vec2(value, 0.0)).r;
}

float getDepth(const in vec2 coords) {
    #ifdef depthTextureSupport
        return texture2D(tDepth, coords).r;
    #else
        return unpackRGBAToDepth(texture2D(tDepth, coords));
    #endif
}

const float gradOffset = 0.5;

/**
 * The offset of a pixel's samples in [0, 1) steps, against the banding of samples at the same depths: interleaved
 * gradient noise (Jimenez 2014), varied by the projection's off-center terms (the sub-pixel offset of multi-sampling,
 * 0 otherwise) so that multi-sampling averages the noise out
 */
float rayJitter(const in vec2 fragCoord) {
    float noise = fract(52.9829189 * fract(dot(fragCoord, vec2(0.06711056, 0.00583715))));
    vec4 offset = vec4(uProjection[2].xy, uProjection[3].xy);
    return fract(noise + fract(sin(dot(offset, vec4(12.9898, 78.233, 37.719, 4.581))) * 43758.5453));
}

// the edge of a brick in cells, \`BrickSize\` of the direct-volume geometry
const float BrickSize = 8.0;

vec3 v3m4(vec3 p, mat4 m) {
    return (m * vec4(p, 1.0)).xyz;
}

float preFogAlphaBlended = 0.0;

// the surface of the marking, pick and depth variants: where the opacity of the samples the pass looks at (all, but
// for marking: unmarked ones in its depth pass, marked ones in its mask pass) first reaches \`uSurfaceOpacity\` (the
// \`surfaceOpacity\` param, low so that it is where the volume looks solid)
float surfaceDepth = -1.0;
float surfaceMarker = 0.0;
float surfaceGroup = 0.0;

/** 1 / v, with components close to 0 as a large positive number (a ray parallel to a slab) */
vec3 safeInverse(const in vec3 v) {
    return vec3(
        abs(v.x) > 1e-8 ? 1.0 / v.x : 1e8,
        abs(v.y) > 1e-8 ? 1.0 / v.y : 1e8,
        abs(v.z) > 1e-8 ? 1.0 / v.z : 1e8
    );
}

vec4 raymarch(vec3 startLoc, vec3 step, vec3 rayDir) {
    mat4 cartnToUnit = uCartnToUnit * inverse4(vTransform);
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        mat4 modelTransform = uModel * vTransform * uTransform;
    #endif
    mat4 modelViewTransform = uModelView * vTransform * uTransform;

    vec3 scaleVol = vec3(1.0) / uGridDim;
    vec4 cell;
    float value = 0.0;
    vec4 src = vec4(0.0);
    vec4 dst = vec4(0.0);
    float fragmentDepth;

    vec3 posMin = vec3(0.0);
    vec3 posMax = vec3(1.0) - vec3(1.0) / uGridDim;

    // the samples of the ray, in unit grid coordinates and in view space, are affine in the step \`t\`
    vec3 unitStart = v3m4(startLoc / uModelScale, cartnToUnit);
    vec3 unitStep = (cartnToUnit * vec4(step / uModelScale, 0.0)).xyz;
    vec3 viewStart = (modelViewTransform * vec4(unitStart * uGridDim, 1.0)).xyz;
    vec3 viewStep = (modelViewTransform * vec4(unitStep * uGridDim, 0.0)).xyz;

    // the steps inside the grid's box (slab test)
    vec3 invStep = safeInverse(unitStep);
    vec3 t0 = (posMin - unitStart) * invStep;
    vec3 t1 = (posMax - unitStart) * invStep;
    vec3 slabNear = min(t0, t1);
    vec3 slabFar = max(t0, t1);
    float enterT = max(max(slabNear.x, slabNear.y), slabNear.z);
    float endT = min(min(slabFar.x, slabFar.y), slabFar.z);

    // before the far plane
    float stepLength = length(step);
    endT = min(endT, (uFar - uNear) / stepLength);

    // before the opaque depth (depth is reversed: a sample is behind it where its view z is smaller)
    float opaqueViewZ = depthToViewZ(uIsOrtho, getDepth(gl_FragCoord.xy / uDrawingBufferSize), uNear, uFar);
    if (viewStep.z < 0.0) {
        endT = min(endT, (opaqueViewZ - viewStart.z) / viewStep.z);
    } else if (viewStart.z < opaqueViewZ) {
        return dst;
    }

    #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
        // the same phase for neighboring rays: a coherent surface, its edges not frayed by noise
        float jitter = 0.5;
    #else
        float jitter = rayJitter(gl_FragCoord.xy);
    #endif
    float t = max(0.0, ceil(enterT - jitter)) + jitter;
    float jumpSteps = uJumpLength / stepLength;

    vec3 unitPos;
    float nextT;
    float currentBrick = -1.0;
    bool brickVisible = true;

    vec4 material;
    vec4 overpaint;
    float metalness = uMetalness;
    float roughness = uRoughness;
    float emissive = uEmissive;

    vec3 gradient = vec3(1.0);
    vec3 dx = vec3(gradOffset * scaleVol.x, 0.0, 0.0);
    vec3 dy = vec3(0.0, gradOffset * scaleVol.y, 0.0);
    vec3 dz = vec3(0.0, 0.0, gradOffset * scaleVol.z);

    #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
        float surfaceOpacity = 0.0;
    #endif

    for (int i = 0; i < uMaxSteps; ++i) {
        if (t > endT) break;

        unitPos = unitStart + t * unitStep;

        // the rest of a brick the transfer function makes transparent is skipped
        if (uBrickDim.x > 0.0) {
            vec3 brick = floor(unitPos * uGridDim / BrickSize);
            float brickIndex = brick.x + uBrickDim.x * (brick.y + uBrickDim.y * brick.z);
            if (brickIndex != currentBrick) {
                currentBrick = brickIndex;
                brickVisible = readFromTexture(tBrickTex, brickIndex, uBrickTexDim).r > 0.0;
            }
            if (!brickVisible) {
                // the first step past the brick's exit face (a step on the face is tested again)
                vec3 exitFace = vec3(
                    invStep.x > 0.0 ? brick.x + 1.0 : brick.x,
                    invStep.y > 0.0 ? brick.y + 1.0 : brick.y,
                    invStep.z > 0.0 ? brick.z + 1.0 : brick.z
                ) * BrickSize / uGridDim;
                vec3 exits = (exitFace - unitStart) * invStep;
                t = max(t + 1.0, ceil(min(min(exits.x, exits.y), exits.z) - jitter - 1e-3) + jitter);
                continue;
            }
        }

        cell = textureVal(unitPos);
        value = cell.a; // current voxel value

        if (uJumpLength > 0.0 && value < 0.01) {
            nextT = t + jumpSteps;
            if (textureVal(unitStart + nextT * unitStep).a < 0.01) {
                t = nextT;
                continue;
            }
        }

        // transparent samples add nothing (markers don't lift them either)
        material.a = transferFunction(value);
        if (material.a <= 0.0) {
            t += 1.0;
            continue;
        }

        #if defined(dClipVariant_pixel) && defined(dClipObjects)
            vec3 modelPosition = v3m4(unitPos * uGridDim, modelTransform);
            if (clipTest(modelPosition)) {
                t += 1.0;
                continue;
            }
        #endif

        vec3 vViewPosition = viewStart + t * viewStep;

        #ifdef dPackedGroup
            // a cell without a group (packed 0: no atom within reach, density about 0) takes group 0, as the clamped
            // reads of index -1 did
            float group = max(unpackRGBToInt(textureGroup(floor(unitPos * uGridDim + 0.5) / uGridDim).rgb), 0.0);
        #else
            vec3 g = floor(unitPos * uGridDim + 0.5);
            // note that we swap x and z because the texture is flipped around y
            #if defined(dAxisOrder_012)
                float group = g.z + g.y * uGridDim.z + g.x * uGridDim.z * uGridDim.y; // 210
            #elif defined(dAxisOrder_021)
                float group = g.y + g.z * uGridDim.y + g.x * uGridDim.y * uGridDim.z; // 120
            #elif defined(dAxisOrder_102)
                float group = g.z + g.x * uGridDim.z + g.y * uGridDim.z * uGridDim.x; // 201
            #elif defined(dAxisOrder_120)
                float group = g.x + g.z * uGridDim.x + g.y * uGridDim.x * uGridDim.z; // 021
            #elif defined(dAxisOrder_201)
                float group = g.y + g.x * uGridDim.y + g.z * uGridDim.y * uGridDim.x; // 102
            #elif defined(dAxisOrder_210)
                float group = g.x + g.y * uGridDim.x + g.z * uGridDim.x * uGridDim.y; // 012
            #endif
        #endif

        #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
            #if defined(dRenderVariant_marking)
                float sampleMarker = uMarker;
                if (uMarker == -1.0) {
                    sampleMarker = read_tMarker((uInstanceGranularity ? vInstance : vInstance * float(uGroupCount) + group)).r;
                    sampleMarker = floor(sampleMarker * 255.0 + 0.5);
                }
                bool surfaceSample = (uMarkingType == 1) == (sampleMarker == 0.0);
            #else
                float sampleMarker = 0.0;
                bool surfaceSample = true;
            #endif
            if (surfaceSample) {
                surfaceOpacity += (1.0 - surfaceOpacity) * (1.0 - pow(1.0 - clamp(material.a * uAlpha, 0.0, 1.0), uTransferScale));
                if (surfaceOpacity >= uSurfaceOpacity) {
                    surfaceDepth = calcDepth(vViewPosition);
                    surfaceMarker = sampleMarker;
                    surfaceGroup = group;
                    break;
                }
            }
            t += 1.0;
        #else

        #if defined(dColorType_direct) && defined(dUsePalette)
            float paletteValue = (value - uPaletteDomain[0]) / (uPaletteDomain[1] - uPaletteDomain[0]);
            material.rgb = texture2D(tPalette, vec2(clamp(paletteValue, 0.0, 1.0), 0.0)).rgb;
        #elif defined(dColorType_uniform)
            material.rgb = uColor;
        #elif defined(dColorType_texture)
            // by the color type (\`ColorTypeIndex\`): instance, group, groupInstance, vertex, vertexInstance
            if (uColorIndexType == 0) {
                material.rgb = read_tColor(vInstance).rgb;
            } else if (uColorIndexType == 1) {
                material.rgb = read_tColor(group).rgb;
            } else if (uColorIndexType == 2) {
                material.rgb = read_tColor(vInstance * float(uGroupCount) + group).rgb;
            } else if (uColorIndexType == 3) {
                material.rgb = texture3dFrom1dTrilinear(0, unitPos, uGridDim, 0.0).rgb;
            } else {
                material.rgb = texture3dFrom1dTrilinear(0, unitPos, uGridDim, vInstance * float(uVertexCount)).rgb;
            }
        #endif

        #ifdef dOverpaint
            #if defined(dOverpaintType_texture)
                overpaint = read_tOverpaint((uInstanceGranularity ? vInstance : vInstance * float(uGroupCount) + group));
            #elif defined(dOverpaintType_vertexInstance)
                overpaint = texture3dFrom1dTrilinear(1, unitPos, uGridDim, vInstance * float(uVertexCount));
            #endif

            material.rgb = mix(material.rgb, overpaint.rgb, overpaint.a);
        #endif

        #if defined(dIgnoreLight)
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
                // a zero gradient (byte normals filtered halfway between 127 and 128) faces the camera, not NaN
                vec3 normal = dot(gradient, gradient) > 1e-12 ? -normalize(vNormalMatrix * gradient) : vec3(0.0, 0.0, 1.0);
                ${apply_light_color}
            } else {
                gl_FragColor.rgb = material.rgb;
            }
        #endif

        // the transfer function's opacity is per cell: a step of 1 / stepsPerCell cells (\`uTransferScale\`) is as opaque
        // as that fraction of a cell, exactly (the linear \`a * uTransferScale\` is only close for small opacities)
        gl_FragColor.a = 1.0 - pow(1.0 - clamp(material.a * uAlpha, 0.0, 1.0), uTransferScale);

        #if defined(dColorMarker)
            float marker = uMarker;
            if (uMarker == -1.0) {
                marker = read_tMarker((uInstanceGranularity ? vInstance : vInstance * float(uGroupCount) + group)).r;
                marker = floor(marker * 255.0 + 0.5); // rounding required to work on some cards on win
            }
        #endif
        ${apply_marker_color}

        preFogAlphaBlended = (1.0 - preFogAlphaBlended) * gl_FragColor.a + preFogAlphaBlended;
        fragmentDepth = calcDepth(vViewPosition);
        ${apply_fog}

        src = gl_FragColor;

        if (!uTransparentBackground || !uFog) {
            // done in 'apply_fog' otherwise
            src.rgb *= src.a;
        }
        dst = (1.0 - dst.a) * src + dst; // standard blending

        // break if the color is opaque enough
        if (dst.a > 0.95)
            break;

        t += 1.0;
        #endif
    }

    return dst;
}

// TODO: support float texture for higher precision values???
// TODO: support clipping exclusion texture support

void main() {
    #if defined(dRenderVariant_tracing) || defined(dRenderVariant_emissive)
        discard;
    #else
        if (gl_FrontFacing)
            discard;

        vec3 rayDir = mix(normalize(vModelPosition - uCameraPosition), uCameraDir, uIsOrtho);
        vec3 step = rayDir * uStepScale * uModelScale;

        float boundingSphereNear = distance(vBoundingSphere.xyz, uCameraPosition) - vBoundingSphere.w;
        float d = max(uNear, boundingSphereNear) - mix(0.0, distance(vModelPosition, uCameraPosition), uIsOrtho);
        vec3 start = mix(uCameraPosition, vModelPosition, uIsOrtho) + (d * rayDir);

        #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick) || defined(dRenderVariant_depth)
            #if defined(enabledFragDepth)
                raymarch(start, step, rayDir);
                if (surfaceDepth < 0.0) discard;
                float fragmentDepth = surfaceDepth;
                gl_FragDepthEXT = fragmentDepth;
                #if defined(dRenderVariant_marking)
                    float marker = surfaceMarker;
                    vec4 material = vec4(0.0, 0.0, 0.0, 1.0);
                    ${apply_marking}
                    gl_FragColor = material;
                #elif defined(dRenderVariant_pick)
                    // too deep in the fog to be picked, as other objects (\`check_picking_alpha\`)
                    float viewZ = depthToViewZ(uIsOrtho, fragmentDepth, uNear, uFar);
                    if ((1.0 - smoothstep(uFogNear, uFogFar, abs(viewZ))) * uAlpha < 0.1) discard;
                    #ifdef requiredDrawBuffers
                        gl_FragColor = vec4(packIntToRGB(float(uObjectId)), 1.0);
                        gl_FragData[1] = vec4(packIntToRGB(vInstance), 1.0);
                        gl_FragData[2] = vec4(packIntToRGB(surfaceGroup), 1.0);
                        gl_FragData[3] = packDepthToRGBA(fragmentDepth);
                    #else
                        if (uPickType == 1) {
                            gl_FragColor = vec4(packIntToRGB(float(uObjectId)), 1.0);
                        } else if (uPickType == 2) {
                            gl_FragColor = vec4(packIntToRGB(vInstance), 1.0);
                        } else {
                            gl_FragColor = vec4(packIntToRGB(surfaceGroup), 1.0);
                        }
                    #endif
                #else
                    gl_FragColor = packDepthToRGBA(fragmentDepth);
                #endif
            #else
                // the surface's depth can't be written
                discard;
            #endif
        #else
            gl_FragColor = raymarch(start, step, rayDir);

            float fragmentDepth = calcDepth((uView * vec4(start, 1.0)).xyz);
            float preFogAlpha = clamp(preFogAlphaBlended, 0.0, 1.0);
            ${wboit_write}
        #endif
    #endif
}
`;
