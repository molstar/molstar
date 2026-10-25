/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { clip_functions_wgsl } from './chunks/clip.wgsl.js';
import { common_wgsl } from './chunks/common.wgsl.js';
import { light_wgsl } from './chunks/light.wgsl.js';
import { fragment_output_wgsl, marking_output_wgsl, wboit_output_wgsl } from './chunks/output.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { texture_sample_wgsl } from './chunks/quad.wgsl.js';
import { read_from_texture_wgsl } from './chunks/read-from-texture.wgsl.js';

/**
 * Direct volume rendering: the back faces of the grid box, ray marched front to back through the 3d grid texture
 * (\`tGridTex\`, values in alpha, gradients in rgb) with a transfer function, lighting, markers and fog; color variants
 * only, as on WebGL
 */
export const direct_volume_wgsl = /* wgsl */ `
${common_wgsl}
${read_from_texture_wgsl}
${light_wgsl}
${clip_functions_wgsl}
${packing_wgsl}
${texture_sample_wgsl('tTransferTex')}
#if defined(dUsePalette)
    ${texture_sample_wgsl('tPalette')}
#endif

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) vModelPosition: vec3<f32>,
    @location(1) @interpolate(flat) vInstance: f32,
    @location(2) @interpolate(flat) vBoundingSphere: vec4<f32>,
    @location(3) @interpolate(flat) vTransform0: vec4<f32>,
    @location(4) @interpolate(flat) vTransform1: vec4<f32>,
    @location(5) @interpolate(flat) vTransform2: vec4<f32>,
    @location(6) @interpolate(flat) vTransform3: vec4<f32>,
    // the normal matrix per instance: per fragment it costs as if it were computed per lit sample
    @location(7) @interpolate(flat) vNormalMatrix0: vec3<f32>,
    @location(8) @interpolate(flat) vNormalMatrix1: vec3<f32>,
    @location(9) @interpolate(flat) vNormalMatrix2: vec3<f32>,
};

${fragment_output_wgsl(`#if defined(dRenderVariant_marking) || defined(dRenderVariant_pick)
    @builtin(frag_depth) fragDepth: f32,
#endif`)}

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let aTransform = mat4x4<f32>(in.aTransform0, in.aTransform1, in.aTransform2, in.aTransform3);
    let unitCoord = vec4<f32>(in.aPosition + vec3<f32>(0.5), 1.0);
    let mvPosition = uGlobal.uModelView * aTransform * uObject.uUnitToCartn * unitCoord;

    out.vModelPosition = (uGlobal.uModel * aTransform * uObject.uUnitToCartn * unitCoord).xyz;
    out.vInstance = in.aInstance;
    let sphere = uObject.uInvariantBoundingSphere;
    out.vBoundingSphere = vec4<f32>((uGlobal.uModel * aTransform * vec4<f32>(sphere.xyz, 1.0)).xyz, uGlobal.uModelScale * sphere.w);
    out.vTransform0 = in.aTransform0;
    out.vTransform1 = in.aTransform1;
    out.vTransform2 = in.aTransform2;
    out.vTransform3 = in.aTransform3;
    let normalMatrix = adjoint(uGlobal.uModelView * aTransform);
    out.vNormalMatrix0 = normalMatrix[0];
    out.vNormalMatrix1 = normalMatrix[1];
    out.vNormalMatrix2 = normalMatrix[2];

    out.position = uGlobal.uProjection * mvPosition;
    // close to the far plane (but not too close to get precision issues), the window depth of WebGL's \`-w + 0.01\` (depth is reversed)
    out.position.z = 0.005;
    return out;
}

fn inverse4(m: mat4x4<f32>) -> mat4x4<f32> {
    let a00 = m[0][0]; let a01 = m[0][1]; let a02 = m[0][2]; let a03 = m[0][3];
    let a10 = m[1][0]; let a11 = m[1][1]; let a12 = m[1][2]; let a13 = m[1][3];
    let a20 = m[2][0]; let a21 = m[2][1]; let a22 = m[2][2]; let a23 = m[2][3];
    let a30 = m[3][0]; let a31 = m[3][1]; let a32 = m[3][2]; let a33 = m[3][3];

    let b00 = a00 * a11 - a01 * a10;
    let b01 = a00 * a12 - a02 * a10;
    let b02 = a00 * a13 - a03 * a10;
    let b03 = a01 * a12 - a02 * a11;
    let b04 = a01 * a13 - a03 * a11;
    let b05 = a02 * a13 - a03 * a12;
    let b06 = a20 * a31 - a21 * a30;
    let b07 = a20 * a32 - a22 * a30;
    let b08 = a20 * a33 - a23 * a30;
    let b09 = a21 * a32 - a22 * a31;
    let b10 = a21 * a33 - a23 * a31;
    let b11 = a22 * a33 - a23 * a32;

    let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

    return mat4x4<f32>(
        a11 * b11 - a12 * b10 + a13 * b09,
        a02 * b10 - a01 * b11 - a03 * b09,
        a31 * b05 - a32 * b04 + a33 * b03,
        a22 * b04 - a21 * b05 - a23 * b03,
        a12 * b08 - a10 * b11 - a13 * b07,
        a00 * b11 - a02 * b08 + a03 * b07,
        a32 * b02 - a30 * b05 - a33 * b01,
        a20 * b05 - a22 * b02 + a23 * b01,
        a10 * b10 - a11 * b08 + a13 * b06,
        a01 * b08 - a00 * b10 - a03 * b06,
        a30 * b04 - a31 * b02 + a33 * b00,
        a21 * b02 - a20 * b04 - a23 * b00,
        a11 * b07 - a10 * b09 - a12 * b06,
        a00 * b09 - a01 * b07 + a02 * b06,
        a31 * b01 - a30 * b03 - a32 * b00,
        a20 * b03 - a21 * b01 + a22 * b00
    ) * (1.0 / det);
}

/** value (alpha) and gradient (rgb) at \`pos\` in unit grid coordinates, filtered where the format allows */
fn textureVal(pos: vec3<f32>) -> vec4<f32> {
    let p = pos + (vec3<f32>(0.5) / uObject.uGridDim);
    #if defined(dLinear_tGridTex)
        return textureSampleLevel(tGridTex, sLinear, p, 0.0);
    #else
        let size = vec3<i32>(textureDimensions(tGridTex));
        return textureLoad(tGridTex, clamp(vec3<i32>(floor(p * vec3<f32>(size))), vec3<i32>(0), size - 1), 0);
    #endif
}

fn textureGroup(pos: vec3<f32>) -> vec4<f32> {
    return textureLoad(tGridTex, vec3<i32>(pos * uObject.uGridDim), 0);
}

/** element \`i\` of the data layer \`layer\`: 0 \`tColor\`, 1 \`tOverpaint\` (generated accessors) */
fn readLayer(layer: u32, i: i32) -> vec4<f32> {
    if (layer == 0u) { return read_tColor(i); }
    return read_tOverpaint(i);
}

/** trilinear interpolation of per-vertex data (stored z fastest) of the data layer \`layer\` from offset \`offset\` */
fn texture3dFrom1dTrilinear(layer: u32, pos: vec3<f32>, gridDim: vec3<f32>, offset: f32) -> vec4<f32> {
    let gdYZ = gridDim.z * gridDim.y;
    let gdZ = gridDim.z;
    let p0 = floor(pos * gridDim);
    let p1 = ceil(pos * gridDim);
    // not \`(pos * gridDim - p0) / (p1 - p0)\`: 0 / 0 on grid coordinates
    let pd = pos * gridDim - p0;
    let s000 = readLayer(layer, i32(offset + p0.z + p0.y * gdZ + p0.x * gdYZ));
    let s100 = readLayer(layer, i32(offset + p0.z + p0.y * gdZ + p1.x * gdYZ));
    let s001 = readLayer(layer, i32(offset + p1.z + p0.y * gdZ + p0.x * gdYZ));
    let s101 = readLayer(layer, i32(offset + p1.z + p0.y * gdZ + p1.x * gdYZ));
    let s010 = readLayer(layer, i32(offset + p0.z + p1.y * gdZ + p0.x * gdYZ));
    let s110 = readLayer(layer, i32(offset + p0.z + p1.y * gdZ + p1.x * gdYZ));
    let s011 = readLayer(layer, i32(offset + p1.z + p1.y * gdZ + p0.x * gdYZ));
    let s111 = readLayer(layer, i32(offset + p1.z + p1.y * gdZ + p1.x * gdYZ));
    let s00 = mix(s000, s100, pd.x);
    let s01 = mix(s001, s101, pd.x);
    let s10 = mix(s010, s110, pd.x);
    let s11 = mix(s011, s111, pd.x);
    let s0 = mix(s00, s10, pd.y);
    let s1 = mix(s01, s11, pd.y);
    return mix(s0, s1, pd.z);
}

fn transferFunction(value: f32) -> f32 {
    return sample_tTransferTex(vec2<f32>(value, 0.0)).r;
}

fn v3m4(p: vec3<f32>, m: mat4x4<f32>) -> vec3<f32> {
    return (m * vec4<f32>(p, 1.0)).xyz;
}

const GradOffset = 0.5;

/**
 * The offset of a pixel's samples in [0, 1) steps, against the banding of samples at the same depths: interleaved
 * gradient noise (Jimenez 2014), varied by the projection's off-center terms (the sub-pixel offset of multi-sampling,
 * 0 otherwise) so that multi-sampling averages the noise out
 */
fn rayJitter(fragCoord: vec2<f32>) -> f32 {
    let noise = fract(52.9829189 * fract(dot(fragCoord, vec2<f32>(0.06711056, 0.00583715))));
    let projection = uGlobal.uProjection;
    let offset = vec4<f32>(projection[2].xy, projection[3].xy);
    return fract(noise + fract(sin(dot(offset, vec4<f32>(12.9898, 78.233, 37.719, 4.581))) * 43758.5453));
}

/** the edge of a brick in cells, \`BrickSize\` of the direct-volume geometry */
const BrickSize = 8.0;

var<private> preFogAlphaBlended = 0.0;

/**
 * the surface of the marking and pick variants: where the opacity of the samples the pass looks at (all, but for
 * marking: unmarked ones in its depth pass, marked ones in its mask pass) first reaches \`uSurfaceOpacity\` (the
 * \`surfaceOpacity\` param, low so that it is where the volume looks solid)
 */
var<private> surfaceDepth = -1.0;
var<private> surfaceMarker = 0.0;
var<private> surfaceGroup = 0.0;

/** 1 / v, with components close to 0 as a large positive number (a ray parallel to a slab) */
fn safeInverse(v: vec3<f32>) -> vec3<f32> {
    return select(vec3<f32>(1e8), vec3<f32>(1.0) / v, abs(v) > vec3<f32>(1e-8));
}

fn raymarch(startLoc: vec3<f32>, step: vec3<f32>, fragCoord: vec2<f32>, instance: f32, transform: mat4x4<f32>, normalMatrix: mat3x3<f32>) -> vec4<f32> {
    let cartnToUnit = uObject.uCartnToUnit * inverse4(transform);
    #if defined(dClipVariant_pixel) && defined(dClipObjects)
        let modelTransform = uGlobal.uModel * transform * uObject.uTransform;
    #endif
    let modelViewTransform = uGlobal.uModelView * transform * uObject.uTransform;

    let gridDim = uObject.uGridDim;
    let scaleVol = vec3<f32>(1.0) / gridDim;
    let modelScale = uGlobal.uModelScale;
    let groupCount = f32(uObject.uGroupCount);
    var value = 0.0;
    var dst = vec4<f32>(0.0);

    let posMin = vec3<f32>(0.0);
    let posMax = vec3<f32>(1.0) - vec3<f32>(1.0) / gridDim;

    // the samples of the ray, in unit grid coordinates and in view space, are affine in the step \`t\`
    let unitStart = v3m4(startLoc / modelScale, cartnToUnit);
    let unitStep = (cartnToUnit * vec4<f32>(step / modelScale, 0.0)).xyz;
    let viewStart = (modelViewTransform * vec4<f32>(unitStart * gridDim, 1.0)).xyz;
    let viewStep = (modelViewTransform * vec4<f32>(unitStep * gridDim, 0.0)).xyz;

    // the steps inside the grid's box (slab test)
    let invStep = safeInverse(unitStep);
    let t0 = (posMin - unitStart) * invStep;
    let t1 = (posMax - unitStart) * invStep;
    let slabNear = min(t0, t1);
    let slabFar = max(t0, t1);
    let enterT = max(max(slabNear.x, slabNear.y), slabNear.z);
    var endT = min(min(slabFar.x, slabFar.y), slabFar.z);

    // before the far plane
    let stepLength = length(step);
    endT = min(endT, (uGlobal.uFar - uGlobal.uNear) / stepLength);

    // before the opaque depth (depth is reversed: a sample is behind it where its view z is smaller)
    // at the fragment's position over the drawing buffer size: volumes may be drawn at a lower resolution than the depth
    let depthSize = vec2<f32>(textureDimensions(tDepth));
    let opaqueDepth = textureLoad(tDepth, vec2<i32>(fragCoord / uGlobal.uDrawingBufferSize * depthSize), 0);
    let opaqueViewZ = depthToViewZ(uGlobal.uIsOrtho, opaqueDepth, uGlobal.uNear, uGlobal.uFar);
    if (viewStep.z < 0.0) {
        endT = min(endT, (opaqueViewZ - viewStart.z) / viewStep.z);
    } else if (viewStart.z < opaqueViewZ) {
        return dst;
    }

    #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick)
        // the same phase for neighboring rays: a coherent surface, its edges not frayed by noise
        let jitter = 0.5;
    #else
        let jitter = rayJitter(fragCoord);
    #endif
    var t = max(0.0, ceil(enterT - jitter)) + jitter;
    let jumpSteps = uObject.uJumpLength / stepLength;
    let brickDim = uObject.uBrickDim;
    var currentBrick = -1;
    var brickVisible = true;

    var material = vec4<f32>(0.0);
    let metalness = uMaterial.uMetalness;
    let roughness = uMaterial.uRoughness;
    let emissive = uMaterial.uEmissive;

    let dx = vec3<f32>(GradOffset * scaleVol.x, 0.0, 0.0);
    let dy = vec3<f32>(0.0, GradOffset * scaleVol.y, 0.0);
    let dz = vec3<f32>(0.0, 0.0, GradOffset * scaleVol.z);

    #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick)
        var surfaceOpacity = 0.0;
    #endif

    let transparentBackground = uGlobal.uTransparentBackground != 0u;
    let fog = uGlobal.uFog != 0u;

    for (var i = 0; i < uObject.uMaxSteps; i++) {
        if (t > endT) {
            break;
        }

        let unitPos = unitStart + t * unitStep;

        // the rest of a brick the transfer function makes transparent is skipped
        if (brickDim.x > 0.0) {
            let brick = floor(unitPos * gridDim / BrickSize);
            let brickIndex = i32(brick.x + brickDim.x * (brick.y + brickDim.y * brick.z));
            if (brickIndex != currentBrick) {
                currentBrick = brickIndex;
                brickVisible = readFromTexture(tBrickTex, brickIndex).r > 0.0;
            }
            if (!brickVisible) {
                // the first step past the brick's exit face (a step on the face is tested again)
                let exitFace = select(brick, brick + 1.0, invStep > vec3<f32>(0.0)) * BrickSize / gridDim;
                let exits = (exitFace - unitStart) * invStep;
                t = max(t + 1.0, ceil(min(min(exits.x, exits.y), exits.z) - jitter - 1e-3) + jitter);
                continue;
            }
        }

        let cell = textureVal(unitPos);
        value = cell.a; // current voxel value

        if (uObject.uJumpLength > 0.0 && value < 0.01) {
            let nextT = t + jumpSteps;
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
            let modelPosition = v3m4(unitPos * gridDim, modelTransform);
            // direct volumes have no per-group clipping flags
            if (clipTest(modelPosition, 0)) {
                t += 1.0;
                continue;
            }
        #endif

        let viewPosition = viewStart + t * viewStep;

        #if defined(dPackedGroup)
            let group = unpackRGBToInt(textureGroup(floor(unitPos * gridDim + 0.5) / gridDim).rgb);
        #else
            let g = floor(unitPos * gridDim + 0.5);
            // x and z are swapped because the texture is flipped around y
            #if defined(dAxisOrder_012)
                let group = g.z + g.y * gridDim.z + g.x * gridDim.z * gridDim.y; // 210
            #elif defined(dAxisOrder_021)
                let group = g.y + g.z * gridDim.y + g.x * gridDim.y * gridDim.z; // 120
            #elif defined(dAxisOrder_102)
                let group = g.z + g.x * gridDim.z + g.y * gridDim.z * gridDim.x; // 201
            #elif defined(dAxisOrder_120)
                let group = g.x + g.z * gridDim.x + g.y * gridDim.x * gridDim.z; // 021
            #elif defined(dAxisOrder_201)
                let group = g.y + g.x * gridDim.y + g.z * gridDim.y * gridDim.x; // 102
            #else
                let group = g.x + g.y * gridDim.x + g.z * gridDim.x * gridDim.y; // 012
            #endif
        #endif
        // the marker and layer data per instance or group instance
        let dataIndex = select(i32(instance * groupCount + group), i32(instance), uObject.uInstanceGranularity != 0u);

        #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick)
            #if defined(dRenderVariant_marking)
                var sampleMarker = uObject.uMarker;
                if (uObject.uMarker == -1.0) {
                    sampleMarker = floor(read_tMarker(dataIndex).r * 255.0 + 0.5);
                }
                let surfaceSample = (uGlobal.uMarkingType == 1) == (sampleMarker == 0.0);
            #else
                let sampleMarker = 0.0;
                let surfaceSample = true;
            #endif
            if (surfaceSample) {
                surfaceOpacity += (1.0 - surfaceOpacity) * (1.0 - pow(1.0 - clamp(material.a * uMaterial.uAlpha, 0.0, 1.0), uMaterial.uTransferScale));
                if (surfaceOpacity >= uObject.uSurfaceOpacity) {
                    surfaceDepth = calcDepth(viewPosition);
                    surfaceMarker = sampleMarker;
                    surfaceGroup = group;
                    break;
                }
            }
            t += 1.0;
        #else

        #if defined(dColorType_direct) && defined(dUsePalette)
            let domain = uObject.uPaletteDomain;
            let paletteValue = (value - domain.x) / (domain.y - domain.x);
            material = vec4<f32>(sample_tPalette(vec2<f32>(clamp(paletteValue, 0.0, 1.0), 0.0)).rgb, material.a);
        #elif defined(dColorType_uniform)
            material = vec4<f32>(uMaterial.uColor, material.a);
        #elif defined(dColorType_texture)
            // by the color type (\`ColorTypeIndex\`): instance, group, groupInstance, vertex, vertexInstance
            switch uObject.uColorIndexType {
                case 1: { material = vec4<f32>(read_tColor(i32(group)).rgb, material.a); }
                case 2: { material = vec4<f32>(read_tColor(i32(instance * groupCount + group)).rgb, material.a); }
                case 3: { material = vec4<f32>(texture3dFrom1dTrilinear(0u, unitPos, gridDim, 0.0).rgb, material.a); }
                case 4: {
                    material = vec4<f32>(texture3dFrom1dTrilinear(0u, unitPos, gridDim, instance * f32(uObject.uVertexCount)).rgb, material.a);
                }
                default: { material = vec4<f32>(read_tColor(i32(instance)).rgb, material.a); }
            }
        #endif

        #if defined(dOverpaint)
            #if defined(dOverpaintType_texture)
                let overpaint = read_tOverpaint(dataIndex);
            #elif defined(dOverpaintType_vertexInstance)
                let overpaint = texture3dFrom1dTrilinear(1u, unitPos, gridDim, instance * f32(uObject.uVertexCount));
            #else
                let overpaint = vec4<f32>(0.0);
            #endif
            material = vec4<f32>(mix(material.rgb, overpaint.rgb, overpaint.a), material.a);
        #endif

        var color = vec4<f32>(material.rgb, 0.0);
        #if !defined(dIgnoreLight)
            if (material.a >= 0.01) {
                #if defined(dPackedGroup)
                    // gradient by central differences
                    let gradient = vec3<f32>(
                        textureVal(unitPos - dx).a - textureVal(unitPos + dx).a,
                        textureVal(unitPos - dy).a - textureVal(unitPos + dy).a,
                        textureVal(unitPos - dz).a - textureVal(unitPos + dz).a
                    );
                #else
                    let gradient = cell.xyz * 2.0 - 1.0;
                #endif
                // a zero gradient (byte normals filtered halfway between 127 and 128) faces the camera, not NaN
                let normal = select(vec3<f32>(0.0, 0.0, 1.0), -normalize(normalMatrix * gradient), dot(gradient, gradient) > 1e-12);
                color = applyLight(material, normal, viewPosition, metalness, roughness, emissive);
                color = vec4<f32>(color.rgb * uGlobal.uExposure, color.a);
            }
        #endif

        // the transfer function's opacity is per cell: a step of 1 / stepsPerCell cells (\`uTransferScale\`) is as opaque
        // as that fraction of a cell, exactly (the linear \`a * uTransferScale\` is only close for small opacities)
        color.a = 1.0 - pow(1.0 - clamp(material.a * uMaterial.uAlpha, 0.0, 1.0), uMaterial.uTransferScale);

        #if defined(dColorMarker)
            var marker = uObject.uMarker;
            if (uObject.uMarker == -1.0) {
                marker = read_tMarker(dataIndex).r;
                marker = floor(marker * 255.0 + 0.5);
            }
            if (marker > 0.0) {
                let priority = uGlobal.uMarkerPriority;
                if ((priority == 1 && marker != 2.0) || (priority != 1 && marker == 1.0)) {
                    color = vec4<f32>(mix(color.rgb, uGlobal.uHighlightColor, uGlobal.uHighlightStrength), color.a);
                } else {
                    color = vec4<f32>(mix(color.rgb, uGlobal.uSelectColor, uGlobal.uSelectStrength), color.a);
                }
            } else if (uGlobal.uMarkerAverage > 0.0) {
                color = vec4<f32>(mix(color.rgb, uGlobal.uDimColor, uGlobal.uDimStrength), color.a);
            }
        #endif

        preFogAlphaBlended = (1.0 - preFogAlphaBlended) * color.a + preFogAlphaBlended;
        let fragmentDepth = calcDepth(viewPosition);
        if (fog) {
            let viewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
            let fogFactor = smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(viewZ));
            let fogAlpha = (1.0 - fogFactor) * color.a;
            if (!transparentBackground) {
                if (color.a < 1.0) {
                    // transparent objects are blended with background color
                    color.a = fogAlpha;
                } else {
                    // mix opaque objects with background color
                    color = vec4<f32>(mix(color.rgb, uGlobal.uFogColor, fogFactor), color.a);
                }
            } else {
                // pre-multiplied alpha expected for transparent background
                color = vec4<f32>(color.rgb * fogAlpha, fogAlpha);
            }
        }

        var src = color;
        if (!transparentBackground || !fog) {
            // done by the fog otherwise
            src = vec4<f32>(src.rgb * src.a, src.a);
        }
        dst = (1.0 - dst.a) * src + dst; // standard blending

        // break if the color is opaque enough
        if (dst.a > 0.95) {
            break;
        }

        t += 1.0;
        #endif
    }

    return dst;
}

@fragment
fn fragmentMain(in: VertexOutput, @builtin(front_facing) frontFacing: bool) -> FragmentOutput {
    var out: FragmentOutput;
    if (frontFacing) {
        discard;
    }

    let cameraPosition = uGlobal.uCameraPosition;
    let isOrtho = uGlobal.uIsOrtho;
    let rayDir = mix(normalize(in.vModelPosition - cameraPosition), uGlobal.uCameraDir, isOrtho);
    let step = rayDir * uObject.uStepScale * uGlobal.uModelScale;

    let boundingSphere = in.vBoundingSphere;
    let boundingSphereNear = distance(boundingSphere.xyz, cameraPosition) - boundingSphere.w;
    let d = max(uGlobal.uNear, boundingSphereNear) - mix(0.0, distance(in.vModelPosition, cameraPosition), isOrtho);
    let start = mix(cameraPosition, in.vModelPosition, isOrtho) + (d * rayDir);
    let transform = mat4x4<f32>(in.vTransform0, in.vTransform1, in.vTransform2, in.vTransform3);
    let normalMatrix = mat3x3<f32>(in.vNormalMatrix0, in.vNormalMatrix1, in.vNormalMatrix2);
    #if defined(dRenderVariant_marking) || defined(dRenderVariant_pick)
        _ = raymarch(start, step, in.position.xy, in.vInstance, transform, normalMatrix);
        if (surfaceDepth < 0.0) {
            discard;
        }
        let fragmentDepth = surfaceDepth;
        out.fragDepth = fragmentDepth;
        #if defined(dRenderVariant_marking)
            let marker = surfaceMarker;
            let material = vec4<f32>(0.0, 0.0, 0.0, 1.0);
            ${marking_output_wgsl}
        #else
            // too deep in the fog to be picked, as other objects (\`pick_output_wgsl\`)
            let viewZ = depthToViewZ(uGlobal.uIsOrtho, fragmentDepth, uGlobal.uNear, uGlobal.uFar);
            if ((1.0 - smoothstep(uGlobal.uFogNear, uGlobal.uFogFar, abs(viewZ))) * uMaterial.uAlpha < 0.1) {
                discard;
            }
            out.object = vec4<f32>(packIntToRGB(f32(uObject.uObjectId)), 1.0);
            out.instance = vec4<f32>(packIntToRGB(in.vInstance), 1.0);
            out.group = vec4<f32>(packIntToRGB(surfaceGroup), 1.0);
            out.depth = packDepthToRGBA(fragmentDepth);
        #endif
        return out;
    #else
        var color = raymarch(start, step, in.position.xy, in.vInstance, transform, normalMatrix);

        let fragmentDepth = calcDepth((uGlobal.uView * vec4<f32>(start, 1.0)).xyz);
        let preFogAlpha = clamp(preFogAlphaBlended, 0.0, 1.0);
        let interior = false;
        ${wboit_output_wgsl}
        #if defined(dRenderVariant_colorDpoit)
            out.dpoitBack = vec4<f32>(0.0);
            out.dpoitDepth = vec4<f32>(0.0);
        #endif
        out.color = color;
        return out;
    #endif
}
`;
