/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { depth_wgsl } from './chunks/depth.wgsl.js';
import { packing_wgsl } from './chunks/packing.wgsl.js';
import { quad_vert_wgsl, texture_read_wgsl } from './chunks/quad.wgsl.js';
import type { DefineLookup } from '../../renderable/schema.js';

/** The blended transparency of the luminosity seed, read with `dBlendTransparency` only */
const luminosity_transparent_wgsl = /* wgsl */ `
${texture_read_wgsl('tColorTransparent')}
${texture_read_wgsl('tDepthTransparent')}

fn getColorTransparent(coords: vec2<f32>) -> vec4<f32> {
    return read_tColorTransparent(coords);
}

fn getDepthTransparent(coords: vec2<f32>) -> f32 {
    return unpackRGBAToDepthWithAlpha(read_tDepthTransparent(coords)).x;
}
`;

/** Without blended transparency: nothing in front of the opaque scene (the far plane, depth reversed) */
const luminosity_opaque_wgsl = /* wgsl */ `
fn getColorTransparent(coords: vec2<f32>) -> vec4<f32> {
    return vec4<f32>(0.0);
}

fn getDepthTransparent(coords: vec2<f32>) -> f32 {
    return 0.0;
}
`;

/**
 * Bloom seed: bright parts of the composited scene (luminosity) or the emissive prepass; the transparent textures are
 * read with `dBlendTransparency` only (`WgslProgram.source`)
 */
export const bloom_luminosity_wgsl = (d: DefineLookup) => /* wgsl */ `
${depth_wgsl}
${packing_wgsl}
${quad_vert_wgsl}
${texture_read_wgsl('tColorOpaque')}
${texture_read_wgsl('tEmissive')}
${texture_read_wgsl('tDepthOpaque')}
${d('dBlendTransparency') ? luminosity_transparent_wgsl : luminosity_opaque_wgsl}

fn getDepthOpaque(coords: vec2<f32>) -> f32 {
    // depth textures are always available on WebGPU
    return read_tDepthOpaque(coords).r;
}

fn isBackground(depth: f32) -> bool {
    // the far plane (depth is reversed), packed or raw
    return depth == 0.0;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, 1.0 / uObject.uTexSizeInv);
    var opaqueTexel = read_tColorOpaque(coords);
    let transparentTexel = getColorTransparent(coords);
    // fog illumination's un-fogged opaque seed by the opaque depth, matching standard's per-layer fog
    if (uObject.uOpaqueFogged == 0u) {
        let viewZ = depthToViewZ(uObject.uIsOrtho, getDepthOpaque(coords), uObject.uNear, uObject.uFar);
        opaqueTexel = vec4<f32>(opaqueTexel.rgb * (1.0 - smoothstep(uObject.uFogNear, uObject.uFogFar, abs(viewZ))), opaqueTexel.a);
    }
    // PMA OVER composite, matches postprocessing's final transparency blend
    let texel = transparentTexel + opaqueTexel * (1.0 - transparentTexel.a);
    let depth = max(getDepthOpaque(coords), getDepthTransparent(coords));

    if (isBackground(depth)) {
        return vec4<f32>(0.0);
    }

    #if defined(dMode_luminosity)
        let outputColor = vec4<f32>(uObject.uDefaultColor, uObject.uDefaultOpacity);
        let luma = vec3<f32>(0.299, 0.587, 0.114);
        let v = dot(texel.xyz, luma);
        let alpha = smoothstep(uObject.uLuminosityThreshold, uObject.uLuminosityThreshold + uObject.uSmoothWidth, v);
        return mix(outputColor, texel, alpha);
    #else
        // the prepass already holds the fogged emissive color, free of lighting/exposure
        return read_tEmissive(coords);
    #endif
}
`;

/** Separable gaussian blur */
export const bloom_blur_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tInput')}

override oKernelRadius: i32 = 3;

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let texSizeInv = uObject.uTexSizeInv;
    let coords = fragmentUv(in.position, 1.0 / texSizeInv);
    var weightSum = uObject.uGaussianCoefficients[0].x;
    var diffuseSum = read_tInput(coords) * weightSum;

    for (var i = 1; i < oKernelRadius; i++) {
        let x = f32(i);
        let w = uObject.uGaussianCoefficients[i].x;
        let offset = uObject.uDirection * texSizeInv * x;
        let sample1 = read_tInput(coords + offset);
        let sample2 = read_tInput(coords - offset);
        diffuseSum += (sample1 + sample2) * w;
        weightSum += 2.0 * w;
    }

    return diffuseSum / weightSum;
}
`;

/** Combines the blur levels */
export const bloom_composite_wgsl = /* wgsl */ `
${quad_vert_wgsl}
${texture_read_wgsl('tBlur1')}
${texture_read_wgsl('tBlur2')}
${texture_read_wgsl('tBlur3')}
${texture_read_wgsl('tBlur4')}
${texture_read_wgsl('tBlur5')}

fn lerpBloomFactor(factor: f32) -> f32 {
    let mirrorFactor = 1.2 - factor;
    return mix(factor, mirrorFactor, uObject.uBloomRadius);
}

fn level(i: i32) -> vec4<f32> {
    return lerpBloomFactor(uObject.uBloomFactors[i].x) * vec4<f32>(uObject.uBloomTints[i].xyz, 1.0);
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let coords = fragmentUv(in.position, 1.0 / uObject.uTexSizeInv);
    return uObject.uBloomStrength * (
        level(0) * read_tBlur1(coords) +
        level(1) * read_tBlur2(coords) +
        level(2) * read_tBlur3(coords) +
        level(3) * read_tBlur4(coords) +
        level(4) * read_tBlur5(coords)
    );
}
`;
