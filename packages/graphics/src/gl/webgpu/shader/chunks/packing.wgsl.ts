/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Packing of depth and unit interval values into 8-bit channels, as in GLSL */
export const packing_wgsl = /* wgsl */ `
const PackUpscale = 256.0 / 255.0; // fraction -> 0..1 (including 1)
const UnpackDownscale = 255.0 / 256.0; // 0..1 -> fraction (excluding 1)
const PackFactors = vec3<f32>(256.0 * 256.0 * 256.0, 256.0 * 256.0, 256.0);
const UnpackFactors = UnpackDownscale / vec4<f32>(PackFactors, 1.0);
const ShiftRight8 = 1.0 / 256.0;

fn packDepthToRGBA(v: f32) -> vec4<f32> {
    let r = vec4<f32>(fract(v * PackFactors), v);
    // tidy overflow
    return vec4<f32>(r.x, r.yzw - r.xyz * ShiftRight8) * PackUpscale;
}

fn unpackRGBAToDepth(v: vec4<f32>) -> f32 {
    return dot(v, UnpackFactors);
}

fn packDepthWithAlphaToRGBA(depth: f32, alpha: f32) -> vec4<f32> {
    let r = vec3<f32>(fract(depth * PackFactors.yz), depth);
    // tidy overflow
    return vec4<f32>(vec3<f32>(r.x, r.yz - r.xy * ShiftRight8) * PackUpscale, alpha);
}

fn unpackRGBAToDepthWithAlpha(v: vec4<f32>) -> vec2<f32> {
    return vec2<f32>(dot(v.xyz, UnpackFactors.yzw), v.w);
}

fn packUnitIntervalToRG(v: f32) -> vec2<f32> {
    let x = fract(v * 256.0);
    return vec2<f32>(x, v - x * (1.0 / 256.0)) * (256.0 / 255.0);
}

fn unpackRGToUnitInterval(enc: vec2<f32>) -> f32 {
    return dot(enc, vec2<f32>(255.0 / (256.0 * 256.0), 255.0 / 256.0));
}

/** occlusion factor from a packed SSAO texel, values close to 0.0 are treated as errors (no occlusion) */
fn unpackSsao(v: vec4<f32>) -> f32 {
    let ssao = unpackRGToUnitInterval(v.xy);
    return select(1.0, ssao, ssao > 0.001 && ssao <= 0.999);
}

/** GL window coordinates (x, y bottom-up in [0, 1], depth in [0, 1]) to view space with a GL inverse projection */
fn screenSpaceToViewSpace(ssPos: vec3<f32>, invProjection: mat4x4<f32>) -> vec3<f32> {
    let p = invProjection * vec4<f32>(ssPos * 2.0 - 1.0, 1.0);
    return p.xyz / p.w;
}
`;
