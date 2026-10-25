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

/** one texel of a packed SSAO target for \`getSsaoUpsampled\`, its occlusion and weight */
fn ssaoUpsampleTexel(tex: texture_2d<f32>, size: vec2<i32>, texel: vec2<i32>, bilinear: f32, viewZ: f32, near: f32, far: f32) -> vec2<f32> {
    let v = textureLoad(tex, clamp(texel, vec2<i32>(0), size - 1), 0);
    let texelViewZ = unpackRGToUnitInterval(v.zw) * (far - near) - far;
    let weight = bilinear / (0.001 + abs(texelViewZ - viewZ));
    return vec2<f32>(unpackSsao(v) * weight, weight);
}

/**
 * occlusion of a packed SSAO target of a lower resolution at \`coords\` (top-left origin), for a pixel at view z
 * \`viewZ\`: the four nearest texels weighted bilinearly and by their depth, so that occlusion does not bleed across
 * depth edges
 */
fn getSsaoUpsampled(tex: texture_2d<f32>, coords: vec2<f32>, viewZ: f32, near: f32, far: f32) -> f32 {
    let size = vec2<i32>(textureDimensions(tex));
    let p = coords * vec2<f32>(size) - 0.5;
    let i = vec2<i32>(floor(p));
    let f = p - floor(p);
    var s = ssaoUpsampleTexel(tex, size, i, (1.0 - f.x) * (1.0 - f.y), viewZ, near, far);
    s += ssaoUpsampleTexel(tex, size, i + vec2<i32>(1, 0), f.x * (1.0 - f.y), viewZ, near, far);
    s += ssaoUpsampleTexel(tex, size, i + vec2<i32>(0, 1), (1.0 - f.x) * f.y, viewZ, near, far);
    s += ssaoUpsampleTexel(tex, size, i + vec2<i32>(1, 1), f.x * f.y, viewZ, near, far);
    return select(1.0, s.x / s.y, s.y > 0.0);
}

/** Window coordinates (x, y in [0, 1] with a top-left origin, depth in [0, 1], reversed) to view space */
fn screenSpaceToViewSpace(ssPos: vec3<f32>, invProjection: mat4x4<f32>) -> vec3<f32> {
    let p = invProjection * vec4<f32>(ssPos.x * 2.0 - 1.0, 1.0 - ssPos.y * 2.0, ssPos.z, 1.0);
    return p.xyz / p.w;
}

// the unit vector \`n\` octahedral-encoded in 2 × 12 bits over 8-bit rgb, \`a\` (unit interval) in alpha; rgb 0 is not
// produced (a zero vector, nothing drawn, see \`unpackRGBAToNormal\`), as in GLSL
fn packNormalToRGBA(n: vec3<f32>, a: f32) -> vec4<f32> {
    var p = n.xy / (abs(n.x) + abs(n.y) + abs(n.z));
    if (n.z < 0.0) {
        p = (1.0 - abs(p.yx)) * select(vec2<f32>(-1.0), vec2<f32>(1.0), p >= vec2<f32>(0.0));
    }
    var q = floor(clamp(p * 0.5 + 0.5, vec2<f32>(0.0), vec2<f32>(1.0)) * 4095.0 + 0.5);
    if (q.x == 0.0 && q.y == 0.0) {
        q.y = 1.0;
    }
    let xHigh = floor(q.x / 16.0);
    let yHigh = floor(q.y / 256.0);
    return vec4<f32>(vec3<f32>(xHigh, (q.x - xHigh * 16.0) * 16.0 + yHigh, q.y - yHigh * 256.0) / 255.0, a);
}

// the unit vector of \`packNormalToRGBA\`, a zero vector for rgb 0 (cleared)
fn unpackRGBAToNormal(v: vec4<f32>) -> vec3<f32> {
    let c = floor(v.rgb * 255.0 + 0.5);
    if (all(c == vec3<f32>(0.0))) {
        return vec3<f32>(0.0);
    }
    let gHigh = floor(c.y / 16.0);
    let p = vec2<f32>(c.x * 16.0 + gHigh, (c.y - gHigh * 16.0) * 256.0 + c.z) / 4095.0 * 2.0 - 1.0;
    var n = vec3<f32>(p, 1.0 - abs(p.x) - abs(p.y));
    let t = max(-n.z, 0.0);
    n = vec3<f32>(n.xy + select(vec2<f32>(t), vec2<f32>(-t), n.xy >= vec2<f32>(0.0)), n.z);
    return normalize(n);
}
`;
