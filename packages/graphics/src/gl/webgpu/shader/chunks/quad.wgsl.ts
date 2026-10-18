/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Vertex stage of full-screen (or scaled) quads of the compute programs, `QuadSchema`, and `fragmentUv`: the quad
 * programs do their math in window space with a top-left origin, as the neutral API and WebGPU's textures
 */
export const quad_vert_wgsl = /* wgsl */ `
struct VertexOutput {
    @builtin(position) position: vec4<f32>,
};

@vertex
fn vertexMain(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    let scale = uObject.uQuadScale;
    out.position = vec4<f32>(in.aPosition * scale - vec2<f32>(1.0) + scale, 0.0, 1.0);
    return out;
}

/** uv of the fragment in a target of \`size\` pixels, top-left origin */
fn fragmentUv(position: vec4<f32>, size: vec2<f32>) -> vec2<f32> {
    return position.xy / size;
}
`;

/**
 * `read_<name>(uv)`: texel of `name` at a uv with a top-left origin (rendered content), filtered with
 * `dLinear_<name>` (the bound texture is linear and filterable), else the nearest texel
 */
export function texture_read_wgsl(name: string) {
  return /* wgsl */ `
fn read_${name}(uv: vec2<f32>) -> vec4<f32> {
    #if defined(dLinear_${name})
        return textureSampleLevel(${name}, sLinear, uv, 0.0);
    #else
        let size = vec2<i32>(textureDimensions(${name}));
        return textureLoad(${name}, clamp(vec2<i32>(floor(uv * vec2<f32>(size))), vec2<i32>(0), size - 1), 0);
    #endif
}
`;
}

/**
 * `sample_<name>(uv)`: texel of `name` at a uv in upload order, for textures loaded from images or data (both backends
 * store their rows in the same order, unlike render targets); filtered with `dLinear_<name>`, else the nearest texel
 */
export function texture_sample_wgsl(name: string) {
  return /* wgsl */ `
fn sample_${name}(uv: vec2<f32>) -> vec4<f32> {
    #if defined(dLinear_${name})
        return textureSampleLevel(${name}, sLinear, uv, 0.0);
    #else
        let size = vec2<i32>(textureDimensions(${name}));
        return textureLoad(${name}, clamp(vec2<i32>(floor(uv * vec2<f32>(size))), vec2<i32>(0), size - 1), 0);
    #endif
}
`;
}
