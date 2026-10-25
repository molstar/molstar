/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** Volume values, the cell grid and the marching cubes case of a cell, shared by the marching cubes kernels */
export const mc_common_wgsl = /* wgsl */ `
// cube corners and the corners of the 12 edges, as the WebGL isosurface
const Corners = array<vec3<i32>, 8>(
    vec3<i32>(0, 0, 0), vec3<i32>(1, 0, 0), vec3<i32>(1, 1, 0), vec3<i32>(0, 1, 0),
    vec3<i32>(0, 0, 1), vec3<i32>(1, 0, 1), vec3<i32>(1, 1, 1), vec3<i32>(0, 1, 1)
);
const EdgeCorners = array<vec2<u32>, 12>(
    vec2<u32>(0u, 1u), vec2<u32>(1u, 2u), vec2<u32>(2u, 3u), vec2<u32>(3u, 0u),
    vec2<u32>(4u, 5u), vec2<u32>(5u, 6u), vec2<u32>(6u, 7u), vec2<u32>(7u, 4u),
    vec2<u32>(0u, 4u), vec2<u32>(1u, 5u), vec2<u32>(2u, 6u), vec2<u32>(3u, 7u)
);

/** the value at grid position \`p\`, clamped to the grid (buffers: wrapped into the data, at most one point larger) */
fn value(p: vec3<i32>) -> f32 {
    let q = clamp(p, vec3<i32>(0), uParams.uGridDim - 1);
    #if defined(dInput_buffer)
        let d = uParams.uDataDim;
        return volume[dot(select(q, q - d, q >= d), uParams.uStrides)];
    #else
        let v = textureLoad(tVolume, q, 0);
        #if defined(dChannel_red)
            return v.r;
        #else
            return v.a;
        #endif
    #endif
}

/** the case of the cell at \`c\`: a bit per corner at or below the iso value, 0 for the inside of the surface */
fn cubeCase(c: vec3<i32>) -> u32 {
    var k = 0u;
    for (var i = 0u; i < 8u; i++) {
        if (value(c + Corners[i]) <= uParams.uIsoValue) {
            k |= 1u << i;
        }
    }
    return select(k, 0u, k == 255u);
}

fn cellIndex(c: vec3<u32>, cells: vec3<u32>) -> u32 {
    return c.x + c.y * cells.x + c.z * cells.x * cells.y;
}
`;

/** The gradient normal and the group packing of the emit kernel */
export const mc_vertex_wgsl = /* wgsl */ `
fn gradientNormal(p: vec3<i32>) -> vec3<f32> {
    let g = vec3<f32>(
        value(p - vec3<i32>(1, 0, 0)) - value(p + vec3<i32>(1, 0, 0)),
        value(p - vec3<i32>(0, 1, 0)) - value(p + vec3<i32>(0, 1, 0)),
        value(p - vec3<i32>(0, 0, 1)) - value(p + vec3<i32>(0, 0, 1))
    );
    // a flat neighbourhood has no gradient (\`normalize\` of a zero vector is NaN): no normal, as the callers check
    if (all(g == vec3<f32>(0.0))) {
        return vec3<f32>(0.0);
    }
    return -normalize(g);
}

/** an integer (group id) packed into rgb, as \`packIntToRGB\` */
fn packGroup(g: u32) -> vec3<f32> {
    if (g >= 16777215u) {
        return vec3<f32>(1.0);
    }
    let v = g + 1u;
    return vec3<f32>(f32((v >> 16u) & 255u), f32((v >> 8u) & 255u), f32(v & 255u)) / 255.0;
}
`;
