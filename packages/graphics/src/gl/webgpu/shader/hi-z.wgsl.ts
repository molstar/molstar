/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { quad_vert_wgsl } from './chunks/quad.wgsl.js';

/**
 * A level of the hierarchical depth buffer: the min (the farthest, depth is reversed) of the texels of the previous
 * level (or the depth texture) the texel covers: 2x2 of a level, up to 3x3 of the depth texture (level 0 maps the
 * viewport onto a power-of-two square, e.g. 1280 pixels onto 1024 texels: a texel's span can touch three pixels)
 */
export const hi_z_wgsl = /* wgsl */ `
${quad_vert_wgsl}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let size = vec2<i32>(textureDimensions(tPreviousLevel));
    let s = vec2<f32>(size);
    // the texel's span in the previous level's texels: those it touches (at most three a side)
    let start = ((in.position.xy - 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let end = ((in.position.xy + 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let lo = clamp(vec2<i32>(floor(start)), vec2<i32>(0), size - 1);
    let hi = clamp(vec2<i32>(ceil(end)) - 1, lo, min(lo + 2, size - 1));

    var depth = 1.0;
    for (var y = lo.y; y <= hi.y; y++) {
        for (var x = lo.x; x <= hi.x; x++) {
            depth = min(depth, textureLoad(tPreviousLevel, vec2<i32>(x, y), 0).r);
        }
    }
    return vec4<f32>(depth);
}
`;

/**
 * The first level primed from an occluder grid (`HiZOccluderSchema`): the depth texture's farthest (as `hi_z_wgsl`), the
 * nearer (larger, depth reversed) where the texel's center ray through the grid reaches `uOpticalDepth`: at the far end
 * of the step it is reached in (a density a lower bound, so a depth not nearer than the content's statistically)
 */
export const hi_z_occluder_wgsl = /* wgsl */ `
${quad_vert_wgsl}

/** the optical density of grid cell \`c\` (its slices along z side by side in \`tGrid\`) */
fn gridDensity(c: vec3<i32>) -> f32 {
    let n = vec3<i32>(uObject.uGridDims);
    let columns = i32(uObject.uGridColumns);
    return textureLoad(tGrid, vec2<i32>(c.x + (c.z % columns) * n.x, c.y + (c.z / columns) * n.y), 0).r;
}

/** the depth where the ray through the texel at \`position\` reaches the threshold in the grid, 0 (the farthest) where not */
fn gridDepth(position: vec2<f32>) -> f32 {
    // its center in the drawing buffer's uv, then the viewport's normalized device coordinates (y up)
    let uv = position * uObject.uInvSize + uObject.uRect.xy;
    let vr = uObject.uViewportRect;
    let ndc = vec2<f32>((uv.x - vr.x) / vr.z * 2.0 - 1.0, 1.0 - (uv.y - vr.y) / vr.w * 2.0);
    // the near plane's point (depth reversed: 1) and one farther along the ray
    let a = uObject.uInvViewProjection * vec4<f32>(ndc, 1.0, 1.0);
    let b = uObject.uInvViewProjection * vec4<f32>(ndc, 0.5, 1.0);
    let origin = a.xyz / a.w;
    let dir = normalize(b.xyz / b.w - origin);
    let lo = uObject.uGridMin;
    let hi = uObject.uGridMax;
    // the slabs' latest entry and earliest exit (an axis the ray runs along: within or not)
    var tn = 0.0;
    var tf = 1e30;
    for (var k = 0; k < 3; k++) {
        if (abs(dir[k]) < 1e-12) {
            if (origin[k] < lo[k] || origin[k] > hi[k]) {
                return 0.0;
            }
            continue;
        }
        let t1 = (lo[k] - origin[k]) / dir[k];
        let t2 = (hi[k] - origin[k]) / dir[k];
        tn = max(tn, min(t1, t2));
        tf = min(tf, max(t1, t2));
    }
    if (!(tf > tn)) {
        return 0.0;
    }
    let n = uObject.uGridDims;
    let cell = (hi - lo) / n;
    let step = min(cell.x, min(cell.y, cell.z)) * 0.5;
    let steps = min(i32(ceil((tf - tn) / step)), 1024);
    var sum = 0.0;
    for (var k = 0; k < steps; k++) {
        let t = tn + f32(k) * step;
        let p = origin + dir * (t + 0.5 * step);
        let c = clamp(vec3<i32>(floor((p - lo) / cell)), vec3<i32>(0), vec3<i32>(n) - 1);
        sum += gridDensity(c) * step;
        if (sum >= uObject.uOpticalDepth) {
            let clip = uObject.uViewProjection * vec4<f32>(origin + dir * (t + step), 1.0);
            return clamp(clip.z / clip.w, 0.0, 1.0);
        }
    }
    return 0.0;
}

@fragment
fn fragmentMain(in: VertexOutput) -> @location(0) vec4<f32> {
    let size = vec2<i32>(textureDimensions(tPreviousLevel));
    let s = vec2<f32>(size);
    let start = ((in.position.xy - 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let end = ((in.position.xy + 0.5) * uObject.uInvSize + uObject.uRect.xy) * s;
    let lo = clamp(vec2<i32>(floor(start)), vec2<i32>(0), size - 1);
    let hi = clamp(vec2<i32>(ceil(end)) - 1, lo, min(lo + 2, size - 1));
    var depth = 1.0;
    for (var y = lo.y; y <= hi.y; y++) {
        for (var x = lo.x; x <= hi.x; x++) {
            depth = min(depth, textureLoad(tPreviousLevel, vec2<i32>(x, y), 0).r);
        }
    }
    return vec4<f32>(max(depth, gridDepth(in.position.xy)));
}
`;
