/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * The 16-bit depth key of each splat of an instance (its view z, `uRow`, within `uDepthMin` and `1 / uDepthScale` of
 * the instance's bounds), counted per key into `histogram`
 */
export const splat_keys_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    // workgroups beyond the per-dimension limit are dispatched in 2d
    let i = (wid.x + wid.y * u32(uParams.uGroupsX)) * WorkgroupSizeX + lid.x;
    if (i >= u32(uParams.uCount)) {
        return;
    }
    let z = dot(uParams.uRow.xyz, positions[i].xyz);
    let key = u32(clamp((z - uParams.uDepthMin) * uParams.uDepthScale, 0.0, 65535.0));
    keys[i] = key;
    atomicAdd(&histogram[key], 1u);
}
`;

/**
 * The centers of SOG's splats read in place (the \`sog\` encoding) from its means: 16 bits an axis within the means'
 * range, out of the log domain, as the splat shaders decode them; a \`vec4\` a splat
 */
export const splat_centers_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let i = (wid.x + wid.y * u32(uParams.uGroupsX)) * WorkgroupSizeX + lid.x;
    if (i >= u32(uParams.uCount)) {
        return;
    }
    let w = textureDimensions(tMeansLow).x;
    let p = vec2<i32>(i32(i % w), i32(i / w));
    let lo = floor(textureLoad(tMeansLow, p, 0).rgb * 255.0 + 0.5);
    let hi = floor(textureLoad(tMeansHigh, p, 0).rgb * 255.0 + 0.5);
    let n = uParams.uMeansMin + (uParams.uMeansMax - uParams.uMeansMin) * ((lo + hi * 256.0) / 65535.0);
    positions[i] = vec4<f32>(sign(n) * (exp(abs(n)) - 1.0), 0.0);
}
`;

/**
 * Each splat of an instance into the instance's block of the order (`uOffset`) at its key's next slot: `offsets` the
 * keys' first slots (the histogram's exclusive prefix sum), advanced as they fill; back to front, equal keys in any order
 */
export const splat_scatter_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(local_invocation_id) lid: vec3<u32>, @builtin(workgroup_id) wid: vec3<u32>) {
    let i = (wid.x + wid.y * u32(uParams.uGroupsX)) * WorkgroupSizeX + lid.x;
    if (i >= u32(uParams.uCount)) {
        return;
    }
    let slot = atomicAdd(&offsets[keys[i]], 1u);
    order[u32(uParams.uOffset) + slot] = i;
}
`;
