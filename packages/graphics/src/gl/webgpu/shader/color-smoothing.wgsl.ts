/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { iso_vertex_wgsl } from './marching-cubes-chunks.wgsl.js';

/** Fixed-point scale of the accumulated weights and weighted colors */
const fixed_point_wgsl = /* wgsl */ `
const Scale = 16384.0;
`;

/**
 * A sample's color weighted by \`2 - distance\` added to the grid cells within a distance of 2, of the region
 * (\`uRegionMin\` to \`uRegionMax\`, exclusive) and the slab of \`uSlab.y\` slices from slice \`uSlab.x\` (within the
 * region); sums are fixed-point (rgba and weight per cell), x fastest over the region's extent in x and y
 */
const color_splat_wgsl = /* wgsl */ `
${fixed_point_wgsl}

fn texelOf(index: u32, width: u32) -> vec2<u32> {
    return vec2<u32>(index % width, index / width);
}

fn colorOf(instance: u32, group: u32) -> vec4<f32> {
    #if defined(dColorType_groupInstance)
        let colorIndex = u32(instances[instance]) * u32(uParams.uGroupCount) + group;
    #else
        let colorIndex = group;
    #endif
    return textureLoad(tColor, texelOf(colorIndex, textureDimensions(tColor).x), 0);
}

/** \`clampToSprite\`: the cells of the point sprite (7 pixels wide, centered at \`sprite\`) only, as WebGL draws them */
fn splat(position: vec3<f32>, color: vec4<f32>, clampToSprite: bool, sprite: vec2<f32>) {
    let p = (position - uParams.uBboxMin) / uParams.uResolution;
    let regionMin = vec3<i32>(uParams.uRegionMin.xy, uParams.uSlab.x);
    let regionMax = min(vec3<i32>(uParams.uRegionMax.xy, uParams.uSlab.x + uParams.uSlab.y), uParams.uRegionMax);
    let lo = max(vec3<i32>(ceil(p - 2.0)), regionMin);
    let hi = min(vec3<i32>(floor(p + 2.0)), regionMax - 1);
    let extent = uParams.uRegionMax.xy - uParams.uRegionMin.xy;
    for (var z = lo.z; z <= hi.z; z++) {
        for (var y = lo.y; y <= hi.y; y++) {
            for (var x = lo.x; x <= hi.x; x++) {
                let c = vec3<f32>(f32(x), f32(y), f32(z));
                if (clampToSprite && any(abs(c.xy + 0.5 - sprite) > vec2<f32>(3.5))) {
                    continue;
                }
                let dist = distance(c, p);
                if (dist > 2.0) {
                    continue;
                }
                let f = 2.0 - dist;
                let r = vec3<i32>(x, y, z) - regionMin;
                let cell = u32(r.x + r.y * extent.x + r.z * extent.x * extent.y) * 5u;
                let weighted = vec4<u32>(round(color * f * Scale));
                atomicAdd(&sums[cell], weighted.r);
                atomicAdd(&sums[cell + 1u], weighted.g);
                atomicAdd(&sums[cell + 2u], weighted.b);
                atomicAdd(&sums[cell + 3u], weighted.a);
                atomicAdd(&sums[cell + 4u], u32(round(f * Scale)));
            }
        }
    }
}
`;

/**
 * Per sampled vertex of the isosurface's arena (texture meshes) and instance: its color splatted (\`splat\`), within
 * the point sprite of the WebGL accumulation
 */
export const color_accumulate_wgsl = /* wgsl */ `
${color_splat_wgsl}

${iso_vertex_wgsl('isoTable', 'isoVertices')}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    let index = id.x + id.y * groups.x * WorkgroupSizeX;
    let sampleCount = u32(uParams.uSampleCount);
    if (index >= sampleCount * u32(uParams.uInstanceCount)) {
        return;
    }
    let instance = index / sampleCount;
    let vertex = (index % sampleCount) * u32(uParams.uStride);
    // in pages of the vertex arena (see \`computeChunkedIsosurface\`), drawn up to the page table's fourth word;
    // \`uVertexCount\` is an upper bound
    if (vertex >= min(u32(uParams.uVertexCount), isoDrawn())) {
        return;
    }
    let slot = isoSlot(vertex);
    let localPosition = isoPosition(slot);
    let packedGroup = isoPackedGroup(slot);
    // a vertex no end of whose edge has a group: no color (the group would wrap to 0xFFFFFFFF)
    if (packedGroup == 0u) {
        return;
    }
    let position = (transforms[instance] * vec4<f32>(localPosition, 1.0)).xyz;
    // pixel center of the point sprite in a slice
    let sprite = (position.xy - uParams.uBboxMin.xy) / uParams.uBboxSize.xy * vec2<f32>(uParams.uGridDim.xy);
    splat(position, colorOf(instance, packedGroup - 1u), true, sprite);
}
`;

/**
 * Per sampled vertex of a mesh (\`positions\` xyz and \`groups\` from \`uVertexOffset\`, every \`uStride\`-th) and
 * instance: its color splatted (\`splat\`), as the CPU splats a mesh's (\`calcMeshColorSmoothing\`)
 */
export const color_accumulate_mesh_wgsl = /* wgsl */ `
${color_splat_wgsl}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    let index = id.x + id.y * groups.x * WorkgroupSizeX;
    let sampleCount = u32(uParams.uSampleCount);
    if (index >= sampleCount * u32(uParams.uInstanceCount)) {
        return;
    }
    let instance = index / sampleCount;
    let vertex = (index % sampleCount) * u32(uParams.uStride);
    if (vertex >= u32(uParams.uVertexCount)) {
        return;
    }
    let v = u32(uParams.uVertexOffset) + vertex;
    let localPosition = vec3<f32>(positions[v * 3u], positions[v * 3u + 1u], positions[v * 3u + 2u]);
    let position = (transforms[instance] * vec4<f32>(localPosition, 1.0)).xyz;
    splat(position, colorOf(instance, vertexGroups[v]), false, vec2<f32>(0.0));
}
`;

/**
 * The weighted colors of the slab \`uSlab\` of the region (\`uRegionMin\` to \`uRegionMax\`) divided by the weights
 * into the 2d atlas of the grid, slices in rows of \`uTexCols\`; texels without weight are zero, and with the region's
 * first slab those outside the region and past the last slice
 */
export const color_normalize_wgsl = /* wgsl */ `
@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let size = textureDimensions(outColor);
    if (any(id.xy >= size)) {
        return;
    }
    let gridDim = vec3<u32>(uParams.uGridDim);
    let texCols = u32(uParams.uTexCols);
    let column = id.x / gridDim.x;
    let z = (id.y / gridDim.y) * texCols + column;
    let x = id.x % gridDim.x;
    let y = id.y % gridDim.y;
    let regionMin = vec3<u32>(uParams.uRegionMin);
    let regionMax = vec3<u32>(uParams.uRegionMax);
    let slabStart = u32(uParams.uSlab.x);
    let inRegion = column < texCols && z < gridDim.z && all(vec3<u32>(x, y, z) >= regionMin) &&
        all(vec3<u32>(x, y, z) < regionMax);
    let inSlab = inRegion && z >= slabStart && z < slabStart + u32(uParams.uSlab.y);
    if (!inSlab && !(slabStart == regionMin.z && !inRegion)) {
        return;
    }
    var color = vec4<f32>(0.0);
    if (inSlab) {
        let extent = regionMax.xy - regionMin.xy;
        let cell = ((x - regionMin.x) + (y - regionMin.y) * extent.x + (z - slabStart) * extent.x * extent.y) * 5u;
        let weight = f32(sums[cell + 4u]);
        if (weight > 0.0) {
            color = vec4<f32>(f32(sums[cell]), f32(sums[cell + 1u]), f32(sums[cell + 2u]), f32(sums[cell + 3u])) / weight;
        }
    }
    textureStore(outColor, id.xy, color);
}
`;
