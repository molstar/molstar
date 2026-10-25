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
 * Per sampled vertex and instance, the color weighted by \`2 - distance\` added to the grid cells within a distance of 2,
 * as the point sprites (7 pixels wide) of the WebGL accumulation; sums are fixed-point (rgba and weight per cell), for
 * the slab of \`uSlab.y\` slices from slice \`uSlab.x\`
 */
export const color_accumulate_wgsl = /* wgsl */ `
${fixed_point_wgsl}

fn unpackRGBToInt(rgb: vec3<f32>) -> u32 {
    let v = vec3<u32>(floor(rgb * 255.0 + 0.5));
    return (v.r << 16u) + (v.g << 8u) + v.b - 1u;
}

${iso_vertex_wgsl('isoTable', 'isoVertices')}

fn texelOf(index: u32, width: u32) -> vec2<u32> {
    return vec2<u32>(index % width, index / width);
}

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
    let group = packedGroup - 1u;
    let position = (transforms[instance] * vec4<f32>(localPosition, 1.0)).xyz;

    #if defined(dColorType_groupInstance)
        let colorIndex = u32(instances[instance]) * u32(uParams.uGroupCount) + group;
    #else
        let colorIndex = group;
    #endif
    let color = textureLoad(tColor, texelOf(colorIndex, textureDimensions(tColor).x), 0);

    let gridDim = uParams.uGridDim;
    let p = (position - uParams.uBboxMin) / uParams.uResolution;
    // pixel center of the point sprite in a slice
    let sprite = (position.xy - uParams.uBboxMin.xy) / uParams.uBboxSize.xy * vec2<f32>(gridDim.xy);
    let slabMin = vec3<i32>(0, 0, uParams.uSlab.x);
    let slabMax = vec3<i32>(gridDim.xy, uParams.uSlab.x + uParams.uSlab.y) - 1;
    let lo = max(vec3<i32>(ceil(p - 2.0)), slabMin);
    let hi = min(vec3<i32>(floor(p + 2.0)), min(gridDim - 1, slabMax));
    for (var z = lo.z; z <= hi.z; z++) {
        for (var y = lo.y; y <= hi.y; y++) {
            for (var x = lo.x; x <= hi.x; x++) {
                let c = vec3<f32>(f32(x), f32(y), f32(z));
                if (any(abs(c.xy + 0.5 - sprite) > vec2<f32>(3.5))) {
                    continue;
                }
                let dist = distance(c, p);
                if (dist > 2.0) {
                    continue;
                }
                let f = 2.0 - dist;
                let cell = u32(x + y * gridDim.x + (z - slabMin.z) * gridDim.x * gridDim.y) * 5u;
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
 * The weighted colors of the slab \`uSlab\` divided by the weights into the 2d atlas of the grid, slices in rows of
 * \`uTexCols\`; texels without weight are zero, those past the last slice too (with the first slab)
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
    let slabStart = u32(uParams.uSlab.x);
    let inSlab = column < texCols && z >= slabStart && z < min(gridDim.z, slabStart + u32(uParams.uSlab.y));
    if (!inSlab && !(slabStart == 0u && z >= gridDim.z)) {
        return;
    }
    var color = vec4<f32>(0.0);
    if (inSlab) {
        let slice = z - slabStart;
        let cell = ((id.x % gridDim.x) + (id.y % gridDim.y) * gridDim.x + slice * gridDim.x * gridDim.y) * 5u;
        let weight = f32(sums[cell + 4u]);
        if (weight > 0.0) {
            color = vec4<f32>(f32(sums[cell]), f32(sums[cell + 1u]), f32(sums[cell + 2u]), f32(sums[cell + 3u])) / weight;
        }
    }
    textureStore(outColor, id.xy, color);
}
`;
