/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * A smoothed layer's value per vertex and instance of a mesh (`vertexInstance`), sampled from its grid as the mesh
 * shaders sample a grid (`texture3dFrom2dLinear`, the texture's linear filter at the cells' corners), written into `outData` from element
 * `uOutOffset`: element `instance * uVertexCount + vertex`, instances not reached (`reached`) zero. `rgba8` one `u32`
 * an element (`dItemSize_4`); `r8` four elements a `u32` (low byte first), a thread a word, the bytes of the words at
 * the range's ends outside it kept
 */
export const layer_sample_vertices_wgsl = /* wgsl */ `
fn sampleGrid(pos: vec3<f32>) -> vec4<f32> {
    let gridDim = uParams.uGridDim;
    let texDim = uParams.uGridTexDim;
    let zSlice0 = floor(pos.z * gridDim.z);
    // the tile of a slice in integers, as \`texture3dFrom2dLinear\`; the values at the cells' corners, a texel read at
    // its center half a cell on in x and y
    let width = i32(texDim.x);
    let start0 = i32(zSlice0) * i32(gridDim.x);
    let start1 = start0 + i32(gridDim.x);
    let column0 = f32(start0 % width);
    let row0 = f32(start0 / width);
    let column1 = f32(start1 % width);
    let row1 = f32(start1 / width);
    let offset = pos.xy * gridDim.xy + 0.5;
    let color0 = textureSampleLevel(tGrid, sLinear, (vec2<f32>(column0, row0 * gridDim.y) + offset) / texDim, 0.0);
    let color1 = textureSampleLevel(tGrid, sLinear, (vec2<f32>(column1, row1 * gridDim.y) + offset) / texDim, 0.0);
    return mix(color0, color1, abs(pos.z * gridDim.z - zSlice0));
}

fn valueAt(element: u32) -> vec4<f32> {
    let vertexCount = u32(uParams.uVertexCount);
    let instance = element / vertexCount;
    if (reached[instance] == 0u) {
        return vec4<f32>(0.0);
    }
    let o = (u32(uParams.uVertexOffset) + element % vertexCount) * 3u;
    let local = vec3<f32>(positions[o], positions[o + 1u], positions[o + 2u]);
    let position = (transforms[u32(uParams.uInstanceOffset) + instance] * vec4<f32>(local, 1.0)).xyz;
    let t = uParams.uGridTransform;
    return sampleGrid((t.w * (position - t.xyz)) / uParams.uGridDim);
}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    let index = id.x + id.y * groups.x * WorkgroupSizeX;
    let count = u32(uParams.uVertexCount) * u32(uParams.uInstanceCount);
    let outOffset = u32(uParams.uOutOffset);
    #if defined(dItemSize_4)
        if (index >= count) {
            return;
        }
        outData[outOffset + index] = pack4x8unorm(valueAt(index));
    #else
        let first = outOffset / 4u;
        if (index >= (outOffset % 4u + count + 3u) / 4u) {
            return;
        }
        let word = first + index;
        let old = outData[word];
        var value = 0u;
        for (var b = 0u; b < 4u; b++) {
            let element = i32(word * 4u + b) - i32(outOffset);
            var byte = (old >> (8u * b)) & 255u;
            if (element >= 0 && u32(element) < count) {
                byte = u32(round(clamp(valueAt(u32(element)).r, 0.0, 1.0) * 255.0));
            }
            value = value | (byte << (8u * b));
        }
        outData[word] = value;
    #endif
}
`;
