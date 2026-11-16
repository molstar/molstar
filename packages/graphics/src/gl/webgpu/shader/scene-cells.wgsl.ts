/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CulledData } from '../../culled-data.js';
import { SceneCellBoxWords } from '../../compute/scene-cells.js';

/** A float's bits ordered as unsigned integers (`atomicMax` on them the float's max), and back */
const ordered_bits_wgsl = /* wgsl */ `
fn orderedBits(x: f32) -> u32 {
    let b = bitcast<u32>(x);
    return select(b | 0x80000000u, ~b, (b & 0x80000000u) != 0u);
}
fn orderedFloat(u: u32) -> f32 {
    return bitcast<f32>(select(~u, u & 0x7fffffffu, (u & 0x80000000u) != 0u));
}
`;

/**
 * Per instance (see `scene-cells.ts`): its sphere (its segment's invariant sphere by its transform, the radius by the
 * largest scale, as the culling's `instanceSphere`), its cell in its segment's class's grid by its center; the cell's
 * box grown by it: the min as the max of the inverted bits (the boxes start zeroed)
 */
export const scene_cell_bin_wgsl = /* wgsl */ `
${ordered_bits_wgsl}
${CulledData.wgsl('data')}
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    // rows of workgroups past the per-dimension limit (\`getDispatchSize\`)
    let i = id.x + id.y * groups.x * WorkgroupSizeX;
    if (i >= u32(uParams.uInstanceCount)) {
        return;
    }
    let segment = u32(culledSegment(i).w);
    let s = segmentData[segment * 3u];
    // the transform's rows (\`affine-rows.ts\`)
    let rows = culledRows(i);
    let r0 = rows[0];
    let r1 = rows[1];
    let r2 = rows[2];
    let center = vec3<f32>(dot(r0.xyz, s.xyz) + r0.w, dot(r1.xyz, s.xyz) + r1.w, dot(r2.xyz, s.xyz) + r2.w);
    let c0 = vec3<f32>(r0.x, r1.x, r2.x);
    let c1 = vec3<f32>(r0.y, r1.y, r2.y);
    let c2 = vec3<f32>(r0.z, r1.z, r2.z);
    let radius = s.w * sqrt(max(dot(c0, c0), max(dot(c1, c1), dot(c2, c2))));

    let k = cellClasses[segment];
    let grid = uParams.uSceneGrids[k];
    let n = vec3<i32>(grid.yzw);
    let g = clamp(vec3<i32>(floor((center - uParams.uSceneMin.xyz) * grid.x)), vec3<i32>(0), n - 1);
    let cell = u32(uParams.uSceneFirsts[k]) + u32(g.x + n.x * (g.y + n.y * g.z));
    instanceCells[i] = vec2<u32>(cell, segment);

    let q = cell * ${SceneCellBoxWords}u;
    atomicMax(&cellBoxes[q], ~orderedBits(center.x));
    atomicMax(&cellBoxes[q + 1u], ~orderedBits(center.y));
    atomicMax(&cellBoxes[q + 2u], ~orderedBits(center.z));
    atomicMax(&cellBoxes[q + 3u], orderedBits(center.x));
    atomicMax(&cellBoxes[q + 4u], orderedBits(center.y));
    atomicMax(&cellBoxes[q + 5u], orderedBits(center.z));
    // radii are not negative: their bits ordered as they are
    atomicMax(&cellBoxes[q + 6u], bitcast<u32>(max(radius, 0.0)));
}
`;

/**
 * Per cell: its sphere from its box (\`Sphere3D.fromBox3D\`) expanded by its largest radius; an empty cell (its max
 * still zero: no float's ordered bits are) a negative radius
 */
export const scene_cell_spheres_wgsl = /* wgsl */ `
${ordered_bits_wgsl}
@compute @workgroup_size(WorkgroupSizeX)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    // rows of workgroups past the per-dimension limit (\`getDispatchSize\`)
    let cell = id.x + id.y * groups.x * WorkgroupSizeX;
    if (cell >= u32(uParams.uCellCount)) {
        return;
    }
    let q = cell * ${SceneCellBoxWords}u;
    if (cellBoxes[q + 3u] == 0u) {
        cellSpheres[cell] = vec4<f32>(0.0, 0.0, 0.0, -1.0);
        return;
    }
    let lo = vec3<f32>(orderedFloat(~cellBoxes[q]), orderedFloat(~cellBoxes[q + 1u]), orderedFloat(~cellBoxes[q + 2u]));
    let hi = vec3<f32>(orderedFloat(cellBoxes[q + 3u]), orderedFloat(cellBoxes[q + 4u]), orderedFloat(cellBoxes[q + 5u]));
    cellSpheres[cell] = vec4<f32>((lo + hi) * 0.5, length((hi - lo) * 0.5) + bitcast<f32>(cellBoxes[q + 6u]));
}
`;
