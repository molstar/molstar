/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import { pack_int_wgsl } from './chunks/common.wgsl.js';

/** Fixed-point scale of the summed density */
const gd_density_scale_wgsl = /* wgsl */ `
const DensityScale = 4194304.0;
`;

/**
 * Per atom and column of its footprint (the quad of the splatting into slices on WebGL, \`uFootprint\` × \`uFootprint\` cells at most),
 * the cells of the slab within \`uZRange\` slices of the atom
 */
const gd_footprint_wgsl = /* wgsl */ `
${gd_density_scale_wgsl}

struct Footprint {
    valid: bool,
    atom: u32,
    /** the atom in grid cells */
    center: vec3<f32>,
    radius: f32,
    x: i32,
    y: i32,
    zMin: i32,
    zMax: i32,
};

fn getFootprint(id: vec3<u32>, workgroupCount: vec3<u32>) -> Footprint {
    var f: Footprint;
    f.valid = false;
    let index = id.x + id.y * workgroupCount.x * WorkgroupSizeX;
    let footprint = u32(uParams.uFootprint);
    let cellsPerAtom = footprint * footprint;
    if (index >= u32(uParams.uAtomCount) * cellsPerAtom) {
        return f;
    }
    f.atom = index / cellsPerAtom;
    let cell = index % cellsPerAtom;
    let atom = atoms[f.atom];
    f.radius = atom.w;
    let resolution = uParams.uResolution;
    f.center = (atom.xyz - uParams.uBboxMin) / resolution;

    // pixel centers within the quad: its size and position on the slice as the WebGL splatting computes them
    let size = ceil(((f.radius * 3.0) / resolution) + resolution);
    let quadCenter = (atom.xy - uParams.uBboxMin.xy) / uParams.uBboxSize.xy * vec2<f32>(uParams.uGridDim.xy);
    let first = vec2<i32>(floor(quadCenter - size * 0.5));
    let xy = first + vec2<i32>(i32(cell % footprint), i32(cell / footprint));
    if (any(abs(vec2<f32>(xy) + 0.5 - quadCenter) >= vec2<f32>(size * 0.5))) {
        return f;
    }
    if (any(xy < vec2<i32>(0)) || any(xy >= uParams.uGridDim.xy)) {
        return f;
    }
    f.x = xy.x;
    f.y = xy.y;
    let zRange = uParams.uZRange;
    f.zMin = max(i32(ceil(f.center.z - zRange)), uParams.uSlab.x);
    f.zMax = min(i32(floor(f.center.z + zRange)), min(uParams.uGridDim.z, uParams.uSlab.x + uParams.uSlab.y) - 1);
    f.valid = f.zMin <= f.zMax;
    return f;
}

fn cellIndex(x: i32, y: i32, z: i32) -> u32 {
    let dim = uParams.uGridDim;
    return u32(x + y * dim.x + (z - uParams.uSlab.x) * dim.x * dim.y);
}

/** distance of the cell to the atom in Å, as the WebGL splatting computes it */
fn cellDistance(f: Footprint, z: i32) -> f32 {
    return distance(vec3<f32>(f32(f.x), f32(f.y), f32(z)), f.center) * uParams.uResolution;
}

/** \`1 - distance / radiusFactor\` in an 8-bit channel, as the minimum distance texture of the WebGL splatting stores it */
fn quantizedNearness(dist: f32) -> u32 {
    return u32(round(clamp(1.0 - dist * uParams.uRadiusFactorInv, 0.0, 1.0) * 255.0));
}
`;

/** The summed density (fixed point) and the nearness of the nearest atom per cell */
export const gd_accumulate_wgsl = /* wgsl */ `
${gd_footprint_wgsl}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) workgroupCount: vec3<u32>) {
    let f = getFootprint(id, workgroupCount);
    if (!f.valid) {
        return;
    }
    let radiusSqInv = 1.0 / (f.radius * f.radius);
    for (var z = f.zMin; z <= f.zMax; z++) {
        let dist = cellDistance(f, z);
        let cell = cellIndex(f.x, f.y, z);
        let d = exp(-uParams.uAlpha * ((dist * dist) * radiusSqInv)) * uParams.uRadiusFactorInv;
        atomicAdd(&density[cell], u32(round(d * DensityScale)));
        let q = quantizedNearness(dist);
        if (q > 0u) {
            atomicMax(&nearness[cell], q);
        }
    }
}
`;

/** Per cell the last atom (in atom order, the draw order of the WebGL splatting) within the distance of the nearest one */
export const gd_group_wgsl = /* wgsl */ `
${gd_footprint_wgsl}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) workgroupCount: vec3<u32>) {
    let f = getFootprint(id, workgroupCount);
    if (!f.valid) {
        return;
    }
    for (var z = f.zMin; z <= f.zMax; z++) {
        let dist = cellDistance(f, z);
        let cell = cellIndex(f.x, f.y, z);
        let minDistance = 1.0 - f32(nearness[cell]) / 255.0;
        if (dist * uParams.uRadiusFactorInv <= minDistance + uParams.uResolution * 0.05) {
            atomicMax(&winner[cell], f.atom + 1u);
        }
    }
}
`;

/** The density in alpha and the group of the nearest atom packed in rgb, per cell of the slab */
export const gd_write_wgsl = /* wgsl */ `
${pack_int_wgsl}
${gd_density_scale_wgsl}

@compute @workgroup_size(WorkgroupSizeX, WorkgroupSizeY, WorkgroupSizeZ)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let dim = vec3<u32>(uParams.uGridDim);
    let slabStart = u32(uParams.uSlab.x);
    let z = id.z + slabStart;
    if (any(id.xy >= dim.xy) || id.z >= u32(uParams.uSlab.y) || z >= dim.z) {
        return;
    }
    let cell = id.x + id.y * dim.x + id.z * dim.x * dim.y;
    let w = winner[cell];
    var group = vec3<f32>(0.0);
    if (w > 0u) {
        group = packIntToRGB(groups[w - 1u]);
    }
    textureStore(outDensity, vec3<u32>(id.xy, z), vec4<f32>(group, f32(density[cell]) / DensityScale));
}
`;
