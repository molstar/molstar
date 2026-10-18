/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Box3D } from '@molstar/core/math/geometry';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import {
  type ComputeContext,
  KernelCode,
  type KernelRequest,
  StorageBufferSpec,
  StorageTextureSpec,
} from '../compute.js';
import type { GraphicsContext } from '../context.js';
import { UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';
import { getDispatchSize } from './prefix-sum.js';
import { getScratchBuffer } from './scratch.js';

/**
 * The Gaussian density of atoms into a 3d texture with compute kernels, as the splatting into its slices on WebGL: per
 * atom the cells of its quad in each slice near the atom, the density summed (fixed-point atomics) and the nearness of
 * the nearest atom (`gd-accumulate`), then the group of the nearest atom (`gd-group`), written as rgb (packed group)
 * and alpha (density) into the texture (`gd-write`). Grids whose cells exceed the storage buffer limit are processed in
 * slabs of slices.
 */

export const GdAccumulateKernel = KernelCode('gd-accumulate', [64]);
export const GdGroupKernel = KernelCode('gd-group', [64]);
export const GdWriteKernel = KernelCode('gd-write', [4, 4, 4]);

const GridUniforms = {
  uGridDim: UniformSpec('iv3'),
  /** first slice and slice count of the slab */
  uSlab: UniformSpec('iv2'),
};

const AtomInput = {
  ...GridUniforms,
  uAtomCount: UniformSpec('i'),
  /** cells per side of the largest quad */
  uFootprint: UniformSpec('i'),
  uBboxMin: UniformSpec('v3'),
  uBboxSize: UniformSpec('v3'),
  uResolution: UniformSpec('f'),
  uRadiusFactorInv: UniformSpec('f'),
  /** slices above and below an atom with a contribution */
  uZRange: UniformSpec('f'),
  /** position and radius */
  atoms: StorageBufferSpec('vec4<f32>', 'read'),
};

export const GdAccumulateSchema = {
  ...AtomInput,
  uAlpha: UniformSpec('f'),
  density: StorageBufferSpec('atomic<u32>', 'read-write'),
  nearness: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const GdGroupSchema = {
  ...AtomInput,
  nearness: StorageBufferSpec('u32', 'read'),
  winner: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const GdWriteSchema = {
  ...GridUniforms,
  density: StorageBufferSpec('u32', 'read'),
  winner: StorageBufferSpec('u32', 'read'),
  groups: StorageBufferSpec('f32', 'read'),
  outDensity: StorageTextureSpec('rgba', 'fp16', '3d', 'write'),
};

/** Per-dimension dispatch limit of the WebGPU core */
const MaxWorkgroups = 65535;

/** Contributions to the density below this fraction of an atom's peak are left out */
const DensityCutoff = 1e-4;

export function isComputeGaussianDensitySupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  const { compute } = graphics;
  return (
    !!compute &&
    compute.isKernelSupported(GdAccumulateKernel) &&
    compute.isKernelSupported(GdGroupKernel) &&
    compute.isKernelSupported(GdWriteKernel)
  );
}

export interface GaussianDensityAtoms {
  count: number;
  /** xyz per atom */
  positions: Float32Array;
  radii: Float32Array;
  groups: Float32Array;
  maxRadius: number;
}

export interface GaussianDensityGrid {
  box: Box3D;
  dim: Vec3;
  resolution: number;
  /** `alpha` of the Gaussian */
  smoothness: number;
  radiusFactor: number;
}

/** The kernels of `computeGaussianDensity`, see `ComputeContext.prepareKernels` */
export const GaussianDensityKernels: readonly KernelRequest[] = [
  { code: GdAccumulateKernel, schema: GdAccumulateSchema },
  { code: GdGroupKernel, schema: GdGroupSchema },
  { code: GdWriteKernel, schema: GdWriteSchema },
];

/**
 * Writes the density of `atoms` on `grid` into the rgba fp16 3d `texture`, defined to the grid size; with `growTexture`
 * only grown (a texture shared by grids of varying size, read by grid coordinates: the cells beyond the grid are left
 * as they were)
 */
export function computeGaussianDensity(
  compute: ComputeContext,
  atoms: GaussianDensityAtoms,
  grid: GaussianDensityGrid,
  texture: Texture,
  options?: { growTexture?: boolean },
) {
  const { dim, box, resolution, smoothness, radiusFactor } = grid;
  const [dx, dy, dz] = dim;
  if (options?.growTexture) {
    const [w, h, d] = [texture.getWidth(), texture.getHeight(), texture.getDepth()];
    if (w < dx || h < dy || d < dz) texture.define(Math.max(w, dx), Math.max(h, dy), Math.max(d, dz));
  } else {
    texture.define(dx, dy, dz);
  }
  if (atoms.count === 0) return;

  const atomData = new Float32Array(atoms.count * 4);
  for (let i = 0; i < atoms.count; ++i) {
    atomData[i * 4] = atoms.positions[i * 3];
    atomData[i * 4 + 1] = atoms.positions[i * 3 + 1];
    atomData[i * 4 + 2] = atoms.positions[i * 3 + 2];
    atomData[i * 4 + 3] = atoms.radii[i];
  }
  const atomBuffer = compute.createBuffer(atomData.byteLength, 'gaussian density atoms');
  atomBuffer.write(atomData);
  const groupBuffer = compute.createBuffer(atoms.groups.byteLength, 'gaussian density groups');
  groupBuffer.write(atoms.groups);

  const sliceSize = dx * dy * 4;
  const slabDepth = Math.min(dz, Math.floor(compute.maxStorageBufferSize / sliceSize));
  if (slabDepth < 1) throw new Error('a gaussian density slice exceeds the storage buffer limit');
  const slabSize = slabDepth * sliceSize;
  // kept between calls (see `getScratchBuffer`), cleared per slab
  const density = getScratchBuffer(compute, 'gaussian density', slabSize);
  const nearness = getScratchBuffer(compute, 'gaussian density nearness', slabSize);
  const winner = getScratchBuffer(compute, 'gaussian density winner', slabSize);

  const { maxRadius } = atoms;
  const footprint = Math.ceil((maxRadius * 3) / resolution + resolution) + 1;
  // the density down to the cutoff, the nearest atom and its group up to the radius factor (with the group tolerance)
  const densityRange = maxRadius * Math.sqrt(-Math.log(DensityCutoff) / smoothness);
  const zRange = Math.max(densityRange, radiusFactor * (1 + resolution * 0.05)) / resolution;

  const gridDim = ValueCell.create(Vec3.clone(dim));
  const slab = ValueCell.create(Vec2.create(0, slabDepth));
  const atomValues = {
    uGridDim: gridDim,
    uSlab: slab,
    uAtomCount: ValueCell.create(atoms.count),
    uFootprint: ValueCell.create(footprint),
    uBboxMin: ValueCell.create(Vec3.clone(box.min)),
    uBboxSize: ValueCell.create(Vec3.sub(Vec3(), box.max, box.min)),
    uResolution: ValueCell.create(resolution),
    uRadiusFactorInv: ValueCell.create(1 / radiusFactor),
    uZRange: ValueCell.create(zRange),
    atoms: ValueCell.create(atomBuffer),
  };
  const accumulateKernel = compute.createKernel(GdAccumulateKernel, GdAccumulateSchema, {
    ...atomValues,
    uAlpha: ValueCell.create(smoothness),
    density: ValueCell.create(density),
    nearness: ValueCell.create(nearness),
  });
  const groupKernel = compute.createKernel(GdGroupKernel, GdGroupSchema, {
    ...atomValues,
    nearness: ValueCell.create(nearness),
    winner: ValueCell.create(winner),
  });
  const writeKernel = compute.createKernel(GdWriteKernel, GdWriteSchema, {
    uGridDim: gridDim,
    uSlab: slab,
    density: ValueCell.create(density),
    winner: ValueCell.create(winner),
    groups: ValueCell.create(groupBuffer),
    outDensity: ValueCell.create(texture),
  });

  const { x, y } = getDispatchSize(Math.ceil((atoms.count * footprint * footprint) / 64), MaxWorkgroups);
  for (let z = 0; z < dz; z += slabDepth) {
    // the clears end the open compute pass
    density.clear(0, slabSize);
    nearness.clear(0, slabSize);
    winner.clear(0, slabSize);
    const depth = Math.min(slabDepth, dz - z);
    ValueCell.update(slab, Vec2.set(slab.ref.value, z, depth));
    compute.beginComputePass({ label: 'gaussian density' });
    accumulateKernel.dispatch(x, y);
    groupKernel.dispatch(x, y);
    writeKernel.dispatch(Math.ceil(dx / 4), Math.ceil(dy / 4), Math.ceil(depth / 4));
    compute.endComputePass();
  }

  accumulateKernel.destroy();
  groupKernel.destroy();
  writeKernel.destroy();
  atomBuffer.destroy();
  groupBuffer.destroy();
}
