/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canComputeGrid3dOnGPU, createGrid3dComputeRenderable } from '@molstar/graphics/gl/compute/grid3d';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { TextureSpec, type UnboxedValues, UniformSpec } from '@molstar/graphics/gl/renderable/schema';
import { RuntimeContext } from '@molstar/core/task';
import { ValueCell } from '@molstar/core/util';
import { arrayMin } from '@molstar/core/util/array';
import type { AlphaOrbital, Basis, CubeGridInfo } from '@molstar/alpha-orbitals-extension/data-model';
import { normalizeBasicOrder, type SphericalBasisOrder } from '@molstar/alpha-orbitals-extension/spherical-functions';
import { MAIN, UTILS } from './shader.frag.js';
import { MAIN_WGSL, UTILS_WGSL } from './shader.wgsl.js';

const Schema = {
  tCenters: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
  tInfo: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
  tCoeff: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
  tAlpha: TextureSpec('image-float32', 'r', 'float', 'nearest'),
  uNCenters: UniformSpec('i'),
  uNAlpha: UniformSpec('i'),
  uNCoeff: UniformSpec('i'),
  uMaxCoeffs: UniformSpec('i'),
};

const Orbitals = createGrid3dComputeRenderable({
  schema: Schema,
  loopBounds: ['uNCenters', 'uMaxCoeffs'],
  mainCode: MAIN,
  utilCode: UTILS,
  returnCode: 'v',
  wgsl: { utilCode: UTILS_WGSL, mainCode: MAIN_WGSL, returnCode: 'v' },
  values(params: { grid: CubeGridInfo; orbital: AlphaOrbital }) {
    return createTextureData(params.grid, params.orbital);
  },
});

const Density = createGrid3dComputeRenderable({
  schema: {
    ...Schema,
    uOccupancy: UniformSpec('f'),
  },
  loopBounds: ['uNCenters', 'uMaxCoeffs'],
  mainCode: MAIN,
  utilCode: UTILS,
  returnCode: 'current + uOccupancy * v * v',
  wgsl: { utilCode: UTILS_WGSL, mainCode: MAIN_WGSL, returnCode: 'current + uParams.uOccupancy * v * v' },
  values(params: { grid: CubeGridInfo; orbitals: AlphaOrbital[] }) {
    return {
      ...createTextureData(params.grid, params.orbitals[0]),
      uOccupancy: 0,
    };
  },
  cumulative: {
    states(params: { grid: CubeGridInfo; orbitals: AlphaOrbital[] }) {
      return params.orbitals.filter((o) => o.occupancy !== 0);
    },
    update({ grid }, state: AlphaOrbital, values) {
      const alpha = getNormalizedAlpha(grid.params.basis, state.alpha, grid.params.sphericalOrder);
      ValueCell.updateIfChanged(values.uOccupancy, state.occupancy);
      ValueCell.update(values.tAlpha, { width: alpha.length, height: 1, array: alpha });
    },
  },
});

/** WebGL, or compute kernels (WebGPU) */
export function canComputeAlphaOrbitalsOnGPU(graphics?: GraphicsContext): graphics is GraphicsContext {
  return canComputeGrid3dOnGPU(graphics, Orbitals) && canComputeGrid3dOnGPU(graphics, Density);
}

export function gpuComputeAlphaOrbitalsGridValues(
  ctx: RuntimeContext,
  graphics: GraphicsContext,
  grid: CubeGridInfo,
  orbital: AlphaOrbital,
) {
  return Orbitals(ctx, graphics, grid, { grid, orbital });
}

export function gpuComputeAlphaOrbitalsDensityGridValues(
  ctx: RuntimeContext,
  graphics: GraphicsContext,
  grid: CubeGridInfo,
  orbitals: AlphaOrbital[],
) {
  return Density(ctx, graphics, grid, { grid, orbitals });
}

function getNormalizedAlpha(basis: Basis, alphaOrbitals: number[], sphericalOrder: SphericalBasisOrder) {
  const alpha = new Float32Array(alphaOrbitals.length);

  let aO = 0;
  for (const atom of basis.atoms) {
    for (const shell of atom.shells) {
      for (const L of shell.angularMomentum) {
        const a0 = normalizeBasicOrder(L, alphaOrbitals.slice(aO, aO + 2 * L + 1), sphericalOrder);
        for (let i = 0; i < a0.length; i++) alpha[aO + i] = a0[i];
        aO += 2 * L + 1;
      }
    }
  }

  return alpha;
}

function createTextureData(grid: CubeGridInfo, orbital: AlphaOrbital): UnboxedValues<typeof Schema> {
  const { basis, sphericalOrder, cutoffThreshold } = grid.params;

  let centerCount = 0;
  let baseCount = 0;
  let coeffCount = 0;
  for (const atom of basis.atoms) {
    for (const shell of atom.shells) {
      for (const L of shell.angularMomentum) {
        if (L > 4) {
          // TODO: will L > 4 be required? Would need to precompute more functions in that case.
          throw new Error('Angular momentum L > 4 not supported.');
        }

        centerCount++;
        baseCount += 2 * L + 1;
        coeffCount += shell.exponents.length;
      }
    }
  }

  const centers = new Float32Array(4 * centerCount);
  // L, alpha_offset, coeff_offset_start, coeff_offset_end
  const info = new Float32Array(4 * centerCount);
  const alpha = new Float32Array(baseCount);
  const coeff = new Float32Array(4 * coeffCount);

  let maxCoeffs = 0;
  let cO = 0,
    aO = 0,
    coeffO = 0;
  for (const atom of basis.atoms) {
    for (const shell of atom.shells) {
      let amIndex = 0;
      for (const L of shell.angularMomentum) {
        const a0 = normalizeBasicOrder(L, orbital.alpha.slice(aO, aO + 2 * L + 1), sphericalOrder);

        const cutoffRadius =
          cutoffThreshold > 0 ? Math.sqrt(-Math.log(cutoffThreshold) / arrayMin(shell.exponents)) : 10000;

        centers[4 * cO + 0] = atom.center[0];
        centers[4 * cO + 1] = atom.center[1];
        centers[4 * cO + 2] = atom.center[2];
        centers[4 * cO + 3] = cutoffRadius * cutoffRadius;

        info[4 * cO + 0] = L;
        info[4 * cO + 1] = aO;
        info[4 * cO + 2] = coeffO;
        info[4 * cO + 3] = coeffO + shell.exponents.length;

        for (let i = 0; i < a0.length; i++) alpha[aO + i] = a0[i];

        const c0 = shell.coefficients[amIndex++];
        for (let i = 0; i < shell.exponents.length; i++) {
          coeff[4 * (coeffO + i) + 0] = c0[i];
          coeff[4 * (coeffO + i) + 1] = shell.exponents[i];
        }

        if (c0.length > maxCoeffs) {
          maxCoeffs = c0.length;
        }

        cO++;
        aO += 2 * L + 1;
        coeffO += shell.exponents.length;
      }
    }
  }

  return {
    uNCenters: centerCount,
    uNAlpha: baseCount,
    uNCoeff: coeffCount,
    uMaxCoeffs: maxCoeffs,
    tCenters: { width: centerCount, height: 1, array: centers },
    tInfo: { width: centerCount, height: 1, array: info },
    tCoeff: { width: coeffCount, height: 1, array: coeff },
    tAlpha: { width: baseCount, height: 1, array: alpha },
  };
}
