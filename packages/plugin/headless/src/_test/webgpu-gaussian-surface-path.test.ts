/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  GaussianSurfaceMeshParams,
  getGaussianSurfaceGpuFallbackReason,
} from '@molstar/graphics/repr/structure/visual/gaussian-surface-mesh';
import type { Structure } from '@molstar/model/model/structure';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

/** `count` carbon atoms on a grid in a ball at protein density (11 Å³ per atom), as PDB */
function globulePdb(count: number) {
  const spacing = Math.cbrt(11);
  const n = Math.ceil(Math.cbrt(count)) + 2;
  const points: number[][] = [];
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) points.push([i - n / 2, j - n / 2, k - n / 2].map((v) => v * spacing));
    }
  }
  points.sort((a, b) => Math.hypot(...a) - Math.hypot(...b));
  const f = (v: number) => v.toFixed(3).padStart(8);
  const lines = points.slice(0, count).map((p, i) => {
    const serial = String((i % 99999) + 1).padStart(5);
    const residue = String((i % 9999) + 1).padStart(4);
    return `ATOM  ${serial}  CA  ALA A${residue}    ${f(p[0])}${f(p[1])}${f(p[2])}  1.00  0.00           C`;
  });
  return `${lines.join('\n')}\nEND\n`;
}

describe.skipIf(!WebGPU.module)('Gaussian surfaces on WebGPU: CPU or GPU', () => {
  it('computes small surfaces on the CPU, others on the GPU at any resolution', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const structure = async (count: number) => {
        const data = await plugin.builders.data.rawData({ data: globulePdb(count), label: `${count}` });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        const model = await plugin.builders.structure.createModel(trajectory);
        return (await plugin.builders.structure.createStructure(model)).data as Structure;
      };
      const reason = (s: Structure, resolution: number) =>
        getGaussianSurfaceGpuFallbackReason(
          s,
          { ...PD.getDefaultValues(GaussianSurfaceMeshParams), resolution },
          graphics,
        );

      const small = await structure(300);
      expect(small.elementCount).toBe(300);
      // under 100k density cells and 5000 atoms: faster on the CPU, at resolution 1 and above
      expect(reason(small, 1)).toMatch(/small/);
      expect(reason(small, 2)).toMatch(/small/);
      // fine grids are not small
      expect(reason(small, 0.5)).toBeUndefined();

      const large = await structure(6000);
      // many atoms: the GPU also at coarse resolutions (the CPU cutoff of the WebGL path doesn't apply)
      expect(reason(large, 1)).toBeUndefined();
      expect(reason(large, 3)).toBeUndefined();
      // as asked for
      const off = { ...PD.getDefaultValues(GaussianSurfaceMeshParams), tryUseGpu: false };
      expect(getGaussianSurfaceGpuFallbackReason(large, off, graphics)).toMatch(/tryUseGpu/);
    } finally {
      plugin.dispose();
    }
  });
});
