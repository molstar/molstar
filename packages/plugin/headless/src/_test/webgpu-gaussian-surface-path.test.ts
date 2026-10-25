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
import type { Structure, Unit } from '@molstar/model/model/structure';
import type { StorageBuffer } from '@molstar/graphics/gl/compute';
import { readChunkedVertices } from '@molstar/graphics/gl/compute/marching-cubes/chunks';
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

/** a globule (`globulePdb`) of each of `counts` atoms, chains `A`, `B`, … on a line `distance` Å apart, as PDB */
function chainsPdb(counts: number[], distance: number) {
  const lines: string[] = [];
  counts.forEach((count, c) => {
    const chain = String.fromCharCode(65 + c);
    for (const line of globulePdb(count).split('\n').slice(0, count)) {
      const x = (Number.parseFloat(line.slice(30, 38)) + c * distance).toFixed(3).padStart(8);
      lines.push(`${line.slice(0, 21)}${chain}${line.slice(22, 30)}${x}${line.slice(38)}`);
    }
  });
  return `${lines.join('\n')}\nEND\n`;
}

/** `n`³ carbon atoms on a lattice `spacing` Å apart, slightly jittered, as PDB */
function latticePdb(n: number, spacing: number) {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const lines: string[] = [];
  for (let i = 0; i < n * n * n; ++i) {
    const p = [i % n, Math.floor(i / n) % n, Math.floor(i / (n * n))].map(
      (v, k) => v * spacing + 0.7 * Math.sin(i * (k + 1.3)),
    );
    const serial = String((i % 99999) + 1).padStart(5);
    const residue = String((i % 9999) + 1).padStart(4);
    lines.push(`ATOM  ${serial}  CA  ALA A${residue}    ${f(p[0])}${f(p[1])}${f(p[2])}  1.00  0.00           C`);
  }
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
      const reason = (s: Structure | Unit, resolution: number) =>
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

      // a unit's surface by the unit's own box and atoms, not the structure's: chains far apart fit the GPU where the
      // structure's box would not, and small chains of a large structure stay on the CPU
      const data = await plugin.builders.data.rawData({
        data: chainsPdb([3000, 3000, 3000, 3000], 1000),
        label: 'chains',
      });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const chains = (await plugin.builders.structure.createStructure(model)).data as Structure;
      expect(chains.units.length).toBe(4);
      expect(reason(chains, 0.5)).toMatch(/along an axis/);
      for (const unit of chains.units) expect(reason(unit, 0.5)).toBeUndefined();
      expect(reason(chains, 2)).toBeUndefined();
      for (const unit of chains.units) expect(reason(unit, 2)).toMatch(/small/);
    } finally {
      plugin.dispose();
    }
  });

  // a cell's group is its first corner's; at coarse resolutions that corner may lie beyond the range the density records
  // the nearest atom in (no group: drawn black), as for the separate atoms of a trace-only surface of a large structure
  it('gives every vertex of a coarse GPU surface a group', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: latticePdb(18, 12), label: 'lattice' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const canvas = plugin.canvas3d!;
      for (const resolution of [3, 8]) {
        const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
          type: 'gaussian-surface',
          typeParams: { quality: 'custom', resolution, tryUseGpu: true, radiusOffset: 2, smoothness: 1 },
        });
        canvas.commit(true);
        await canvas.graphics.resources.pending.whenReady();
        const [surface] = canvas.getRenderObjects().filter((o) => o.type === 'texture-mesh');
        const values = surface.values as typeof surface.values & {
          sIsoTable: { ref: { value: StorageBuffer } };
          sIsoVertices: { ref: { value: StorageBuffer } };
          uGroupCount: { ref: { value: number } };
        };
        const table = values.sIsoTable.ref.value;
        const count = new Uint32Array(await table.readAsync(0, 16))[3];
        const chunks = { table, vertices: values.sIsoVertices.ref.value };
        const { groups, vertexCount } = await readChunkedVertices(canvas.graphics.compute!, chunks, count);
        expect(vertexCount).toBeGreaterThan(0);
        let missing = 0;
        for (let i = 0; i < vertexCount; ++i) {
          // packed as group + 1, 0 for none
          const packed = groups[i * 4] * 65536 + groups[i * 4 + 1] * 256 + groups[i * 4 + 2];
          if (packed === 0 || packed > values.uGroupCount.ref.value) missing += 1;
        }
        expect(missing).toBe(0);
        await plugin.build().delete(repr!).commit();
      }
    } finally {
      plugin.dispose();
    }
  });

  it('chooses the GPU or the CPU per unit, and again when the props change', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: chainsPdb([6000, 300, 300], 100), label: 'chains' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const canvas = plugin.canvas3d!;
      const types = async () => {
        canvas.commit(true);
        await canvas.graphics.resources.pending.whenReady();
        return canvas
          .getRenderObjects()
          .map((o) => o.type)
          .sort();
      };
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'gaussian-surface',
        typeParams: { quality: 'custom', resolution: 1, tryUseGpu: true },
      });
      // the large chain on the GPU, the small ones on the CPU
      expect(await types()).toEqual(['mesh', 'mesh', 'texture-mesh']);
      // on a fine grid the small chains are not small: recreated on the GPU
      await plugin
        .build()
        .to(repr!)
        .update((old) => {
          old.type.params = { ...old.type.params, resolution: 0.5 };
        })
        .commit();
      expect(await types()).toEqual(['texture-mesh', 'texture-mesh', 'texture-mesh']);
    } finally {
      plugin.dispose();
    }
  });
});
