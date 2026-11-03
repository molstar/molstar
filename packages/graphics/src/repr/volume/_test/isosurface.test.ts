/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { perVertexLines } from '../../../geo/geometry/_test/item-vertices.js';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { type Grid, Volume } from '@molstar/model/model/volume';
import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import {
  type AtlasIsosurfaceExtraction,
  registerAtlasIsosurfaceExtraction,
} from '../../../gl/compute/marching-cubes/atlas.js';
import type { GraphicsContext } from '../../../gl/context.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  createVolumeIsosurfaceMesh,
  createVolumeIsosurfaceWireframe,
  IsosurfaceMeshParams,
  IsosurfaceVisual,
  type VolumeIsosurfaceProps,
} from '../isosurface.js';

const N = 12;
/** Above every cell id (N³ = 1728): a floodfilled `get` would turn the ids of filled cells into `2 * Iso` */
const Iso = 2000;

/** A shell around a cavity: floodfill 'inside' fills the cavity, 'outside' the outside, each leaving one surface */
function createTestVolume(n = N): Volume {
  const space = Tensor.Space([n, n, n], [2, 1, 0], Float32Array);
  const data = space.create();
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const r = Math.hypot(i - 5.5, j - 5.5, k - 5.5);
        space.set(data, i, j, k, 4000 * (1 - Math.abs(r - 3.5) / 1.5));
      }
    }
  }
  return {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: -10000, max: 4000, mean: 0, sigma: 2000 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

/** The group of a vertex is the cell at one end of its (axis-aligned) edge */
function expectEdgeCellGroups(volume: Volume, positions: ArrayLike<number>, groups: ArrayLike<number>, count: number) {
  const { dataOffset } = volume.grid.cells.space;
  let checked = 0;
  for (let v = 0; v < count; v++) {
    const x = positions[3 * v],
      y = positions[3 * v + 1],
      z = positions[3 * v + 2];
    const lo = dataOffset(Math.floor(x), Math.floor(y), Math.floor(z));
    const hi = dataOffset(Math.ceil(x), Math.ceil(y), Math.ceil(z));
    expect([lo, hi], `vertex ${v} at ${x} ${y} ${z}`).toContain(groups[v]);
    checked++;
  }
  expect(checked).toBeGreaterThan(0);
}

const ctx = { runtime: SyncRuntimeContext } as any;

function props(floodfill: VolumeIsosurfaceProps['floodfill']) {
  return { isoValue: Volume.IsoValue.absolute(Iso), wrap: 'off', floodfill } as VolumeIsosurfaceProps;
}

describe('volume isosurface', () => {
  for (const floodfill of ['off', 'inside', 'outside'] as const) {
    it(`groups vertices by their cells, floodfill ${floodfill}`, async () => {
      const volume = createTestVolume();
      const mesh = await createVolumeIsosurfaceMesh(ctx, volume, 0, undefined as any, props(floodfill));
      expectEdgeCellGroups(volume, mesh.vertexBuffer.ref.value, mesh.groupBuffer.ref.value, mesh.vertexCount);
    });
  }

  for (const floodfill of ['off', 'inside', 'outside'] as const) {
    it(`groups wireframe lines by their cells, floodfill ${floodfill}`, async () => {
      const volume = createTestVolume();
      const lines = await createVolumeIsosurfaceWireframe(ctx, volume, 0, undefined as any, props(floodfill));
      // each line has 4 vertices with its start and end and one group
      const { lineCount } = lines;
      const { start: starts, group: groups } = perVertexLines(lines.lineData.ref.value.array, lineCount);
      const { dataOffset } = volume.grid.cells.space;
      expect(lineCount).toBeGreaterThan(0);
      for (let l = 0; l < lineCount; l++) {
        const x = starts[12 * l],
          y = starts[12 * l + 1],
          z = starts[12 * l + 2];
        const lo = dataOffset(Math.floor(x), Math.floor(y), Math.floor(z));
        const hi = dataOffset(Math.ceil(x), Math.ceil(y), Math.ceil(z));
        expect([lo, hi], `line ${l}`).toContain(groups[4 * l]);
      }
    });
  }

  it('meshes one surface less with a floodfill', async () => {
    const volume = createTestVolume();
    const counts: number[] = [];
    for (const floodfill of ['off', 'inside', 'outside'] as const) {
      counts.push((await createVolumeIsosurfaceMesh(ctx, volume, 0, undefined as any, props(floodfill))).vertexCount);
    }
    // both surfaces without, the outer one with 'inside', the cavity's with 'outside'
    expect(counts[1]).toBeLessThan(counts[0]);
    expect(counts[2]).toBeLessThan(counts[0]);
    expect(counts[1] + counts[2]).toBe(counts[0]);
  });
});

describe('volume isosurface gpu fallback', () => {
  // the 2d atlas of 14³ fits in 64², wrapped (15³) it needs 128²
  const graphics = {
    capabilities: { floatRenderTarget: true, multipleRenderTargets: true },
    resources: { limits: { maxTextureSize: 128 } },
  } as unknown as GraphicsContext;
  registerAtlasIsosurfaceExtraction({ isSupported: (g) => g === graphics } as AtlasIsosurfaceExtraction);
  const key = { volume: createTestVolume(14), key: 0 };
  const props = (wrap: VolumeIsosurfaceProps['wrap']) => ({ ...PD.getDefaultValues(IsosurfaceMeshParams), wrap });

  it('falls back to the cpu once the wrapped volume is too large for the atlas', () => {
    const gpu = IsosurfaceVisual(0, key, props('off'), graphics);
    expect(gpu.mustRecreate?.(key, props('off'), graphics)).toBe(false);
    expect(gpu.mustRecreate?.(key, props('on'), graphics)).toBe(true);
    const cpu = IsosurfaceVisual(0, key, props('on'), graphics);
    expect(cpu.mustRecreate?.(key, props('on'), graphics)).toBe(false);
    expect(cpu.mustRecreate?.(key, props('off'), graphics)).toBe(true);
  });
});
