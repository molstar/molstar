/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Box3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { type Grid, Volume } from '@molstar/model/model/volume';
import { createVolumeSegmentMesh } from '../segment.js';
import { createSegmentTexture2d, getSegmentMaskValue } from '../util.js';

const N = 8;

/** Segment 1 (label 1) is a slab through the whole grid along x: it touches the x = 0 and x = N - 1 faces */
function createSlab(): Volume {
  const space = Tensor.Space([N, N, N], [2, 1, 0], Float32Array);
  const data = space.create();
  for (let k = 2; k <= 5; k++) {
    for (let j = 2; j <= 5; j++) {
      for (let i = 0; i < N; i++) space.set(data, i, j, k, 1);
    }
  }
  const volume: Volume = {
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, data),
      stats: { min: 0, max: 1, mean: 0, sigma: 1 },
    } satisfies Grid,
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
  const segment = 1 as Volume.SegmentIndex;
  Volume.Segmentation.set(volume, {
    segments: new Map([[segment, new Set([1])]]),
    sets: new Map([[1, new Set([segment])]]),
    bounds: { [segment]: Box3D.create(Vec3.create(0, 2, 2), Vec3.create(N - 1, 5, 5)) },
    labels: { [segment]: 'slab' },
  });
  return volume;
}

describe('volume segment mesh', () => {
  it('reads nothing past the grid faces (perf-audit-2 07/N1)', async () => {
    const mesh = await createVolumeSegmentMesh(
      { runtime: SyncRuntimeContext } as any,
      createSlab(),
      1 as Volume.SegmentIndex,
      undefined as any,
      {} as any,
    );
    expect(mesh.vertexCount).toBeGreaterThan(0);

    // grid space = model space here (identity grid transform): the surface closes within a voxel of each face
    const v = mesh.vertexBuffer.ref.value;
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < mesh.vertexCount; i++) {
      for (let c = 0; c < 3; c++) {
        min[c] = Math.min(min[c], v[3 * i + c]);
        max[c] = Math.max(max[c], v[3 * i + c]);
      }
    }
    expect(min[0]).toBeGreaterThanOrEqual(-1);
    expect(max[0]).toBeLessThanOrEqual(N);
    expect(min[1]).toBeGreaterThanOrEqual(1);
    expect(max[1]).toBeLessThanOrEqual(6);
  });

  it('the WebGL segment texture reads nothing past the grid faces either (perf-audit-2 07/N1)', () => {
    const volume = createSlab();
    // the slab's box, expanded by 2 as the visuals do: x from -2 to N + 1
    const bbox = Box3D.create(Vec3.create(-2, 0, 0), Vec3.create(N + 2, 8, 8));
    const { array, width } = createSegmentTexture2d(volume, [1], bbox);
    const [xn, yn, zn] = Box3D.size(Vec3(), bbox);
    let inside = 0;
    for (let z = 0; z < zn; ++z) {
      for (let y = 0; y < yn; ++y) {
        for (let x = 0; x < xn; ++x) {
          const column = Math.floor(((z * xn) % width) / xn);
          const row = Math.floor((z * xn) / width);
          const v = array[row * yn * width + y * width + column * xn + x];
          const i = x - 2;
          // one voxel past a face holds the closing value; two past, nothing
          if (i < -1 || i > N) expect(v, `x ${i}`).toBe(0);
          else if (v > 0) inside++;
        }
      }
    }
    expect(inside).toBeGreaterThan(0);
  });
});

describe('segment mask', () => {
  it('tests each voxel of the box once, with the same values (perf-audit-2 07/P8)', () => {
    const { cells } = createSlab().grid;
    // past every face, as the visuals' expanded boxes are
    const bbox = Box3D.create(Vec3.create(-2, 0, 1), Vec3.create(N + 2, 7, 8));
    const set = [1];
    const includes = vi.spyOn(set, 'includes');
    const boxed = getSegmentMaskValue(cells, set, bbox);
    const plain = getSegmentMaskValue(cells, [1]);
    let n = 0;
    for (let z = bbox.min[2]; z < bbox.max[2]; z++) {
      for (let y = bbox.min[1]; y < bbox.max[1]; y++) {
        for (let x = bbox.min[0]; x < bbox.max[0]; x++) {
          expect(boxed(x, y, z)).toBe(plain(x, y, z));
          n++;
        }
      }
    }
    // the box and its border, inside the grid; 7 per voxel before
    expect(includes.mock.calls.length).toBeLessThanOrEqual(n);
  });
});
