/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { CIF } from '@molstar/io/reader/cif';
import { Volume } from '../../../model/volume.js';
import { volumeFromSegmentationData } from '../segmentation.js';

// a 2×2×2 grid: sets 1 and 2 each hold one segment; label 3 is in the grid but in no set
const SegCif = `data_seg
_volume_data_3d_info.name seg
_volume_data_3d_info.axis_order[0] 0
_volume_data_3d_info.axis_order[1] 1
_volume_data_3d_info.axis_order[2] 2
_volume_data_3d_info.origin[0] 0
_volume_data_3d_info.origin[1] 0
_volume_data_3d_info.origin[2] 0
_volume_data_3d_info.dimensions[0] 1
_volume_data_3d_info.dimensions[1] 1
_volume_data_3d_info.dimensions[2] 1
_volume_data_3d_info.sample_rate 1
_volume_data_3d_info.sample_count[0] 2
_volume_data_3d_info.sample_count[1] 2
_volume_data_3d_info.sample_count[2] 2
_volume_data_3d_info.spacegroup_number 1
_volume_data_3d_info.spacegroup_cell_size[0] 2
_volume_data_3d_info.spacegroup_cell_size[1] 2
_volume_data_3d_info.spacegroup_cell_size[2] 2
_volume_data_3d_info.spacegroup_cell_angles[0] 90
_volume_data_3d_info.spacegroup_cell_angles[1] 90
_volume_data_3d_info.spacegroup_cell_angles[2] 90
#
loop_
_segmentation_data_table.set_id
_segmentation_data_table.segment_id
1 1
2 2
#
loop_
_segmentation_data_3d.values
0 1 1 0 2 2 0 3
`;

describe('volumeFromSegmentationData', () => {
  it('skips voxel labels in no set (perf-audit-2 04/B13)', async () => {
    const parsed = await CIF.parse(SegCif).run();
    if (parsed.isError) throw new Error(parsed.message);
    const db = CIF.schema.segmentation(parsed.result.blocks[0]);
    const volume = await volumeFromSegmentationData(db).run();
    const segmentation = Volume.Segmentation.get(volume)!;
    expect([...segmentation.segments.keys()].sort()).toEqual([1, 2]);
    // segment 1: voxels 1 and 2 (x 1, then x 0 at y 1); segment 2: voxels 4 and 5 (z 1)
    expect([...segmentation.bounds[1 as Volume.SegmentIndex].min]).toEqual([0, 0, 0]);
    expect([...segmentation.bounds[1 as Volume.SegmentIndex].max]).toEqual([1, 1, 0]);
    expect([...segmentation.bounds[2 as Volume.SegmentIndex].min]).toEqual([0, 0, 1]);
    expect([...segmentation.bounds[2 as Volume.SegmentIndex].max]).toEqual([1, 0, 1]);
  });

  it('bounds sets in memory order, without getCoords per voxel (perf-audit-2 04/12/B1)', async () => {
    // a 3×4×5 grid stored with axis order 2, 0, 1 (fastest first); labels 1 and 2 scattered
    const [nx, ny, nz] = [3, 4, 5];
    const order = [2, 0, 1];
    const counts = order.map((a) => [nx, ny, nz][a]);
    const n = nx * ny * nz;
    const labels = Array.from({ length: n }, (_, i) => (i % 7 === 3 ? 1 : i % 11 === 5 ? 2 : 0));
    const info = (k: string, v: number[]) => v.map((x, i) => `_volume_data_3d_info.${k}[${i}] ${x}`).join('\n');
    const text = `data_seg
_volume_data_3d_info.name seg
${info('axis_order', order)}
${info('origin', [0, 0, 0])}
${info('dimensions', [1, 1, 1])}
_volume_data_3d_info.sample_rate 1
${info('sample_count', counts)}
_volume_data_3d_info.spacegroup_number 1
${info('spacegroup_cell_size', [nx, ny, nz])}
${info('spacegroup_cell_angles', [90, 90, 90])}
#
loop_
_segmentation_data_table.set_id
_segmentation_data_table.segment_id
1 1
2 2
#
loop_
_segmentation_data_3d.values
${labels.join(' ')}
`;
    const parsed = await CIF.parse(text).run();
    if (parsed.isError) throw new Error(parsed.message);
    const db = CIF.schema.segmentation(parsed.result.blocks[0]);

    const space = Tensor.Space;
    const getCoords: ReturnType<typeof vi.spyOn>[] = [];
    const spy = vi.spyOn(Tensor, 'Space').mockImplementation((...args: Parameters<typeof space>) => {
      const s = space(...args);
      getCoords.push(vi.spyOn(s, 'getCoords'));
      return s;
    });
    let volume;
    try {
      volume = await volumeFromSegmentationData(db).run();
    } finally {
      spy.mockRestore();
    }
    expect(getCoords.length).toBe(1);
    expect(getCoords[0]).not.toHaveBeenCalled();

    // the bounds by a scan through getCoords
    const { space: s, data } = volume.grid.cells;
    const c = [0, 0, 0];
    for (const seg of [1, 2] as Volume.SegmentIndex[]) {
      const min = [Infinity, Infinity, Infinity],
        max = [-1, -1, -1];
      for (let i = 0; i < data.length; i++) {
        if (data[i] !== seg) continue;
        s.getCoords(i, c);
        for (let a = 0; a < 3; a++) {
          min[a] = Math.min(min[a], c[a]);
          max[a] = Math.max(max[a], c[a]);
        }
      }
      expect([...Volume.Segmentation.get(volume)!.bounds[seg].min]).toEqual(min);
      expect([...Volume.Segmentation.get(volume)!.bounds[seg].max]).toEqual(max);
    }
  });
});
