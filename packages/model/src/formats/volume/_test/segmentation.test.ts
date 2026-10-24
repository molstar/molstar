/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
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
});
