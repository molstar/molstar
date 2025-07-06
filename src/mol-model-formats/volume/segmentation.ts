/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { SpacegroupCell, Box3D } from '../../mol-math/geometry';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';
import { Segmentation_Data_Database } from '../../mol-io/reader/cif/schema/segmentation';
import { objectForEach } from '../../mol-util/object';

export function volumeFromSegmentationData(source: Segmentation_Data_Database, params?: Partial<{ label: string, segmentLabels: { [id: number]: string }, ownerId: string }>): Task<Volume> {
    return Task.create<Volume>('Create Segmentation Volume', async ctx => {
        const { volume_data_3d_info: info, segmentation_data_3d: values } = source;
        const cell = SpacegroupCell.create(
            info.spacegroup_number.value(0),
            Vec3.ofArray(info.spacegroup_cell_size.value(0)),
            Vec3.scale(Vec3(), Vec3.ofArray(info.spacegroup_cell_angles.value(0)), Math.PI / 180)
        );

        const axis_order_fast_to_slow = info.axis_order.value(0);

        const normalizeOrder = Tensor.convertToCanonicalAxisIndicesFastToSlow(axis_order_fast_to_slow);

        // sample count is in "axis order" and needs to be reordered
        const sample_count = normalizeOrder(info.sample_count.value(0));
        const tensorSpace = Tensor.Space(sample_count, Tensor.invertAxisOrder(axis_order_fast_to_slow), Float32Array);

        const t = Tensor.create(tensorSpace, Tensor.Data1(values.values.toArray({ array: Float32Array })));

        // origin and dimensions are in "axis order" and need to be reordered
        const origin = Vec3.ofArray(normalizeOrder(info.origin.value(0)));
        const dimensions = Vec3.ofArray(normalizeOrder(info.dimensions.value(0)));

        const v: Volume = {
            label: params?.label,
            entryId: undefined,
            grid: {
                transform: {
                    kind: 'spacegroup',
                    cell,
                    fractionalBox: Box3D.create(origin, Vec3.add(Vec3(), origin, dimensions))
                },
                cells: t,
                stats: {
                    min: 0, max: 1, mean: 0, sigma: 1
                },
            },
            instances: [{ transform: Mat4.identity() }],
            sourceData: SegcifFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: { ownerId: params?.ownerId },
        };

        Volume.PickingGranularity.set(v, 'object');

        const segments = new Map<Volume.SegmentIndex, Set<number>>();
        const sets = new Map<number, Set<Volume.SegmentIndex>>();
        const { segment_id, set_id } = source.segmentation_data_table;
        for (let i = 0, il = segment_id.rowCount; i < il; ++i) {
            const segment = segment_id.value(i) as Volume.SegmentIndex;
            const set = set_id.value(i);
            if (set === 0 || segment === 0) continue;

            if (!sets.has(set)) sets.set(set, new Set());
            sets.get(set)!.add(segment);
        }
        sets.forEach((segs, set) => {
            segs.forEach(seg => {
                if (!segments.has(seg)) segments.set(seg, new Set());
                segments.get(seg)!.add(set);
            });
        });

        const c = [0, 0, 0];
        const getCoords = t.space.getCoords;
        const d = t.data;
        const [xn, yn, zn] = v.grid.cells.space.dimensions;
        const xn1 = xn - 1;
        const yn1 = yn - 1;
        const zn1 = zn - 1;

        const setBounds: { [k: number]: [number, number, number, number, number, number] } = {};
        sets.forEach((v, k) => {
            setBounds[k] = [xn1, yn1, zn1, -1, -1, -1];
        });

        for (let i = 0, il = d.length; i < il; ++i) {
            const v = d[i];
            if (v === 0) continue;

            getCoords(i, c);
            const b = setBounds[v];
            if (c[0] < b[0]) b[0] = c[0];
            if (c[1] < b[1]) b[1] = c[1];
            if (c[2] < b[2]) b[2] = c[2];
            if (c[0] > b[3]) b[3] = c[0];
            if (c[1] > b[4]) b[4] = c[1];
            if (c[2] > b[5]) b[5] = c[2];
        }

        const bounds: { [k: Volume.SegmentIndex]: Box3D } = {};
        segments.forEach((v, k) => {
            bounds[k] = Box3D.create(Vec3.create(xn1, yn1, zn1), Vec3.create(-1, -1, -1));
        });

        objectForEach(setBounds, (b, s) => {
            sets.get(parseInt(s))!.forEach(seg => {
                const sb = bounds[seg];
                if (b[0] < sb.min[0]) sb.min[0] = b[0];
                if (b[1] < sb.min[1]) sb.min[1] = b[1];
                if (b[2] < sb.min[2]) sb.min[2] = b[2];
                if (b[3] > sb.max[0]) sb.max[0] = b[3];
                if (b[4] > sb.max[1]) sb.max[1] = b[4];
                if (b[5] > sb.max[2]) sb.max[2] = b[5];
            });
        });

        Volume.Segmentation.set(v, { segments, sets, bounds, labels: params?.segmentLabels ?? {} });

        return v;
    });
}

//

export { SegcifFormat };

type SegcifFormat = ModelFormat<Segmentation_Data_Database>

namespace SegcifFormat {
    export function is(x?: ModelFormat): x is SegcifFormat {
        return x?.kind === 'segcif';
    }

    export function create(segcif: Segmentation_Data_Database): SegcifFormat {
        return { kind: 'segcif', name: segcif._name, data: segcif };
    }
}