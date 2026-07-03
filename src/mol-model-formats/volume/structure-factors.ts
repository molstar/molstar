/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SF_Database } from '../../mol-io/reader/cif/schema/sf';
import { Volume } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { SpacegroupCell, Box3D, Spacegroup } from '../../mol-math/geometry';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { ModelFormat } from '../format';
import { CustomProperties } from '../../mol-model/custom-property';
import { degToRad } from '../../mol-math/misc';
import { computeElectronDensityFromReflections } from './shared/reflections';

export type StructureFactorMapType = '2fo-fc' | 'fo-fc';

export function volumeFromStructureFactors(
    source: SF_Database,
    params?: Partial<{ label: string, entryId: string, mapType: StructureFactorMapType }>
): Task<Volume> {
    return Task.create<Volume>('Compute Structure Factor Map', async ctx => {
        const mapType = params?.mapType ?? '2fo-fc';

        // 1. Parse cell parameters
        const { cell, symmetry, refln } = source;

        const spaceGroupId: number | string =
            symmetry.Int_Tables_number.value(0) ||
            symmetry['space_group_name_H-M'].value(0) ||
            'P 1';

        const cellSize = Vec3.create(
            cell.length_a.value(0),
            cell.length_b.value(0),
            cell.length_c.value(0)
        );
        const cellAngles = Vec3.create(
            degToRad(cell.angle_alpha.value(0)),
            degToRad(cell.angle_beta.value(0)),
            degToRad(cell.angle_gamma.value(0))
        );

        const sgCell = SpacegroupCell.create(spaceGroupId, cellSize, cellAngles);
        const sg = Spacegroup.create(sgCell);

        // 2. Extract reflection arrays
        await ctx.update({ message: 'Reading reflections...' });

        const count = refln.index_h.rowCount;
        if (count === 0) throw new Error('No reflection data found in the _refln category.');

        const hArr = refln.index_h.toArray({ array: Int16Array });
        const kArr = refln.index_k.toArray({ array: Int16Array });
        const lArr = refln.index_l.toArray({ array: Int16Array });

        const ampCol = (mapType === '2fo-fc' ? refln.pdbx_FWT : refln.pdbx_DELFWT).toArray({ array: Float32Array });
        const phiCol = (mapType === '2fo-fc' ? refln.pdbx_PHWT : refln.pdbx_DELPHWT).toArray({ array: Float32Array });

        // 3–5. Compute density via shared helper
        const { density, N0, N1, N2 } = await computeElectronDensityFromReflections(
            hArr, kArr, lArr, ampCol, phiCol, sg, ctx
        );

        // 6. Compute statistics
        await ctx.update({ message: 'Finalizing volume...' });
        const totalSize = N0 * N1 * N2;
        let min = Infinity, max = -Infinity, sum = 0, sum2 = 0;
        for (let i = 0; i < totalSize; i++) {
            const v = density[i];
            if (v < min) min = v;
            if (v > max) max = v;
            sum += v;
            sum2 += v * v;
        }
        const mean = sum / totalSize;
        const sigma = Math.sqrt(Math.max(0, sum2 / totalSize - mean * mean));

        // 7. Create Tensor and Volume
        // axisOrderSlowToFast = [0, 1, 2]: axis 0 is slowest (h), axis 2 is fastest (l)
        const tensorSpace = Tensor.Space([N0, N1, N2], [0, 1, 2], Float32Array);
        const data = Tensor.create(tensorSpace, Tensor.Data1(density));

        return {
            label: params?.label,
            entryId: params?.entryId,
            grid: {
                transform: {
                    kind: 'spacegroup',
                    cell: sgCell,
                    // The FFT grid covers the full unit cell [0, 1) in each direction
                    fractionalBox: Box3D.create(Vec3.zero(), Vec3.create(1, 1, 1)),
                },
                cells: data,
                stats: { min, max, mean, sigma },
                periodicity: 'xyz',
            },
            instances: [{ transform: Mat4.identity() }],
            sourceData: SfFormat.create(source),
            customProperties: new CustomProperties(),
            _propertyData: Object.create(null),
            _localPropertyData: Object.create(null),
        };
    });
}

//

export { SfFormat };

type SfFormat = ModelFormat<SF_Database>

namespace SfFormat {
    export function is(x?: ModelFormat): x is SfFormat {
        return x?.kind === 'sf';
    }

    export function create(db: SF_Database): SfFormat {
        return { kind: 'sf', name: db._name, data: db };
    }
}
