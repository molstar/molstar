/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Tensor } from '../../../mol-math/linear-algebra/tensor';
import { Mat4 } from '../../../mol-math/linear-algebra/3d/mat4';
import { Grid, Volume } from '../../../mol-model/volume';
import { ColorNames } from '../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { SyncRuntimeContext } from '../../../mol-task/execution/synchronous';
import { VtpFile } from '../../../mol-io/reader/vtp/schema';
import { createVtpShapeParams, shapeFromVtp } from '../vtp';

// ─── helpers ─────────────────────────────────────────────────────────────────

// One triangle whose vertices sit at x = 0, 1, 0 so volume sampling differs between them.
function makeTriangleVtpFile(): VtpFile {
    return {
        numberOfPoints: 3,
        numberOfCells: 1,
        positions: new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]),
        connectivity: new Int32Array([0, 1, 2]),
        numberOfTriangles: 1,
        triangleCellIndex: new Int32Array([0]),
        pointData: new Map(),
        cellData: new Map(),
    };
}

// 2×2×2 grid with an identity transform whose value equals the x index, so trilinear
// sampling at Cartesian x returns x on [0, 1].
function makeGradientVolume(): Volume {
    const space = Tensor.Space([2, 2, 2], [0, 1, 2]);
    const data = Tensor.Data1(new Float32Array(8));
    for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 2; j++) {
            for (let k = 0; k < 2; k++) {
                space.set(data, i, j, k, i);
            }
        }
    }
    const grid: Grid = {
        transform: { kind: 'matrix', matrix: Mat4.identity() },
        cells: Tensor.create(space, data),
        stats: { min: 0, max: 1, mean: 0.5, sigma: 0.5 },
    };
    // Only `grid` is read by the volume sampler; the rest of Volume is irrelevant here.
    return { grid } as Volume;
}

async function getVtpShape(vtpFile: VtpFile, volume?: Volume) {
    const params = createVtpShapeParams(vtpFile);
    const volumeParams = PD.getDefaultValues((params.colorTheme.map('external-volume') as PD.Group<any>).params);
    // Mock the resolved ValueRef the way state reconciliation would bind it.
    if (volume) volumeParams.volume = { ref: 'test', getValue: () => volume };
    const props = { ...PD.getDefaultValues(params), colorTheme: { name: 'external-volume' as const, params: volumeParams } };
    const provider = await shapeFromVtp(vtpFile).run();
    return provider.getShape(SyncRuntimeContext, provider.data, props as any);
}

// ─── tests ───────────────────────────────────────────────────────────────────

describe('VTP external-volume color option', () => {
    it('exposes an `external-volume` color option', () => {
        const names = createVtpShapeParams(makeTriangleVtpFile()).colorTheme.select.options.map(o => o[0]);
        expect(names).toContain('external-volume');
        expect(names).toContain('attribute');
    });

    it('colors vertices from the sampled volume value', async () => {
        const shape = await getVtpShape(makeTriangleVtpFile(), makeGradientVolume());

        const c0 = shape.getColor(0, 0); // x = 0
        const c1 = shape.getColor(1, 0); // x = 1
        // Distinct volume values ⇒ distinct colors, neither the grey fallback.
        expect(c0).not.toBe(c1);
        expect(c0).not.toBe(ColorNames.grey);
        expect(c1).not.toBe(ColorNames.grey);
    });

    it('falls back to grey when no volume is bound', async () => {
        const shape = await getVtpShape(makeTriangleVtpFile(), undefined);
        expect(shape.getColor(0, 0)).toBe(ColorNames.grey);
        expect(shape.getColor(1, 0)).toBe(ColorNames.grey);
    });
});
