/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../../mol-canvas3d/camera';
import { Sphere3D } from '../../../mol-math/geometry';
import { calcInstanceGrid } from '../../../mol-math/geometry/instance-grid';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { ValueCell } from '../../../mol-util/value-cell';
import { RenderableValues } from '../../renderable/schema';
import { getSphereLodInstanceRanges } from '../sphere-lod';

function fixture(depths: number[]) {
    const camera = new Camera({ position: Vec3.create(0, 0, 20), target: Vec3(), radius: 10, radiusMax: 10, fog: 0 }); camera.update();
    const transforms = new Float32Array(depths.length * 16);
    depths.forEach((depth, i) => { transforms.set(Mat4.identity(), i * 16); transforms[i * 16 + 14] = 20 - depth; });
    const values: RenderableValues = {
        aTransform: ValueCell.create(transforms), aInstance: ValueCell.create(Float32Array.from(depths, (_, i) => i)),
        uInvariantBoundingSphere: ValueCell.create([0, 0, 0, 1]),
    };
    return { camera, transforms, values };
}

describe('native sphere instance distance culling', () => {
    it('coalesces contiguous physical instances and preserves gaps', () => {
        const { camera, values } = fixture([20, 20, 140, 20, 40]);
        expect(getSphereLodInstanceRanges(values, 5, camera, 0, 25, false)).toEqual([{ first: 0, count: 2 }, { first: 3, count: 1 }]);
        expect(getSphereLodInstanceRanges(values, 5, camera, 25, 100, false)).toEqual([{ first: 4, count: 1 }]);
        expect(getSphereLodInstanceRanges(values, 5, camera, 200, 300, false)).toEqual([]);
    });

    it('uses reordered grid cells and falls back when transforms invalidate the grid', () => {
        const { camera, transforms, values } = fixture([20, 140, 40, 20]);
        const grid = calcInstanceGrid({ instanceCount: 4, instance: new Float32Array([0, 1, 2, 3]), transform: transforms, invariantBoundingSphere: Sphere3D.create(Vec3(), 1) }, 4, 20);
        values.instanceGrid = ValueCell.create(grid); ValueCell.update(values.aTransform, grid.cellTransform); ValueCell.update(values.aInstance, grid.cellInstance);
        const ranges = getSphereLodInstanceRanges(values, 4, camera, 0, 25, false);
        const ids = ranges.flatMap(r => Array.from(grid.cellInstance.subarray(r.first, r.first + r.count))).sort();
        expect(ids).toEqual([0, 3]);
        const changed = grid.cellTransform.slice();
        for (let i = 0; i < 4; i++) changed[i * 16 + 14] = 0;
        ValueCell.update(values.aTransform, changed);
        expect(getSphereLodInstanceRanges(values, 4, camera, 0, 25, false)).toEqual([{ first: 0, count: 4 }]);
    });

    it('keeps conservative animation and shear bounds even when grid spheres are insufficient', () => {
        const { camera, values } = fixture([27]);
        expect(getSphereLodInstanceRanges(values, 1, camera, 0, 25, false)).toEqual([]);
        values.uWiggleAmplitude = ValueCell.create(2);
        expect(getSphereLodInstanceRanges(values, 1, camera, 0, 25, true)).toEqual([{ first: 0, count: 1 }]);
        ValueCell.update(values.uWiggleAmplitude, 0); values.uTumbleAmplitude = ValueCell.create(4);
        expect(getSphereLodInstanceRanges(values, 1, camera, 0, 25, true)).toEqual([{ first: 0, count: 1 }]);
        const sheared = fixture([35]); sheared.transforms[6] = 10;
        const grid = calcInstanceGrid({ instanceCount: 1, instance: new Float32Array([0]), transform: sheared.transforms, invariantBoundingSphere: Sphere3D.create(Vec3(), 1) }, 4, 20);
        sheared.values.instanceGrid = ValueCell.create(grid); ValueCell.update(sheared.values.aTransform, grid.cellTransform); ValueCell.update(sheared.values.aInstance, grid.cellInstance);
        expect(getSphereLodInstanceRanges(sheared.values, 1, sheared.camera, 0, 25, false)).toEqual([{ first: 0, count: 1 }]);
    });
});
