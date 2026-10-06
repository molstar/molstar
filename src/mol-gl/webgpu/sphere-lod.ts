/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */
import { Camera } from '../../mol-canvas3d/camera';
import { InstanceGrid } from '../../mol-math/geometry/instance-grid';
import { Vec3 } from '../../mol-math/linear-algebra';
import { RenderableValues } from '../renderable/schema';
import { getWebGPUModelView } from './camera';
import { value } from './geometry';

export interface SphereInstanceRange { first: number, count: number }

/** Physical instance ranges retain builtin instance indices and logical picking IDs. */
export function getSphereLodInstanceRanges(values: RenderableValues, instances: number, camera: Camera, min: number, max: number, animation: boolean): SphereInstanceRange[] {
    const transforms = value(values, 'aTransform', new Float32Array(0));
    const invariant = value(values, 'uInvariantBoundingSphere', [0, 0, 0, 0]);
    const grid = value<InstanceGrid | undefined>(values, 'instanceGrid', undefined);
    const modelView = getWebGPUModelView(camera), point = Vec3(), ranges: SphereInstanceRange[] = [];
    const wiggle = animation ? Math.max(0, value(values, 'uWiggleAmplitude', 0)) + (value(values, 'dWiggle', false) ? Math.max(0, value(values, 'uWiggleStrength', 1)) : 0) : 0;
    const tumble = animation ? Math.sqrt(3) * Math.max(0, value(values, 'uTumbleAmplitude', 0)) / Math.max(invariant[3], 1) : 0;
    const emit = (first: number, count: number) => {
        const previous = ranges[ranges.length - 1];
        if (previous && previous.first + previous.count === first) previous.count += count;
        else ranges.push({ first, count });
    };
    const intersects = (x: number, y: number, z: number, radius: number) => {
        const distance = -(modelView[2] * x + modelView[6] * y + modelView[10] * z + modelView[14]) / camera.scale;
        return distance + radius >= min && distance - radius <= max;
    };
    const inspect = (first: number, end: number) => {
        for (let i = first; i < end; i++) {
            const o = i * 16;
            Vec3.set(point, invariant[0], invariant[1], invariant[2]);
            Vec3.transformMat4Offset(point, point, transforms, 0, 0, o);
            // sqrt(||A||1 ||A||inf) bounds the largest singular value, including shear.
            let column = 0, row = 0;
            for (let a = 0; a < 3; a++) {
                column = Math.max(column, Math.abs(transforms[o + a * 4]) + Math.abs(transforms[o + a * 4 + 1]) + Math.abs(transforms[o + a * 4 + 2]));
                row = Math.max(row, Math.abs(transforms[o + a]) + Math.abs(transforms[o + a + 4]) + Math.abs(transforms[o + a + 8]));
            }
            const scale = Math.sqrt(column * row);
            const radius = (invariant[3] + Math.sqrt(3) * wiggle) * scale + tumble;
            if (intersects(point[0], point[1], point[2], radius)) emit(i, 1);
        }
    };
    let rigid = true;
    for (let i = 0; i < instances && rigid; i++) {
        const o = i * 16;
        for (let a = 0; a < 3; a++) for (let b = a; b < 3; b++) {
            let dot = 0;
            for (let c = 0; c < 3; c++) dot += transforms[o + a * 4 + c] * transforms[o + b * 4 + c];
            if (a === b ? dot > 1.000001 : Math.abs(dot) > 0.000001) rigid = false;
        }
    }
    const sourceBounds = grid?.invariantBoundingSphere;
    const sameBounds = sourceBounds && sourceBounds.radius === invariant[3] && sourceBounds.center.every((v, i) => v === invariant[i]);
    const validGrid = sameBounds && rigid && !wiggle && !tumble && grid && grid.cellCount > 0 && grid.cellTransform === transforms && grid.cellInstance === value(values, 'aInstance', undefined) && grid.cellOffsets[grid.cellCount] === instances;
    if (validGrid) {
        for (let cell = 0; cell < grid.cellCount; cell++) {
            const o = cell * 4;
            if (intersects(grid.cellSpheres[o], grid.cellSpheres[o + 1], grid.cellSpheres[o + 2], grid.cellSpheres[o + 3] * 1.000004)) inspect(grid.cellOffsets[cell], grid.cellOffsets[cell + 1]);
        }
    } else inspect(0, instances);
    return ranges;
}
