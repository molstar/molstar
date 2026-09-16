/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Projecting voxels into the camera a `ViewMask` was drawn in — shared by the volume tools.
 */

import { Camera } from '../../mol-canvas3d/camera';
import { Viewport, cameraProject } from '../../mol-canvas3d/camera/util';
import { Mat4, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import type { ViewMask } from './types';

const tmpVec4 = Vec4();

/** Pre-computed camera data for a single ViewMask — built once before the voxel loop. */
export interface PreparedMask {
    projectionView: Mat4;
    viewport: Viewport;
    normPolygon: [number, number][];
    inverted: boolean;
}

export function prepareMask(mask: ViewMask): PreparedMask {
    const viewport = Viewport.create(0, 0, mask.viewportWidth, mask.viewportHeight);
    const cam = new Camera(mask.cameraSnapshot, viewport);
    cam.update();
    // Deep-copy the matrix (Camera keeps it as a mutable property)
    const projectionView = Mat4.copy(Mat4(), cam.projectionView);
    const w = mask.canvasWidth, h = mask.canvasHeight;
    const normPolygon = mask.polygon.map(([px, py]) => [px / w, py / h] as [number, number]);
    return { projectionView, viewport, normPolygon, inverted: !!mask.inverted };
}

/** Projects a world-space point to normalised [0,1] canvas coords (Y flipped), written to `out`. */
export function projectToNormInPlace(worldPos: Vec3, prepared: PreparedMask, out: [number, number]): void {
    cameraProject(tmpVec4, worldPos, prepared.viewport, prepared.projectionView);
    out[0] = tmpVec4[0] / prepared.viewport.width;
    out[1] = 1 - tmpVec4[1] / prepared.viewport.height;
}
