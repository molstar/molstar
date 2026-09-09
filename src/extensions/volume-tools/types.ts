/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Types shared by the volume tools.
 */

import { Camera } from '../../mol-canvas3d/camera';

export type Point2D = [number, number];

/**
 * A 2D polygon drawn on top of the viewport at a specific camera orientation.
 * Everything needed to project any voxel into this view is captured at draw time.
 */
export interface ViewMask {
    id: string;
    label: string;
    /** Polygon vertices in CSS pixel coords on the overlay canvas (y=0 at top). */
    polygon: Point2D[];
    /** CSS pixel dimensions of the overlay canvas at draw time. */
    canvasWidth: number;
    canvasHeight: number;
    /** Physical pixel dimensions of the WebGL canvas at draw time (may differ by dpr). */
    viewportWidth: number;
    viewportHeight: number;
    /** Full camera state frozen at draw time. */
    cameraSnapshot: Camera.Snapshot;
    /** When true, voxels OUTSIDE this polygon are selected instead of inside. */
    inverted?: boolean;
}
