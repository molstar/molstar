/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Camera } from '../mol-canvas3d/camera';
import { Vec3 } from '../mol-math/linear-algebra';

/** Return the distance adjustment ratio for conversion from the "reference camera"
 * to a camera with an arbitrary field of view `fov`. */
function distanceAdjustment(mode: Camera.Mode, fov: number) {
    if (mode === 'orthographic') return 1 / (2 * Math.tan(fov / 2));
    else return 1 / (2 * Math.sin(fov / 2));
}

/** Return the position for a camera with an arbitrary field of view `fov`
 * necessary to just fit into view the same sphere (with center at `target`)
 * as the "reference camera" placed at `refPosition` would fit, while keeping the camera orientation.
 * The "reference camera" is a camera which can just fit into view a sphere of radius R with center at distance 2R
 * (this corresponds to FOV = 2 * asin(1/2) in perspective mode or FOV = 2 * atan(1/2) in orthographic mode). */
export function fovAdjustedPosition(target: Vec3, refPosition: Vec3, mode: Camera.Mode, fov: number) {
    const delta = Vec3.sub(Vec3(), refPosition, target);
    const adjustment = distanceAdjustment(mode, fov);
    return Vec3.scaleAndAdd(delta, target, delta, adjustment); // return target + delta * adjustment
}

/** Return the inverse of fovAdjustedPosition to be able to store invariant camera position,
 * e.g., in MolViewSpec snapshots.
 */
export function fovNormalizedCameraPosition(target: Vec3, refPosition: Vec3, mode: Camera.Mode, fov: number) {
    const delta = Vec3.sub(Vec3(), refPosition, target);
    const adjustment = distanceAdjustment(mode, fov) || 1;
    return Vec3.scaleAndAdd(delta, target, delta, 1 / adjustment); // return target + delta / adjustment
}