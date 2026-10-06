/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Camera } from '../../mol-canvas3d/camera';
import { Mat4 } from '../../mol-math/linear-algebra';

/** Transform molecular coordinates into the scaled camera's view space. */
export function getWebGPUModelView(camera: Camera) {
    return Mat4.scaleUniformly(Mat4(), camera.view, camera.scale);
}

/** Camera snapshots omit runtime properties used by rendering and picking. */
export function copyWebGPUCameraState(target: Camera, source: Camera) {
    target.setState(source.getSnapshot(), 0);
    target.scale = source.scale;
    target.forceFull = source.forceFull;
    target.minTargetDistance = source.minTargetDistance;
    Mat4.copy(target.headRotation, source.headRotation);
    Mat4.copy(target.viewEye, source.viewEye);
}
