/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Camera, ICamera } from '../../mol-canvas3d/camera';
import { StereoCamera, StereoCameraProps } from '../../mol-canvas3d/camera/stereo';
import { Mat4 } from '../../mol-math/linear-algebra';
import { copyWebGPUCameraState } from './camera';

/** An asymmetric eye projection that remains intact when supersampling updates the camera. */
class WebGPUEyeCamera extends Camera {
    override readonly isAsymmetricProjection = true;
    private readonly baseView = Mat4.identity();
    private readonly baseProjection = Mat4.identity();

    setEye(eye: ICamera, offset: Camera.ViewOffset) {
        Camera.copySnapshot(this.state, eye.state);
        Object.assign(this.viewport, eye.viewport);
        Object.assign(this.viewOffset, offset);
        Mat4.copy(this.baseView, eye.view); Mat4.copy(this.baseProjection, eye.projection);
        Mat4.copy(this.headRotation, eye.headRotation); Mat4.copy(this.viewEye, eye.viewEye);
        this.scale = eye.scale; this.forceFull = eye.forceFull; this.minTargetDistance = eye.minTargetDistance;
        this.near = eye.near; this.far = eye.far; this.fogNear = eye.fogNear; this.fogFar = eye.fogFar;
        this.update();
    }

    override update() {
        const previousView = Mat4.clone(this.view), previousProjection = Mat4.clone(this.projection);
        Mat4.copy(this.view, this.baseView); Mat4.copy(this.projection, this.baseProjection);
        const v = this.viewOffset;
        if (v.enabled) {
            const sx = v.fullWidth / v.width, sy = v.fullHeight / v.height;
            const tx = (v.fullWidth - v.width - 2 * v.offsetX) / v.width;
            const ty = (2 * v.offsetY + v.height - v.fullHeight) / v.height;
            for (let column = 0; column < 4; column++) {
                const i = column * 4;
                this.projection[i] = sx * this.baseProjection[i] + tx * this.baseProjection[i + 3];
                this.projection[i + 1] = sy * this.baseProjection[i + 1] + ty * this.baseProjection[i + 3];
            }
        }
        Mat4.mul(this.projectionView, this.projection, this.view);
        Mat4.invert(this.inverseProjectionView, this.projectionView);
        const changed = !Mat4.areEqual(previousView, this.view, 1e-6) || !Mat4.areEqual(previousProjection, this.projection, 1e-6);
        if (changed) this.changed.next();
        return changed;
    }
}

export class WebGPUStereoCamera {
    readonly left = new WebGPUEyeCamera();
    readonly right = new WebGPUEyeCamera();
    private readonly parent = new Camera();
    private readonly stereo = new StereoCamera(this.parent);

    update(source: Camera, props: StereoCameraProps) {
        copyWebGPUCameraState(this.parent, source);
        Object.assign(this.parent.viewport, source.viewport);
        this.parent.viewOffset.enabled = false;
        this.stereo.setProps(props); this.stereo.update();
        this.left.setEye(this.stereo.left, source.viewOffset);
        this.right.setEye(this.stereo.right, source.viewOffset);
    }
}
