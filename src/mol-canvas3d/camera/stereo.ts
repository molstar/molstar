/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from three.js, The MIT License, Copyright © 2010-2020 three.js authors
 */

import { Ray3D } from '../../mol-math/geometry/primitives/ray3d';
import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera, ICamera } from '../camera';
import { cameraUnproject, Viewport } from './util';

export const StereoCameraParams = {
    eyeSeparation: PD.Numeric(0.062, { min: 0.02, max: 0.1, step: 0.001 }, { description: 'Distance between left and right camera.' }),
    focus: PD.Numeric(10, { min: 1, max: 20, step: 0.1 }, { description: 'Apparent object distance.' }),
};
export const DefaultStereoCameraProps = PD.getDefaultValues(StereoCameraParams);
export type StereoCameraProps = PD.Values<typeof StereoCameraParams>

export { StereoCamera };

class StereoCamera {
    readonly left = new EyeCamera();
    readonly right = new EyeCamera();

    get viewport() {
        return this.parent.viewport;
    }

    get viewOffset() {
        return this.parent.viewOffset;
    }

    private props: StereoCameraProps;

    constructor(private parent: Camera, props: Partial<StereoCameraProps> = {}) {
        this.props = { ...DefaultStereoCameraProps, ...props };
    }

    setProps(props: Partial<StereoCameraProps>) {
        Object.assign(this.props, props);
    }

    update(xr?: { pose: XRViewerPose, layer: XRWebGLLayer }) {
        this.parent.update();
        if (xr) {
            xrUpdate(this.parent, this.left, this.right, xr);
        } else {
            update(this.parent, this.props, this.left, this.right);
        }
    }
}

namespace StereoCamera {
    export function is(camera: Camera | StereoCamera): camera is StereoCamera {
        return 'left' in camera && 'right' in camera;
    }
}

class EyeCamera implements ICamera {
    viewport = Viewport.create(0, 0, 0, 0);
    view = Mat4();
    projection = Mat4();
    projectionView = Mat4();
    inverseProjectionView = Mat4();
    headRotation = Mat4();
    viewEye = Mat4();
    isAsymmetricProjection = true;

    state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    viewOffset: Readonly<Camera.ViewOffset> = Camera.ViewOffset();
    far: number = 0;
    near: number = 0;
    fogFar: number = 0;
    fogNear: number = 0;

    forceFull: boolean = false;
    scale: number = 0;
    minTargetDistance: number = 0;

    disabled = false;

    getRay(out: Ray3D, x: number, y: number) {
        Mat4.getTranslation(out.origin, Mat4.invert(Mat4(), this.view));
        Vec3.set(out.direction, x, y, 0.5);
        cameraUnproject(out.direction, out.direction, this.viewport, this.inverseProjectionView);
        Vec3.normalize(out.direction, Vec3.sub(out.direction, out.direction, out.origin));
        return out;
    }
}

const tmpEyeLeft = Mat4.identity();
const tmpEyeRight = Mat4.identity();

function copyStates(parent: Camera, eye: EyeCamera) {
    Viewport.copy(eye.viewport, parent.viewport);
    Mat4.copy(eye.view, parent.view);
    Mat4.copy(eye.projection, parent.projection);
    Mat4.copy(eye.headRotation, parent.headRotation);
    Camera.copySnapshot(eye.state, parent.state);
    Camera.copyViewOffset(eye.viewOffset, parent.viewOffset);
    eye.far = parent.far;
    eye.near = parent.near;
    eye.fogFar = parent.fogFar;
    eye.fogNear = parent.fogNear;

    eye.forceFull = parent.forceFull;
    eye.scale = parent.scale;
    eye.minTargetDistance = parent.minTargetDistance;
}

//

function update(camera: Camera, props: StereoCameraProps, left: EyeCamera, right: EyeCamera) {
    // Copy the states

    copyStates(camera, left);
    copyStates(camera, right);

    // update the view offsets

    const w = Math.floor(camera.viewport.width / 2);
    const aspect = w / camera.viewport.height;

    left.viewport.width = w;
    right.viewport.x += w;
    right.viewport.width -= w;

    // update the projection and view matrices

    const eyeSepHalf = props.eyeSeparation / 2;
    const eyeSepOnProjection = eyeSepHalf * camera.near / props.focus;
    const ymax = camera.near * Math.tan(camera.state.fov * 0.5);
    let xmin, xmax;

    // translate xOffset

    tmpEyeLeft[12] = -eyeSepHalf;
    tmpEyeRight[12] = eyeSepHalf;

    // for left eye

    xmin = -ymax * aspect + eyeSepOnProjection;
    xmax = ymax * aspect + eyeSepOnProjection;

    left.projection[0] = 2 * camera.near / (xmax - xmin);
    left.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(left.view, left.view, tmpEyeLeft);
    Mat4.mul(left.projectionView, left.projection, left.view);
    Mat4.invert(left.inverseProjectionView, left.projectionView);

    // for right eye

    xmin = -ymax * aspect - eyeSepOnProjection;
    xmax = ymax * aspect - eyeSepOnProjection;

    right.projection[0] = 2 * camera.near / (xmax - xmin);
    right.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(right.view, right.view, tmpEyeRight);
    Mat4.mul(right.projectionView, right.projection, right.view);
    Mat4.invert(right.inverseProjectionView, right.projectionView);

    // ensure enabled

    left.disabled = false;
    right.disabled = false;
}

//

function xrUpdate(camera: Camera, left: EyeCamera, right: EyeCamera, xr: { pose: XRViewerPose, layer: XRWebGLLayer }) {
    _xrUpdate(camera, left, xr.pose.views[0], xr.layer);
    if (xr.pose.views.length === 1) {
        right.disabled = true;
    } else {
        _xrUpdate(camera, right, xr.pose.views[1], xr.layer);
    }
}

function _xrUpdate(camera: Camera, eye: EyeCamera, view: XRView, layer: XRWebGLLayer) {
    copyStates(camera, eye);
    const lvp = layer.getViewport(view)!;
    Viewport.set(eye.viewport, lvp.x, lvp.y, lvp.width, lvp.height);
    Mat4.fromArray(eye.projection, view.projectionMatrix, 0);
    Mat4.fromArray(eye.view, view.transform.inverse.matrix, 0);
    Mat4.mul(eye.projectionView, eye.projection, eye.view);
    Mat4.invert(eye.inverseProjectionView, eye.projectionView);
    eye.disabled = false;
}
