/**
 * Copyright (c) 2020-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * Adapted from three.js, The MIT License, Copyright © 2010-2020 three.js authors
 */

import { Mat4 } from '../../mol-math/linear-algebra';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Camera, ICamera } from '../camera';
import { Viewport } from './util';

export const StereoCameraParams = {
    eyeSeparation: PD.Numeric(0.062, { min: 0.02, max: 0.1, step: 0.001 }, { description: 'Distance between left and right camera.' }),
    focus: PD.Numeric(10, { min: 1, max: 20, step: 0.1 }, { description: 'Apparent object distance.' }),
};
export const DefaultStereoCameraProps = PD.getDefaultValues(StereoCameraParams);
export type StereoCameraProps = PD.Values<typeof StereoCameraParams>

export { StereoCamera };

class StereoCamera {
    readonly left: ICamera = new EyeCamera();
    readonly right: ICamera = new EyeCamera();

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

    update() {
        this.parent.update();
        update(this.parent, this.props, this.left as EyeCamera, this.right as EyeCamera);
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
    state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    viewOffset: Readonly<Camera.ViewOffset> = Camera.ViewOffset();
    far: number = 0;
    near: number = 0;
    fogFar: number = 0;
    fogNear: number = 0;
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
}