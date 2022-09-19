/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
    state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    viewOffset: Readonly<Camera.ViewOffset> = Camera.ViewOffset();
    far: number = 0;
    near: number = 0;
    fogFar: number = 0;
    fogNear: number = 0;
}

const eyeLeft = Mat4.identity(), eyeRight = Mat4.identity();

function update(camera: Camera, props: StereoCameraProps, left: EyeCamera, right: EyeCamera) {
    // Copy the states

    Viewport.copy(left.viewport, camera.viewport);
    Mat4.copy(left.view, camera.view);
    Mat4.copy(left.projection, camera.projection);
    Camera.copySnapshot(left.state, camera.state);
    Camera.copyViewOffset(left.viewOffset, camera.viewOffset);
    left.far = camera.far;
    left.near = camera.near;
    left.fogFar = camera.fogFar;
    left.fogNear = camera.fogNear;

    Viewport.copy(right.viewport, camera.viewport);
    Mat4.copy(right.view, camera.view);
    Mat4.copy(right.projection, camera.projection);
    Camera.copySnapshot(right.state, camera.state);
    Camera.copyViewOffset(right.viewOffset, camera.viewOffset);
    right.far = camera.far;
    right.near = camera.near;
    right.fogFar = camera.fogFar;
    right.fogNear = camera.fogNear;

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

    eyeLeft[12] = -eyeSepHalf;
    eyeRight[12] = eyeSepHalf;

    // for left eye

    xmin = -ymax * aspect + eyeSepOnProjection;
    xmax = ymax * aspect + eyeSepOnProjection;

    left.projection[0] = 2 * camera.near / (xmax - xmin);
    left.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(left.view, left.view, eyeLeft);
    Mat4.mul(left.projectionView, left.projection, left.view);
    Mat4.invert(left.inverseProjectionView, left.projectionView);

    // for right eye

    xmin = -ymax * aspect - eyeSepOnProjection;
    xmax = ymax * aspect - eyeSepOnProjection;

    right.projection[0] = 2 * camera.near / (xmax - xmin);
    right.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(right.view, right.view, eyeRight);
    Mat4.mul(right.projectionView, right.projection, right.view);
    Mat4.invert(right.inverseProjectionView, right.projectionView);
}