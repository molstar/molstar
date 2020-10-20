/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from three.js, The MIT License, Copyright © 2010-2020 three.js authors
 */

import { Mat4 } from '../../mol-math/linear-algebra';
import { ParamDefinition } from '../../mol-util/param-definition';
import { Camera, ICamera } from '../camera';
import { Viewport } from './util';

export class StereoCamera {
    readonly left: ICamera = new EyeCamera();
    readonly right: ICamera = new EyeCamera();

    update(camera: Camera, params: StereoCameraParams) {
        update(camera, params, this.left as EyeCamera, this.right as EyeCamera);
    }
}

class EyeCamera implements ICamera {
    viewport = Viewport.create(0, 0, 0, 0);
    view = Mat4();
    projection = Mat4();
    state: Readonly<Camera.Snapshot> = Camera.createDefaultSnapshot();
    viewOffset: Readonly<Camera.ViewOffset> = Camera.ViewOffset();
    far: number = 0;
    near: number = 0;
    fogFar: number = 0;
    fogNear: number = 0;
}

export const StereoCameraParams = {
    aspect: ParamDefinition.Numeric(1, { min: 0.1, max: 3, step: 0.01 }),
    eyeSeparation: ParamDefinition.Numeric(0.064, { min: 0.01, max: 0.5, step: 0.001 }),
    focus: ParamDefinition.Numeric(10, { min: 1, max: 100, step: 0.01 }),
};
export type StereoCameraParams = ParamDefinition.Values<typeof StereoCameraParams>

const eyeLeft = Mat4.identity(), eyeRight = Mat4.identity();

function update(camera: ICamera, params: StereoCameraParams, left: EyeCamera, right: EyeCamera) {
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
    let w = (camera.viewport.width / 2) | 0;

    left.viewport.width = w;
    right.viewport.x = w;
    right.viewport.width -= w;

    // update the projection and view matrices

    const eyeSepHalf = params.eyeSeparation / 2;
    const eyeSepOnProjection = eyeSepHalf * camera.near / params.focus;
    const ymax = (camera.near * Math.tan(camera.state.fov * 0.5)) / /* cache.zoom */ 1;
    let xmin, xmax;

    // translate xOffset

    eyeLeft[12] = - eyeSepHalf;
    eyeRight[12] = eyeSepHalf;

    // for left eye

    xmin = - ymax * params.aspect + eyeSepOnProjection;
    xmax = ymax * params.aspect + eyeSepOnProjection;

    left.projection[0] = 2 * camera.near / (xmax - xmin);
    left.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(left.view, left.view, eyeLeft);

    // for right eye

    xmin = - ymax * params.aspect - eyeSepOnProjection;
    xmax = ymax * params.aspect - eyeSepOnProjection;

    right.projection[0] = 2 * camera.near / (xmax - xmin);
    right.projection[8] = (xmax + xmin) / (xmax - xmin);

    Mat4.mul(right.view, right.view, eyeRight);
}