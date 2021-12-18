/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BoundingSphereHelper, DebugHelperParams } from './bounding-sphere-helper';
import { CameraHelper, CameraHelperParams } from './camera-helper';
import { HandleHelper, HandleHelperParams } from './handle-helper';

export const HelperParams = {
    debug: PD.Group(DebugHelperParams),
    camera: PD.Group({
        helper: PD.Group(CameraHelperParams)
    }),
    handle: PD.Group(HandleHelperParams),
};
export const DefaultHelperProps = PD.getDefaultValues(HelperParams);
export type HelperProps = PD.Values<typeof HelperParams>


export class Helper {
    readonly debug: BoundingSphereHelper;
    readonly camera: CameraHelper;
    readonly handle: HandleHelper;

    constructor(webgl: WebGLContext, scene: Scene, props: Partial<HelperProps> = {}) {
        const p = { ...DefaultHelperProps, ...props };

        this.debug = new BoundingSphereHelper(webgl, scene, p.debug);
        this.camera = new CameraHelper(webgl, p.camera.helper);
        this.handle = new HandleHelper(webgl, p.handle);
    }
}