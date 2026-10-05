/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Scene } from '@molstar/graphics/gl/scene';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DebugRegistry } from './debug-registry.js';
import { CameraHelper, CameraHelperParams } from './camera-helper.js';
import { HandleHelper, HandleHelperParams } from './handle-helper.js';
import { PointerHelper, PointerHelperParams } from './pointer-helper.js';

export const HelperParams = {
    camera: PD.Group({
        helper: PD.Group(CameraHelperParams)
    }),
    handle: PD.Group(HandleHelperParams),
    pointer: PD.Group(PointerHelperParams),
};
export const DefaultHelperProps = PD.getDefaultValues(HelperParams);
export type HelperProps = PD.Values<typeof HelperParams>


export class Helper {
    readonly debug: DebugRegistry;
    readonly camera: CameraHelper;
    readonly handle: HandleHelper;
    readonly pointer: PointerHelper;

    constructor(webgl: WebGLContext, scene: Scene, props: Partial<HelperProps> = {}) {
        const p = { ...DefaultHelperProps, ...props };

        this.debug = new DebugRegistry(webgl, scene);
        this.camera = new CameraHelper(webgl, p.camera.helper);
        this.handle = new HandleHelper(webgl, p.handle);
        this.pointer = new PointerHelper(webgl, p.pointer);
    }
}