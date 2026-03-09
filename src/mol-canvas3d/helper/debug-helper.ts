/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BoundingSphereHelper, BoundingSphereHelperParams } from './bounding-sphere-helper';

export const DebugHelperParams = {
    ...BoundingSphereHelperParams,
};
export type DebugHelperParams = typeof DebugHelperParams;
export type DebugHelperProps = PD.Values<DebugHelperParams>;

export class DebugHelper {
    readonly boundingSphereHelper: BoundingSphereHelper;

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<DebugHelperProps>) {
        this.boundingSphereHelper = new BoundingSphereHelper(ctx, parent, props);
    }

    get boundingSphereScene() { return this.boundingSphereHelper.scene; }

    update() {
        if (this.boundingSphereHelper.isEnabled) this.boundingSphereHelper.update();
    }

    syncVisibility() {
        this.boundingSphereHelper.syncVisibility();
    }

    clear() {
        this.boundingSphereHelper.clear();
    }

    get isEnabled() {
        return this.boundingSphereHelper.isEnabled;
    }

    get props(): Readonly<DebugHelperProps> {
        return {
            ...this.boundingSphereHelper.props,
        };
    }

    setProps(props: Partial<DebugHelperProps>) {
        this.boundingSphereHelper.setProps(props);
    }
}
