/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Scene } from '../../mol-gl/scene';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BoundingSphereHelper, BoundingSphereHelperParams } from './bounding-sphere-helper';
import { ClipObjectHelper, ClipObjectHelperParams } from './clip-object-helper';
import { DirectVolumeHelper, DirectVolumeHelperParams } from './direct-volume-helper';
import { ImageHelper, ImageHelperParams } from './image-helper';
import { MeshHelper, MeshHelperParams } from './mesh-helper';

export const DebugHelperParams = {
    ...BoundingSphereHelperParams,
    ...ClipObjectHelperParams,
    ...MeshHelperParams,
    ...ImageHelperParams,
    ...DirectVolumeHelperParams,
};
export type DebugHelperParams = typeof DebugHelperParams;
export type DebugHelperProps = PD.Values<DebugHelperParams>;

export class DebugHelper {
    readonly boundingSphereHelper: BoundingSphereHelper;
    readonly clipObjectHelper: ClipObjectHelper;
    readonly meshHelper: MeshHelper;
    readonly imageHelper: ImageHelper;
    readonly directVolumeHelper: DirectVolumeHelper;

    constructor(ctx: WebGLContext, parent: Scene, props: Partial<DebugHelperProps>) {
        this.boundingSphereHelper = new BoundingSphereHelper(ctx, parent, props);
        this.clipObjectHelper = new ClipObjectHelper(ctx, parent, props);
        this.meshHelper = new MeshHelper(ctx, parent, props);
        this.imageHelper = new ImageHelper(ctx, parent, props);
        this.directVolumeHelper = new DirectVolumeHelper(ctx, parent, props);
    }

    get boundingSphereScene() { return this.boundingSphereHelper.scene; }
    get clipScene() { return this.clipObjectHelper.scene; }
    get meshScene() { return this.meshHelper.scene; }
    get imageScene() { return this.imageHelper.scene; }
    get directVolumeScene() { return this.directVolumeHelper.scene; }

    update() {
        if (this.boundingSphereHelper.isEnabled) this.boundingSphereHelper.update();
        if (this.clipObjectHelper.isEnabled) this.clipObjectHelper.update();
        if (this.meshHelper.isEnabled) this.meshHelper.update();
        if (this.imageHelper.isEnabled) this.imageHelper.update();
        if (this.directVolumeHelper.isEnabled) this.directVolumeHelper.update();
    }

    syncVisibility() {
        this.boundingSphereHelper.syncVisibility();
        this.clipObjectHelper.syncVisibility();
        this.meshHelper.syncVisibility();
        this.imageHelper.syncVisibility();
        this.directVolumeHelper.syncVisibility();
    }

    clear() {
        this.boundingSphereHelper.clear();
        this.clipObjectHelper.clear();
        this.meshHelper.clear();
        this.imageHelper.clear();
        this.directVolumeHelper.clear();
    }

    get isEnabled() {
        return this.boundingSphereHelper.isEnabled || this.clipObjectHelper.isEnabled || this.meshHelper.isEnabled || this.imageHelper.isEnabled || this.directVolumeHelper.isEnabled;
    }

    get props(): Readonly<DebugHelperProps> {
        return {
            ...this.boundingSphereHelper.props,
            ...this.clipObjectHelper.props,
            ...this.meshHelper.props,
            ...this.imageHelper.props,
            ...this.directVolumeHelper.props,
        };
    }

    setProps(props: Partial<DebugHelperProps>) {
        this.boundingSphereHelper.setProps(props);
        this.clipObjectHelper.setProps(props);
        this.meshHelper.setProps(props);
        this.imageHelper.setProps(props);
        this.directVolumeHelper.setProps(props);
    }
}
