/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../mol-plugin/behavior/behavior';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BoundingSphereHelper, BoundingSphereHelperParams } from './bounding-sphere-helper';
import { ClipObjectHelper, ClipObjectHelperParams } from './clip-object-helper';
import { DirectVolumeHelper, DirectVolumeHelperParams } from './direct-volume-helper';
import { ImageHelper, ImageHelperParams } from './image-helper';
import { MeshHelper, MeshHelperParams } from './mesh-helper';
import { WebGPUDebugRegistry } from '../../mol-canvas3d/helper/debug-registry';
import { WebGPUHelperScene } from '../../mol-canvas3d/helper/webgpu-scene';

const DebugHelpersParams = {
    ...BoundingSphereHelperParams,
    ...ClipObjectHelperParams,
    ...MeshHelperParams,
    ...ImageHelperParams,
    ...DirectVolumeHelperParams,
};
type DebugHelpersParams = typeof DebugHelpersParams;
type DebugHelpersProps = PD.Values<DebugHelpersParams>;

export const DebugHelpers = PluginBehavior.create<DebugHelpersProps>({
    name: 'extension-debug-helpers',
    category: 'misc',
    display: {
        name: 'Debug Helpers'
    },
    ctor: class extends PluginBehavior.Handler<DebugHelpersProps> {
        async register(): Promise<void> {
            await this.ctx.canvas3dInitialized;
            this.registerHelpers();
            this.subscribeObservable(this.ctx.events.canvas3d.recreated, () => this.registerHelpers());
        }

        private registerHelpers() {
            const canvas3d = this.ctx.canvas3d;
            if (!canvas3d) return;

            const dr = canvas3d.debugRegistry;
            if (!dr) return;
            if (dr instanceof WebGPUDebugRegistry) {
                dr.register('bounding-sphere', new BoundingSphereHelper(undefined, dr.parent, this.params, new WebGPUHelperScene()));
                dr.register('clip-object', new ClipObjectHelper(undefined, dr.parent, this.params, new WebGPUHelperScene()));
                dr.register('mesh', new MeshHelper(undefined, dr.parent, this.params, new WebGPUHelperScene()));
                dr.register('image', new ImageHelper(undefined, dr.parent, this.params, new WebGPUHelperScene()));
                dr.register('direct-volume', new DirectVolumeHelper(undefined, dr.parent, this.params, new WebGPUHelperScene()));
                canvas3d.requestDraw();
                return;
            }
            const { ctx, parent } = dr;

            dr.register('bounding-sphere', new BoundingSphereHelper(ctx, parent, this.params));
            dr.register('clip-object', new ClipObjectHelper(ctx, parent, this.params));
            dr.register('mesh', new MeshHelper(ctx, parent, this.params));
            dr.register('image', new ImageHelper(ctx, parent, this.params));
            dr.register('direct-volume', new DirectVolumeHelper(ctx, parent, this.params));
            canvas3d.requestDraw();
        }

        update(params: DebugHelpersProps) {
            const changed = super.update(params);
            const canvas3d = this.ctx.canvas3d;
            if (changed && canvas3d) {
                canvas3d.debugRegistry?.setProps(params);
                canvas3d.requestDraw();
            }
            return changed;
        }

        unregister() {
            const canvas3d = this.ctx.canvas3d;
            if (!canvas3d) return;

            const dr = canvas3d.debugRegistry;
            if (!dr) return;
            dr.unregister('bounding-sphere');
            dr.unregister('clip-object');
            dr.unregister('mesh');
            dr.unregister('image');
            dr.unregister('direct-volume');
            canvas3d.requestDraw();
        }
    },
    params: () => DebugHelpersParams,
    canAutoUpdate: () => true,
});
