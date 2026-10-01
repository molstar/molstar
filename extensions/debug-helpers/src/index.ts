/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { BoundingSphereHelper, BoundingSphereHelperParams } from '@molstar/debug-helpers-extension/bounding-sphere-helper';
import { ClipObjectHelper, ClipObjectHelperParams } from '@molstar/debug-helpers-extension/clip-object-helper';
import { DirectVolumeHelper, DirectVolumeHelperParams } from '@molstar/debug-helpers-extension/direct-volume-helper';
import { ImageHelper, ImageHelperParams } from '@molstar/debug-helpers-extension/image-helper';
import { MeshHelper, MeshHelperParams } from '@molstar/debug-helpers-extension/mesh-helper';

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
            const canvas3d = this.ctx.canvas3d;
            if (!canvas3d) return;

            const dr = canvas3d.debugRegistry;
            const { ctx, parent } = dr;

            dr.register('bounding-sphere', new BoundingSphereHelper(ctx, parent, this.params));
            dr.register('clip-object', new ClipObjectHelper(ctx, parent, this.params));
            dr.register('mesh', new MeshHelper(ctx, parent, this.params));
            dr.register('image', new ImageHelper(ctx, parent, this.params));
            dr.register('direct-volume', new DirectVolumeHelper(ctx, parent, this.params));
        }

        update(params: DebugHelpersProps) {
            const changed = super.update(params);
            const canvas3d = this.ctx.canvas3d;
            if (changed && canvas3d) {
                canvas3d.debugRegistry.setProps(params);
                canvas3d.requestDraw();
            }
            return changed;
        }

        unregister() {
            const canvas3d = this.ctx.canvas3d;
            if (!canvas3d) return;

            const dr = canvas3d.debugRegistry;
            dr.unregister('bounding-sphere');
            dr.unregister('clip-object');
            dr.unregister('mesh');
            dr.unregister('image');
            dr.unregister('direct-volume');
        }
    },
    params: () => DebugHelpersParams,
    canAutoUpdate: () => true,
});
