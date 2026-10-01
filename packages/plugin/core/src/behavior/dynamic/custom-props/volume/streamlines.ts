/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StreamlinesProvider } from '@molstar/model/props/volume/streamlines';
import { StreamlinesRepresentationProvider } from '@molstar/graphics/props/volume/streamlines/representation';

export const Streamlines = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'streamlines-volume-prop',
    category: 'custom-props',
    display: { name: 'Streamlines' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = StreamlinesProvider;

        update(p: { autoAttach: boolean, showTooltip: boolean }) {
            const updated = (
                this.params.autoAttach !== p.autoAttach
            );
            this.params.autoAttach = p.autoAttach;
            this.ctx.customVolumeProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        register(): void {
            this.ctx.customVolumeProperties.register(this.provider, this.params.autoAttach);
            this.ctx.representation.volume.registry.add(StreamlinesRepresentationProvider);
        }

        unregister() {
            this.ctx.customVolumeProperties.unregister(this.provider.descriptor.name);
            this.ctx.representation.volume.registry.remove(StreamlinesRepresentationProvider);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});
