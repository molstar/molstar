/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginBehavior } from '../../../behavior';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { InteractionsProvider } from '../../../../../mol-model-props/computed/interactions';

export const Interactions = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'computed-interactions-prop',
    category: 'custom-props',
    display: { name: 'Interactions' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = InteractionsProvider

        update(p: { autoAttach: boolean, showTooltip: boolean }) {
            let updated = (
                this.params.autoAttach !== p.autoAttach
            )
            this.params.autoAttach = p.autoAttach;
            this.ctx.customStructureProperties.setDefaultAutoAttach(this.provider.descriptor.name, this.params.autoAttach);
            return updated;
        }

        register(): void {
            this.ctx.customStructureProperties.register(this.provider, this.params.autoAttach);
        }

        unregister() {
            this.ctx.customStructureProperties.unregister(this.provider.descriptor.name);
        }
    },
    params: () => ({
        autoAttach: PD.Boolean(false)
    })
});