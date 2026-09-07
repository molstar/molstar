/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { PluginBehavior } from '../../../behavior';
import { ParamDefinition as PD } from '../../../../../mol-util/param-definition';
import { BondOrderProvider } from '../../../../../mol-model-props/computed/bond-orders';

export const BondOrders = PluginBehavior.create<{ autoAttach: boolean }>({
    name: 'computed-bond-orders-prop',
    category: 'custom-props',
    display: { name: 'Bond Orders' },
    ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
        private provider = BondOrderProvider;

        update(p: { autoAttach: boolean }) {
            const updated = (
                this.params.autoAttach !== p.autoAttach
            );
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
