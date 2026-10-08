/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior';
import { createModelBondProviderProperty } from '../../mol-model-props/common/model-bond-provider';
import { BondOrdersTrajectoryPreset } from './preset';
import { BondOrderProvider } from './provider';

export const BondOrders = PluginBehavior.create({
    name: 'bond-orders-prop',
    category: 'custom-props',
    display: {
        name: 'Bond Orders',
        description: 'Perceive missing bond orders and expose them through unit.bonds.'
    },
    ctor: class extends PluginBehavior.Handler {
        private readonly property = createModelBondProviderProperty(this.ctx.model.bondProviderRegistry);

        register() {
            this.ctx.model.bondProviderRegistry.add(BondOrderProvider);
            this.ctx.customModelProperties.register(this.property, true);
            this.ctx.builders.structure.hierarchy.registerPreset(BondOrdersTrajectoryPreset);
        }
        unregister() {
            this.ctx.builders.structure.hierarchy.unregisterPreset(BondOrdersTrajectoryPreset);
            this.ctx.customModelProperties.unregister(this.property.descriptor.name);
            this.ctx.model.bondProviderRegistry.remove(BondOrderProvider);
        }
    },
    params: () => ({})
});
