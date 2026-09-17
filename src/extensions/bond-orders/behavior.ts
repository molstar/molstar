/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { PluginBehavior } from '../../mol-plugin/behavior';
import { BondOrdersTrajectoryPreset } from './preset';
import './transforms';

export const BondOrders = PluginBehavior.create({
    name: 'bond-orders-prop',
    category: 'custom-props',
    display: {
        name: 'Bond Orders',
        description: 'Perceive missing bond orders and expose them through unit.bonds.'
    },
    ctor: class extends PluginBehavior.Handler {
        register() {
            this.ctx.builders.structure.hierarchy.registerPreset(BondOrdersTrajectoryPreset);
        }
        unregister() {
            this.ctx.builders.structure.hierarchy.unregisterPreset(BondOrdersTrajectoryPreset);
        }
    },
    params: () => ({})
});
