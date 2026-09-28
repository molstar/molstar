/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { PluginStateObject as SO, PluginStateTransform } from '../../mol-plugin-state/objects';
import { Task } from '../../mol-task';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BondOrdersMode } from './perceiver';
import { RegisteredBondOrderProvider, registerBondOrderProviders, unregisterBondOrderProviders } from './provider';

export { PerceiveBondOrders };

interface PerceiveBondOrdersCache {
    registered: RegisteredBondOrderProvider[]
}

type PerceiveBondOrders = typeof PerceiveBondOrders
const PerceiveBondOrders = PluginStateTransform.BuiltIn({
    name: 'perceive-bond-orders',
    display: { name: 'Perceive Bond Orders', description: 'Perceive missing bond orders and expose them through unit.bonds.' },
    isDecorator: true,
    from: SO.Molecule.Structure,
    to: SO.Molecule.Structure,
    params: {
        mode: PD.Select<BondOrdersMode>('model', [
            ['model', 'Model (fill unknown orders)'],
            ['force', 'Force (re-perceive all)'],
        ]),
    },
})({
    apply({ a, params, cache }) {
        return Task.create('Perceive Bond Orders', async () => {
            const registered = registerBondOrderProviders(a.data, params.mode);
            (cache as PerceiveBondOrdersCache).registered = registered;
            return new SO.Molecule.Structure(a.data, { label: a.label, description: a.description });
        });
    },
    dispose({ cache }) {
        const registered = (cache as PerceiveBondOrdersCache | undefined)?.registered;
        if (registered) unregisterBondOrderProviders(registered);
    }
});
