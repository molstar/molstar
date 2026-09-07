/**
 * Copyright (c) 2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Structure, Unit } from '../../mol-model/structure';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { BondOrdersValue, BondOrdersMode, calcBondOrders } from './chemistry/bond-orders';

function getBondOrdersParams(_data?: Structure) {
    return {
        // `auto`: perceive only Computed (from distance, without dict) bonds.
        // `model`: no perception — use only the file/table orders already in the data.
        // `force`: re-perceive every bond, overriding even authoritative file/table orders.
        mode: PD.MappedStatic('auto', {
            'auto': PD.EmptyGroup({ label: 'Automatic' }),
            'model': PD.EmptyGroup({ label: 'Model' }),
            'force': PD.EmptyGroup({ label: 'Force' }),
        }, { options: [['auto', 'Automatic'], ['model', 'Model'], ['force', 'Force']] })
    };
}

export const BondOrdersParams = getBondOrdersParams();
export type BondOrdersParams = typeof BondOrdersParams
export type BondOrdersProps = PD.Values<BondOrdersParams>

export const BondOrderProvider: CustomStructureProperty.Provider<BondOrdersParams, BondOrdersValue> = CustomStructureProperty.createProvider({
    label: 'Bond Orders',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_bond_orders',
        // TODO `cifExport` and `symbol`
    }),
    type: 'root',
    defaultParams: BondOrdersParams,
    getParams: getBondOrdersParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, _data: Structure, props: Partial<BondOrdersProps>) => {
        const p = { ...PD.getDefaultValues(BondOrdersParams), ...props };
        return { value: calcBondOrders(p.mode.name as BondOrdersMode) };
    }
});

/**
 * Resolve a unit's per-edge bond order/flags: the perceived overrides when `BondOrderProvider` is
 * attached, else the file/table orders from the bond graph. Resolve once per unit (the returned
 * arrays are parallel to `unit.bonds` edge indexing), not per edge, to keep render loops cheap.
 */
export function getUnitBondProps(structure: Structure, unit: Unit.Atomic): { order: ArrayLike<number>, flags: ArrayLike<number> } {
    const ov = BondOrderProvider.get(structure).value?.getUnit(structure, unit);
    if (ov) return ov;
    const { order, flags } = unit.bonds.edgeProps;
    return { order, flags };
}
