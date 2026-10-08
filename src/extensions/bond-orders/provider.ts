/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { IndexPairBonds } from '../../mol-model-formats/structure/property/bonds/index-pair';
import { BondProvider } from '../../mol-model/structure/structure/unit/bonds/bond-provider';
import { Model } from '../../mol-model/structure/model/model';
import { hasIntraBondOrderFromTable } from '../../mol-model/structure/model/properties/atomic/bonds';
import { BondType, WaterNames } from '../../mol-model/structure/model/types';
import { Unit } from '../../mol-model/structure/structure/unit';
import { DefaultBondComputationProps } from '../../mol-model/structure/structure/unit/bonds/common';
import { IntraUnitBonds } from '../../mol-model/structure/structure/unit/bonds/data';
import { ElementSetIntraBondCache } from '../../mol-model/structure/structure/unit/bonds/element-set-intra-bond-cache';
import { findBonds } from '../../mol-model/structure/structure/unit/bonds/intra-compute';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { BondOrdersMode, perceiveIntra } from './perceiver';

export const BondOrderProviderName = 'bond-order-perception';

export const BondOrderProviderParams = {
    mode: PD.Select<BondOrdersMode>('model', [
        ['model', 'Model'],
        ['auto', 'Auto'],
        ['forceCompute', 'Force Compute'],
    ]),
};

class BondOrderProviderInstance implements BondProvider {
    private readonly cache = new ElementSetIntraBondCache();
    private readonly computing = new WeakMap<Unit.Atomic, IntraUnitBonds>();

    constructor(readonly context: BondProvider.Context, readonly model: Model, readonly mode: BondOrdersMode) {
    }

    getBonds(unit: Unit.Atomic): IntraUnitBonds | undefined {
        if (unit.model !== this.model || this.mode === 'model') return undefined;
        if (IndexPairBonds.Provider.get(unit.model)) return undefined;

        const structure = this.context.structure;
        if (!structure) return undefined;

        const inProgress = this.computing.get(unit);
        if (inProgress) return inProgress;

        const cached = this.cache.get(unit.elements);
        if (cached) return cached;

        const bonds = unit.elements.length <= 1
            ? IntraUnitBonds.Empty
            : findBonds(unit, DefaultBondComputationProps);
        if (!hasPerceivableBond(unit, bonds, this.mode)) {
            this.cache.set(unit.elements, bonds);
            return bonds;
        }

        this.computing.set(unit, bonds);
        try {
            const perceived = perceiveIntra(structure, unit, bonds, unit.rings, this.mode);
            this.cache.set(unit.elements, perceived);
            return perceived;
        } finally {
            this.computing.delete(unit);
        }
    }
}

export const BondOrderProvider: BondProvider.Provider<typeof BondOrderProviderParams> = {
    name: BondOrderProviderName,
    label: 'Bond Order Perception',
    getParams: () => BondOrderProviderParams,
    isApplicable: model => !IndexPairBonds.Provider.get(model),
    factory: (model, props, context) => props.mode === 'model'
        ? undefined
        : new BondOrderProviderInstance(context, model, props.mode),
};

function hasPerceivableBond(unit: Unit.Atomic, bonds: IntraUnitBonds, mode: BondOrdersMode) {
    const { elements, residueIndex } = unit;
    const { type_symbol, label_comp_id } = unit.model.atomicHierarchy.atoms;
    const { offset, b, edgeProps: { order, flags } } = bonds;

    for (let u = 0, ul = elements.length; u < ul; u++) {
        const eU = elements[u];
        const compId = label_comp_id.value(eU);
        if (WaterNames.has(compId) || hasIntraBondOrderFromTable(compId) || type_symbol.value(eU) !== 'C') continue;

        for (let i = offset[u], il = offset[u + 1]; i < il; i++) {
            const v = b[i];
            if (u >= v || type_symbol.value(elements[v]) !== 'C') continue;
            if (residueIndex[elements[v]] !== residueIndex[eU] || !BondType.isCovalent(flags[i])) continue;
            if (mode === 'forceCompute' || (order[i] === 1 && BondType.is(flags[i], BondType.Flag.OrderUnknown))) return true;
        }
    }
    return false;
}
