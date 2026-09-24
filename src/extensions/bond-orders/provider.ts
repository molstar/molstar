/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { IndexPairBonds } from '../../mol-model-formats/structure/property/bonds/index-pair';
import { BondProvider, BondProviderRegistry } from '../../mol-model/structure/structure/unit/bonds/bond-provider';
import { Model } from '../../mol-model/structure/model/model';
import { hasIntraBondOrderFromTable } from '../../mol-model/structure/model/properties/atomic/bonds';
import { BondType, WaterNames } from '../../mol-model/structure/model/types';
import { Structure } from '../../mol-model/structure/structure/structure';
import { Unit } from '../../mol-model/structure/structure/unit';
import { DefaultBondComputationProps } from '../../mol-model/structure/structure/unit/bonds/common';
import { IntraUnitBonds } from '../../mol-model/structure/structure/unit/bonds/data';
import { ElementSetIntraBondCache } from '../../mol-model/structure/structure/unit/bonds/element-set-intra-bond-cache';
import { findBonds } from '../../mol-model/structure/structure/unit/bonds/intra-compute';
import { BondOrdersMode, perceiveIntra } from './perceiver';

export const BondOrderProviderName = 'bond-order-perception';

/**
 * Lazy bond-order customization registered on a model.
 *
 * The provider owns structure-specific context while the registry integration is
 * generic. A future implementation can replace `perceiveIntra` without changing
 * `unit.bonds` or the registry contract.
 */
export class BondOrderProvider implements BondProvider {
    readonly name = BondOrderProviderName;
    private readonly cache = new ElementSetIntraBondCache();
    private readonly computing = new WeakMap<Unit.Atomic, IntraUnitBonds>();

    constructor(readonly structure: Structure, readonly mode: BondOrdersMode) {
    }

    getBonds(unit: Unit.Atomic): IntraUnitBonds | undefined {
        if (IndexPairBonds.Provider.get(unit.model)) return undefined;

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
            const perceived = perceiveIntra(this.structure, unit, bonds, unit.rings, this.mode);
            this.cache.set(unit.elements, perceived);
            return perceived;
        } finally {
            this.computing.delete(unit);
        }
    }
}

export interface RegisteredBondOrderProvider {
    readonly model: Model
    readonly provider: BondOrderProvider
}

export function registerBondOrderProviders(structure: Structure, mode: BondOrdersMode): RegisteredBondOrderProvider[] {
    if (structure.models.length !== 1) return [];

    const registered: RegisteredBondOrderProvider[] = [];
    const model = structure.models[0];
    const registry = BondProviderRegistry.get(model);
    const provider = new BondOrderProvider(structure, mode);
    registry.add(provider);
    registry.select(provider.name);
    registered.push({ model, provider });
    return registered;
}

export function unregisterBondOrderProviders(registered: ReadonlyArray<RegisteredBondOrderProvider>) {
    for (const { model, provider } of registered) {
        BondProviderRegistry.get(model).remove(provider);
    }
}

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
            if (mode === 'force' || (order[i] === 1 && BondType.is(flags[i], BondType.Flag.Computed))) return true;
        }
    }
    return false;
}
