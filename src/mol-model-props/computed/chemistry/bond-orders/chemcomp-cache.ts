/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { Model, Unit } from '../../../../mol-model/structure';
import { BondType } from '../../../../mol-model/structure/model/types';
import { UnitIndex } from '../../../../mol-model/structure/structure/element/element';
import { IntraUnitBonds } from '../../../../mol-model/structure/structure/unit/bonds';
import { eachIntraBondedAtom } from '../util';
import { State, getFlags, getOrder, isPerceivable, setBond } from './common';

type BondingPattern = Map<string, { order: number, flags: number }>;

const PerceivedCache = new WeakMap<Model, Map<string, BondingPattern>>();

function getModelCache(model: Model) {
    let c = PerceivedCache.get(model);
    if (!c) { c = new Map(); PerceivedCache.set(model, c); }
    return c;
}

function getCompHash(unit: Unit.Atomic, start: number, end: number, compId: string) {
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const names: string[] = [], allNames = new Set<string>();
    for (let i = start; i < end; i++) {
        const name = label_atom_id.value(unit.elements[i]);
        if (allNames.has(name)) continue;
        names.push(name);
        allNames.add(name);
    }
    names.sort();
    return `${compId}|${names.join(',')}`;
}

function extractPattern(state: State): BondingPattern {
    const { unit, unitIndices, n, bonds } = state;
    const { label_atom_id } = unit.model.atomicHierarchy.atoms;
    const pattern: BondingPattern = new Map();
    for (let i = 0; i < n; i++) {
        const u = unitIndices[i];
        for (const j of state.heavyNeighbours[i]) {
            if (i >= j) continue;
            const v = unitIndices[j];
            const order = getOrder(bonds, u, v);
            const flags = getFlags(bonds, u, v);
            if (order > 1 || (flags & BondType.Flag.AromaticHuckel)) {
                const a = label_atom_id.value(unit.elements[u]);
                const b = label_atom_id.value(unit.elements[v]);
                const key = a < b ? a + '|' + b : b + '|' + a;
                pattern.set(key, { order, flags: flags & (BondType.Flag.AromaticHuckel) });
            }
        }
    }
    return pattern;
}

export function applyCachedChemCompPattern(unit: Unit.Atomic, bonds: IntraUnitBonds, start: UnitIndex, end: UnitIndex, cachePrefix = ''): boolean {
    const model = unit.model;
    const cache = getModelCache(model);
    const { label_comp_id, label_atom_id } = model.atomicHierarchy.atoms;
    const compId = label_comp_id.value(unit.elements[start]);
    const compHash = cachePrefix + getCompHash(unit, start, end, compId);
    const pattern = cache.get(compHash);
    if (!pattern || pattern.size === 0) return false;

    for (let u = start; u < end; u++) {
        const nameA = label_atom_id.value(unit.elements[u]);
        eachIntraBondedAtom(unit, u, (otherUnit, v) => {
            if (otherUnit.id !== unit.id || u > v) return;
            const nameB = label_atom_id.value(unit.elements[v]);
            const key = nameA < nameB ? nameA + '|' + nameB : nameB + '|' + nameA;
            const p = pattern.get(key);
            if (!p) return;
            const flags = getFlags(bonds, u, v);
            if (!isPerceivable(flags, getOrder(bonds, u, v))) return;
            setBond(bonds, u, v, p.order, BondType.Flag.Computed | p.flags);
        });
    }
    return true;
}

export function cacheChemCompPattern(unit: Unit.Atomic, start: UnitIndex, end: UnitIndex, cachePrefix = '', state: State) {
    const model = unit.model;
    const cache = getModelCache(model);
    const { label_comp_id } = model.atomicHierarchy.atoms;
    const compId = label_comp_id.value(unit.elements[start]);
    const compHash = cachePrefix + getCompHash(unit, start, end, compId);
    if (cache.has(compHash)) return;
    const pattern = extractPattern(state);
    cache.set(compHash, pattern);
}
