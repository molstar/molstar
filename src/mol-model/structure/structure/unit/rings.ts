/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeRings, getFingerprint, createIndex } from './rings/compute'
import Unit from '../unit';
import StructureElement from '../element';
import { SortedArray } from '../../../../mol-data/int';
import { ResidueIndex } from '../../model';
import { ElementSymbol } from '../../model/types';

type UnitRing = SortedArray<StructureElement.UnitIndex>

class UnitRings {
    /** Each ring is specified as an array of indices in Unit.elements. */
    readonly all: ReadonlyArray<UnitRing>;

    private _byFingerprint?: ReadonlyMap<UnitRing.Fingerprint, ReadonlyArray<UnitRings.Index>>;
    private _index?: {
        readonly elementRingIndices: ReadonlyMap<StructureElement.UnitIndex, UnitRings.Index[]>,
        readonly ringComponentIndex: ReadonlyArray<UnitRings.ComponentIndex>,
        readonly ringComponents: ReadonlyArray<ReadonlyArray<UnitRings.Index>>
    };

    private get index() {
        if (this._index) return this._index;
        this._index = createIndex(this.all);
        return this._index;
    }

    get byFingerprint() {
        if (this._byFingerprint) return this._byFingerprint;
        this._byFingerprint = createByFingerprint(this.unit, this.all);
        return this._byFingerprint;
    }

    /** Maps atom index inside a Unit to the smallest ring index (an atom can be part of more than one ring) */
    get elementRingIndices() {
        return this.index.elementRingIndices;
    }

    /** Maps UnitRings.Index to index to ringComponents */
    get ringComponentIndex() {
        return this.index.ringComponentIndex;
    }

    get ringComponents() {
        return this.index.ringComponents;
    }

    constructor(all: ReadonlyArray<UnitRing>, public unit: Unit.Atomic) {
        this.all = all;
    }
}

namespace UnitRing {
    export type Fingerprint = { readonly '@type': 'unit-ring-fingerprint' } & string

    export function fingerprint(unit: Unit.Atomic, ring: UnitRing): Fingerprint {
        const { elements } = unit;
        const { type_symbol } = unit.model.atomicHierarchy.atoms;

        const symbols: ElementSymbol[] = [];
        for (let i = 0, _i = ring.length; i < _i; i++) symbols[symbols.length] = type_symbol.value(elements[ring[i]]);
        return elementFingerprint(symbols);
    }

    export function elementFingerprint(elements: ArrayLike<ElementSymbol>) {
        return getFingerprint(elements as ArrayLike<string> as string[]) as Fingerprint;
    }
}

namespace UnitRings {
    /** Index into UnitRings.all */
    export type Index = { readonly '@type': 'unit-ring-index' } & number
    export type ComponentIndex = { readonly '@type': 'unit-ring-component-index' } & number

    export function create(unit: Unit.Atomic): UnitRings {
        const rings = computeRings(unit);
        return new UnitRings(rings, unit);
    }

    /** Creates a mapping ResidueIndex -> list or rings that are on that residue and have one of the specified fingerprints. */
    export function byFingerprintAndResidue(rings: UnitRings, fingerprints: ReadonlyArray<UnitRing.Fingerprint>) {
        const map = new Map<ResidueIndex, Index[]>();
        for (const fp of fingerprints) {
            addSingleResidueRings(rings, fp, map);
        }
        return map;
    }
}

function createByFingerprint(unit: Unit.Atomic, rings: ReadonlyArray<UnitRing>) {
    const byFingerprint = new Map<UnitRing.Fingerprint, UnitRings.Index[]>();
    let idx = 0 as UnitRings.Index;
    for (const r of rings) {
        const fp = UnitRing.fingerprint(unit, r);
        if (byFingerprint.has(fp)) byFingerprint.get(fp)!.push(idx);
        else byFingerprint.set(fp, [idx]);
        idx++;
    }
    return byFingerprint;
}

function ringResidueIdx(unit: Unit.Atomic, ring: ArrayLike<StructureElement.UnitIndex>): ResidueIndex {
    const { elements } = unit;
    const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index;
    const idx = residueIndex[elements[ring[0]]];
    for (let rI = 1, _rI = ring.length; rI < _rI; rI++) {
        if (idx !== residueIndex[elements[ring[rI]]]) return -1 as ResidueIndex;
    }
    return idx;
}

function addSingleResidueRings(rings: UnitRings, fp: UnitRing.Fingerprint, map: Map<ResidueIndex, UnitRings.Index[]>) {
    const byFp = rings.byFingerprint.get(fp);
    if (!byFp) return;
    for (const r of byFp) {
        const idx = ringResidueIdx(rings.unit, rings.all[r]);
        if (idx >= 0) {
            if (map.has(idx)) map.get(idx)!.push(r);
            else map.set(idx, [r]);
        }
    }
}


export { UnitRing, UnitRings }