/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { computeRings, getFingerprint, createIndex } from './rings/compute';
import Unit from '../unit';
import StructureElement from '../element';
import { SortedArray } from '../../../../mol-data/int';
import { ResidueIndex } from '../../model';
import { ElementSymbol, BondType } from '../../model/types';
import { Elements } from '../../model/properties/atomic/types';
import { getPositions } from '../../util';
import { PrincipalAxes } from '../../../../mol-math/linear-algebra/matrix/principal-axes';
import { Vec3 } from '../../../../mol-math/linear-algebra';

type UnitRing = SortedArray<StructureElement.UnitIndex>

class UnitRings {
    /** Each ring is specified as an array of indices in Unit.elements. */
    readonly all: ReadonlyArray<UnitRing>;

    private _byFingerprint?: ReadonlyMap<UnitRing.Fingerprint, ReadonlyArray<UnitRings.Index>>;
    private _index?: {
        readonly elementRingIndices: ReadonlyMap<StructureElement.UnitIndex, UnitRings.Index[]>,
        readonly elementAromaticRingIndices: ReadonlyMap<StructureElement.UnitIndex, UnitRings.Index[]>,
        readonly ringComponentIndex: ReadonlyArray<UnitRings.ComponentIndex>,
        readonly ringComponents: ReadonlyArray<ReadonlyArray<UnitRings.Index>>
    };
    private _aromaticRings?: ReadonlyArray<UnitRings.Index>

    private get index() {
        if (this._index) return this._index;
        this._index = createIndex(this.all, this.aromaticRings);
        return this._index;
    }

    get byFingerprint() {
        if (this._byFingerprint) return this._byFingerprint;
        this._byFingerprint = createByFingerprint(this.unit, this.all);
        return this._byFingerprint;
    }

    /** Maps atom index inside a Unit to ring indices (an atom can be part of more than one ring) */
    get elementRingIndices() {
        return this.index.elementRingIndices;
    }

    get elementAromaticRingIndices() {
        return this.index.elementAromaticRingIndices;
    }

    /** Maps UnitRings.Index to index to ringComponents */
    get ringComponentIndex() {
        return this.index.ringComponentIndex;
    }

    get ringComponents() {
        return this.index.ringComponents;
    }

    get aromaticRings() {
        if (this._aromaticRings) return this._aromaticRings;
        this._aromaticRings = getAromaticRings(this.unit, this.all);
        return this._aromaticRings;
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

    const AromaticRingElements = new Set([
        Elements.B, Elements.C, Elements.N, Elements.O,
        Elements.SI, Elements.P, Elements.S,
        Elements.GE, Elements.AS,
        Elements.SN, Elements.SB,
        Elements.BI
    ] as ElementSymbol[]);
    const AromaticRingPlanarityThreshold = 0.05;

    export function isAromatic(unit: Unit.Atomic, ring: UnitRing): boolean {
        const { elements, bonds: { b, offset, edgeProps: { flags } } } = unit;
        const { type_symbol, label_comp_id } = unit.model.atomicHierarchy.atoms;

        // ignore Proline (can be flat because of bad geometry)
        if (label_comp_id.value(unit.elements[ring[0]]) === 'PRO') return false;

        let aromaticBondCount = 0;
        let hasAromaticRingElement = false;

        for (let i = 0, il = ring.length; i < il; ++i) {
            const aI = ring[i];
            if (!hasAromaticRingElement && AromaticRingElements.has(type_symbol.value(elements[aI]))) {
                hasAromaticRingElement = true;
            }

            for (let j = offset[aI], jl = offset[aI + 1]; j < jl; ++j) {
                // comes e.g. from `chem_comp_bond.pdbx_aromatic_flag`
                if (BondType.is(BondType.Flag.Aromatic, flags[j])) {
                    if (SortedArray.has(ring, b[j])) aromaticBondCount += 1;

                }
            }
        }
        if (aromaticBondCount === 2 * ring.length) return true;
        if (!hasAromaticRingElement) return false;

        const ma = PrincipalAxes.calculateMomentsAxes(getPositions(unit, ring));
        return Vec3.magnitude(ma.dirC) < AromaticRingPlanarityThreshold;
    }

    /** Get the alternate location of the 1st non '' alt loc atom. */
    export function getAltId(unit: Unit.Atomic, ring: UnitRing) {
        const { label_alt_id } = unit.model.atomicHierarchy.atoms;
        const { elements } = unit;

        for (let i = 0, il = ring.length; i < il; ++i) {
            const eI = elements[ring[i]];
            const altId = label_alt_id.value(eI);
            if (altId) return altId;
        }

        return '';
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

        for (let fI = 0, _fI = fingerprints.length; fI < _fI; fI++) {
            const fp = fingerprints[fI];
            addSingleResidueRings(rings, fp, map);
        }
        return map;
    }
}

function createByFingerprint(unit: Unit.Atomic, rings: ReadonlyArray<UnitRing>) {
    const byFingerprint = new Map<UnitRing.Fingerprint, UnitRings.Index[]>();
    let idx = 0 as UnitRings.Index;
    for (let rI = 0, _rI = rings.length; rI < _rI; rI++) {
        const r = rings[rI];
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
    for (let rI = 0, _rI = byFp.length; rI < _rI; rI++) {
        const r = byFp[rI];
        const idx = ringResidueIdx(rings.unit, rings.all[r]);
        if (idx >= 0) {
            if (map.has(idx)) map.get(idx)!.push(r);
            else map.set(idx, [r]);
        }
    }
}

function getAromaticRings(unit: Unit.Atomic, rings: ReadonlyArray<UnitRing>): ReadonlyArray<UnitRings.Index> {
    const aromaticRings: UnitRings.Index[] = [];
    for (let i = 0 as UnitRings.Index, il = rings.length; i < il; ++i) {
        if (UnitRing.isAromatic(unit, rings[i])) aromaticRings.push(i);
    }
    return aromaticRings;
}

export { UnitRing, UnitRings };