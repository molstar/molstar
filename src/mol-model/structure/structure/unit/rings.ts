/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { computeRings, getFingerprint, createIndex } from './rings/compute'
import Unit from '../unit';
import StructureElement from '../element';

type UnitRing = ReadonlyArray<StructureElement.UnitIndex>

interface UnitRings {
    /** Each ring is specified as an array of indices in Unit.elements. */
    readonly all: ReadonlyArray<UnitRing>,
    readonly byFingerprint: ReadonlyMap<string, ReadonlyArray<UnitRings.Index>>,

    readonly index: {
        /** Maps atom index inside a Unit to the smallest ring index (an atom can be part of more than one ring) */
        readonly elementRingIndices: ReadonlyMap<StructureElement.UnitIndex, UnitRings.Index[]>,

        /** Maps UnitRings.Index to index to ringComponents */
        readonly ringComponentIndex: ReadonlyArray<UnitRings.ComponentIndex>,
        readonly ringComponents: ReadonlyArray<ReadonlyArray<UnitRings.Index>>
    }
}

namespace UnitRings {
    /** Index into UnitRings.all */
    export type Index = { readonly '@type': 'unit-ring-index' } & number
    export type ComponentIndex = { readonly '@type': 'unit-ring-component-index' } & number

    export function getRingFingerprint(unit: Unit.Atomic, ring: UnitRing) {
        const { elements } = unit;
        const { type_symbol } = unit.model.atomicHierarchy.atoms;

        const symbols: string[] = [];
        for (let i = 0, _i = ring.length; i < _i; i++) symbols[symbols.length] = type_symbol.value(elements[ring[i]]) as String as string;
        return getFingerprint(symbols);
    }

    export function create(unit: Unit.Atomic): UnitRings {
        const rings = computeRings(unit);

        let _byFingerprint: Map<string, Index[]> | undefined = void 0;
        let _index: UnitRings['index'] | undefined = void 0;
        return {
            all: rings,
            get byFingerprint() {
                if (_byFingerprint) return _byFingerprint;
                _byFingerprint = new Map();
                let idx = 0 as Index;
                for (const r of rings) {
                    const fp = getRingFingerprint(unit, r);
                    if (_byFingerprint.has(fp)) _byFingerprint.get(fp)!.push(idx);
                    else _byFingerprint.set(fp, [idx]);
                    idx++;
                }
                return _byFingerprint;
            },
            get index() {
                if (_index) return _index;
                _index = createIndex(rings);
                return _index;
            }
        };
    }
}

export { UnitRing, UnitRings }