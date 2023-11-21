/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Column } from '../../../mol-data/db';
import { SortedArray } from '../../../mol-data/int';
import { ChainIndex, ElementIndex, Model, ResidueIndex } from '../../../mol-model/structure';
import { filterInPlace, range, sortIfNeeded } from '../../../mol-util/array';
import { Mapping, MultiMap, NumberMap } from './utils';


/** Auxiliary data structure for efficiently finding chains/residues/atoms in a model by their properties */
export interface IndicesAndSortings {
    chainsByLabelEntityId: Mapping<string, readonly ChainIndex[]>,
    chainsByLabelAsymId: Mapping<string, readonly ChainIndex[]>,
    chainsByAuthAsymId: Mapping<string, readonly ChainIndex[]>,
    residuesSortedByLabelSeqId: Mapping<ChainIndex, Sorting<ResidueIndex, number>>,
    residuesSortedByAuthSeqId: Mapping<ChainIndex, Sorting<ResidueIndex, number>>,
    residuesByInsCode: Mapping<ChainIndex, Mapping<string, readonly ResidueIndex[]>>,
    atomsById: Mapping<number, ElementIndex>,
    atomsByIndex: Mapping<number, ElementIndex>,
}

export const IndicesAndSortings = {
    /** Get `IndicesAndSortings` for a model (use a cached value or create if not available yet) */
    get(model: Model): IndicesAndSortings {
        return model._dynamicPropertyData['indices-and-sortings'] ??= IndicesAndSortings.create(model);
    },

    /** Create `IndicesAndSortings` for a model */
    create(model: Model): IndicesAndSortings {
        const h = model.atomicHierarchy;
        const nAtoms = h.atoms._rowCount;
        const nChains = h.chains._rowCount;
        const { label_entity_id, label_asym_id, auth_asym_id } = h.chains;
        const { label_seq_id, auth_seq_id, pdbx_PDB_ins_code } = h.residues;
        const { Present } = Column.ValueKind;

        const chainsByLabelEntityId = new MultiMap<string, ChainIndex>();
        const chainsByLabelAsymId = new MultiMap<string, ChainIndex>();
        const chainsByAuthAsymId = new MultiMap<string, ChainIndex>();
        const residuesSortedByLabelSeqId = new Map<ChainIndex, Sorting<ResidueIndex, number>>();
        const residuesSortedByAuthSeqId = new Map<ChainIndex, Sorting<ResidueIndex, number>>();
        const residuesByInsCode = new Map<ChainIndex, MultiMap<string, ResidueIndex>>();
        const atomsById = new NumberMap<number, ElementIndex>(nAtoms + 1);
        const atomsByIndex = new NumberMap<number, ElementIndex>(nAtoms);

        for (let iChain = 0 as ChainIndex; iChain < nChains; iChain++) {
            chainsByLabelEntityId.add(label_entity_id.value(iChain), iChain);
            chainsByLabelAsymId.add(label_asym_id.value(iChain), iChain);
            chainsByAuthAsymId.add(auth_asym_id.value(iChain), iChain);

            const iResFrom = h.residueAtomSegments.index[h.chainAtomSegments.offsets[iChain]];
            const iResTo = h.residueAtomSegments.index[h.chainAtomSegments.offsets[iChain + 1] - 1] + 1;

            const residuesWithLabelSeqId = filterInPlace(range(iResFrom, iResTo) as ResidueIndex[], iRes => label_seq_id.valueKind(iRes) === Present);
            residuesSortedByLabelSeqId.set(iChain, Sorting.create(residuesWithLabelSeqId, label_seq_id.value));

            const residuesWithAuthSeqId = filterInPlace(range(iResFrom, iResTo) as ResidueIndex[], iRes => auth_seq_id.valueKind(iRes) === Present);
            residuesSortedByAuthSeqId.set(iChain, Sorting.create(residuesWithAuthSeqId, auth_seq_id.value));

            const residuesHereByInsCode = new MultiMap<string, ResidueIndex>();
            for (let iRes = iResFrom; iRes < iResTo; iRes++) {
                if (pdbx_PDB_ins_code.valueKind(iRes) === Present) {
                    residuesHereByInsCode.add(pdbx_PDB_ins_code.value(iRes), iRes);
                }
            }
            residuesByInsCode.set(iChain, residuesHereByInsCode);
        }

        const atomId = model.atomicConformation.atomId.value;
        const atomIndex = h.atomSourceIndex.value;
        for (let iAtom = 0 as ElementIndex; iAtom < nAtoms; iAtom++) {
            atomsById.set(atomId(iAtom), iAtom);
            atomsByIndex.set(atomIndex(iAtom), iAtom);
        }

        return {
            chainsByLabelEntityId, chainsByLabelAsymId, chainsByAuthAsymId,
            residuesSortedByLabelSeqId, residuesSortedByAuthSeqId, residuesByInsCode,
            atomsById, atomsByIndex,
        };
    },
};


/** Represents a set of things (keys) of type `K`, sorted by some property (value) of type `V` */
export interface Sorting<K, V extends number> {
    /** Keys sorted by their corresponding values */
    keys: readonly K[],
    /** Sorted values corresponding to each key (value for `keys[i]` is `values[i]`) */
    values: SortedArray<V>,
}

export const Sorting = {
    /** Create a `Sorting` from an array of keys and a function returning their corresponding values.
     * If two keys have the same value, the smaller key will come first.
     * This function modifies `keys` - create a copy if you need the original order! */
    create<K extends number, V extends number>(keys: K[], valueFunction: (k: K) => V): Sorting<K, V> {
        sortIfNeeded(keys, (a, b) => valueFunction(a) - valueFunction(b) || a - b);
        const values: SortedArray<V> = SortedArray.ofSortedArray(keys.map(valueFunction));
        return { keys, values };
    },

    /** Return a newly allocated array of keys which have value equal to `target`.
     * The returned keys are sorted by their value. */
    getKeysWithValue<K, V extends number>(sorting: Sorting<K, V>, target: V): K[] {
        return Sorting.getKeysWithValueInRange(sorting, target, target);
    },

    /** Return a newly allocated array of keys which have value within interval `[min, max]` (inclusive).
     * The returned keys are sorted by their value.
     * Undefined `min` is interpreted as negative infitity, undefined `max` is interpreted as positive infinity. */
    getKeysWithValueInRange<K, V extends number>(sorting: Sorting<K, V>, min: V | undefined, max: V | undefined): K[] {
        const { keys, values } = sorting;
        if (!keys) return [];
        const n = keys.length;
        const from = (min !== undefined) ? SortedArray.findPredecessorIndex(values, min) : 0;
        let to: number;
        if (max !== undefined) {
            to = from;
            while (to < n && values[to] <= max) to++;
        } else {
            to = n;
        }
        return keys.slice(from, to);
    },
};
