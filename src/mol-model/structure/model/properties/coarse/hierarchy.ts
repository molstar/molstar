/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../../../mol-data/db';
import { Segmentation } from '../../../../../mol-data/int';
import { ElementIndex, ChainIndex, EntityIndex } from '../../indexing';
import { SortedRanges } from '../../../../../mol-data/int/sorted-ranges';
import { EmptyCoarseIndex } from '../utils/coarse-index';

export interface CoarsedElementKeys {
    /** Assign a key to each element */
    chainKey: ArrayLike<ChainIndex>,
    /** Assign a key to each element, index to the Model.entities.data table */
    entityKey: ArrayLike<EntityIndex>,

    /** Find index of the residue/feature element where seq_id is included */
    findSequenceKey(entityId: string, asym_id: string, seq_id: number): ElementIndex
    findChainKey(entityId: string, asym_id: string): ChainIndex

    /** Returns index or -1 if not present. */
    getEntityFromChain(cI: ChainIndex): EntityIndex
}

export interface CoarseElementData {
    count: number,
    /**
     * The entity identifier corresponding to this coarse object.
     * In mmCIF files, this points to entity_poly_seq.entity_id in the ENTITY_POLY category.
     */
    entity_id: Column<string>,
    /**
     * An asym/strand identifier corresponding to this coarse object.
     * In mmCIF files, this points to struct_asym.id in the STRUCT_ASYM category
     */
    asym_id: Column<string>,
    /**
     * The leading sequence index corresponding to this coarse object.
     * In mmCIF files, this points to entity_poly_seq.num in the ENTITY_POLY category.
     */
    seq_id_begin: Column<number>,
    /**
     * The trailing sequence index corresponding to this coarse object.
     * In mmCIF files, this points to entity_poly_seq.num in the ENTITY_POLY category.
     */
    seq_id_end: Column<number>,

    chainElementSegments: Segmentation<ElementIndex, ChainIndex>,
}

export interface CoarseRanges {
    polymerRanges: SortedRanges<ElementIndex>
    gapRanges: SortedRanges<ElementIndex>
}

export type CoarseElements = CoarsedElementKeys & CoarseElementData & CoarseRanges

export interface CoarseHierarchy {
    isDefined: boolean,
    spheres: CoarseElements,
    gaussians: CoarseElements,
    index: CoarseIndex
}

const EmptyCoarseElements: CoarseElements = {
    chainKey: [],
    entityKey: [],
    findSequenceKey: () => -1 as ElementIndex,
    findChainKey: () => -1 as ChainIndex,
    getEntityFromChain: () => -1 as EntityIndex,

    count: 0,
    entity_id: Column.Undefined(0, Column.Schema.str),
    asym_id: Column.Undefined(0, Column.Schema.str),
    seq_id_begin: Column.Undefined(0, Column.Schema.int),
    seq_id_end: Column.Undefined(0, Column.Schema.int),
    chainElementSegments: Segmentation.create([]),

    polymerRanges: SortedRanges.ofSortedRanges([]),
    gapRanges: SortedRanges.ofSortedRanges([]),
};

export interface CoarseIndex {
    /**
     * Find element index of a sphere
     * @param key
     * @returns index or -1 if the atom is not present.
     */
    findSphereElement(key: CoarseElementKey): ElementIndex

    /**
     * Find element index of a gaussian
     * @param key
     * @returns index or -1 if the atom is not present.
     */
    findGaussianElement(key: CoarseElementKey): ElementIndex

    /**
     * Finds coarse element and assigns a reference to it.
     * @param key
     */
    findElement(key: CoarseElementKey, out: CoarseElementReference): boolean
}

export interface CoarseElementReference { kind?: 'spheres' | 'gaussians', index: ElementIndex }
export function CoarseElementReference(): CoarseElementReference { return { kind: undefined, index: -1 as ElementIndex }; }

export interface CoarseElementKey { label_entity_id: string, label_asym_id: string, label_seq_id: number }
export function CoarseElementKey(): CoarseElementKey { return { label_entity_id: '', label_asym_id: '', label_seq_id: -1 }; }

export namespace CoarseHierarchy {
    export const Empty: CoarseHierarchy = {
        isDefined: false,
        spheres: EmptyCoarseElements,
        gaussians: EmptyCoarseElements,
        index: EmptyCoarseIndex,
    };
}