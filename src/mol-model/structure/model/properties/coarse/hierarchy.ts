/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from 'mol-data/db'
import { Segmentation } from 'mol-data/int';
import { Element } from '../../../structure'
import SortedRanges from 'mol-data/int/sorted-ranges';

export interface CoarsedElementKeys {
    // assign a key to each element
    chainKey: ArrayLike<number>,
    // assign a key to each element, index to the Model.entities.data table
    entityKey: ArrayLike<number>,

    /** find index of the residue/feature element where seq_id is included */
    findSequenceKey(entityId: string, asym_id: string, seq_id: number): number
    findChainKey(entityId: string, asym_id: string): number
}

export interface CoarseElementData {
    count: number,
    entity_id: Column<string>,
    asym_id: Column<string>,
    seq_id_begin: Column<number>,
    seq_id_end: Column<number>,

    chainSegments: Segmentation<Element>,
}

export interface CoarseRanges {
    polymerRanges: SortedRanges<Element>
    gapRanges: SortedRanges<Element>
}

export type CoarseElements = CoarsedElementKeys & CoarseElementData & CoarseRanges

export interface CoarseHierarchy {
    isDefined: boolean,
    spheres: CoarseElements,
    gaussians: CoarseElements
}

export namespace CoarseHierarchy {
    export const Empty: CoarseHierarchy = { isDefined: false } as any;
}