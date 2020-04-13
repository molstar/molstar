/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Entities } from '../common';
import { CoarseElementData, CoarsedElementKeys } from '../coarse';
import { ChainIndex, ElementIndex, EntityIndex } from '../../indexing';
import SortedRanges from '../../../../../mol-data/int/sorted-ranges';
import { OrderedSet } from '../../../../../mol-data/int';

function getElementKey(map: Map<string, number>, key: string, counter: { index: number }) {
    if (map.has(key)) return map.get(key)!;
    const ret = counter.index++;
    map.set(key, ret);
    return ret;
}

function getElementSubstructureKeyMap(map: Map<number, Map<string, number>>, key: number) {
    if (map.has(key)) return map.get(key)!;
    const ret = new Map<string, number>();
    map.set(key, ret);
    return ret;
}

function createLookUp(entities: Entities, chain: Map<number, Map<string, number>>, seq: Map<number, SeqMap>) {
    const getEntKey = entities.getEntityIndex;
    const findChainKey: CoarsedElementKeys['findChainKey'] = (e, c) => {
        const eKey = getEntKey(e);
        if (eKey < 0) return -1 as ChainIndex;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1 as ChainIndex;
        return cm.get(c)! as ChainIndex;
    };
    const findSequenceKey: CoarsedElementKeys['findSequenceKey'] = (e, c, s) => {
        const eKey = getEntKey(e);
        if (eKey < 0) return -1 as ElementIndex;
        const cm = chain.get(eKey);
        if (cm === undefined) return -1 as ElementIndex;
        const cKey = cm.get(c);
        if (cKey === undefined) return -1 as ElementIndex;
        const sm = seq.get(cKey)!;
        const { elementIndices, seqRanges } = sm;
        const idx = SortedRanges.firstIntersectionIndex(seqRanges, OrderedSet.ofSingleton(s));
        return (idx !== -1 ? elementIndices[idx] : -1)  as ElementIndex;
    };
    return { findChainKey, findSequenceKey };
}

function missingEntity(k: string) {
    throw new Error(`Missing entity entry for entity id '${k}'.`);
}

type SeqMap = { elementIndices: number[], seqRanges: SortedRanges }

export function getCoarseKeys(data: CoarseElementData, entities: Entities): CoarsedElementKeys {
    const { entity_id, asym_id, seq_id_begin, seq_id_end, count, chainElementSegments } = data;

    const seqMaps = new Map<number, SeqMap>();
    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };

    const chainKey = new Int32Array(count) as any as ChainIndex[];
    const entityKey = new Int32Array(count) as any as EntityIndex[];

    const chainToEntity = new Int32Array(chainElementSegments.count) as any as EntityIndex[];

    for (let i = 0; i < count; i++) {
        entityKey[i] = entities.getEntityIndex(entity_id.value(i));
        if (entityKey[i] < 0) missingEntity(entity_id.value(i));
    }

    for (let cI = 0; cI < chainElementSegments.count; cI++) {
        const start = chainElementSegments.offsets[cI];
        const end = chainElementSegments.offsets[cI + 1];
        const eK = entityKey[start];

        chainToEntity[cI] = eK;

        const map = getElementSubstructureKeyMap(chainMaps, eK);
        const key = getElementKey(map, asym_id.value(start), chainCounter) as ChainIndex;
        for (let i = start; i < end; i++) chainKey[i] = key;

        // create seq_id map for the ranges defined by seq_id_begin and seq_id_end
        const elementIndices: number[] = [];
        const seqRanges: number[] = [];
        for (let i = start; i < end; i++) {
            const seqStart = seq_id_begin.value(i);
            const seqEnd = seq_id_end.value(i);
            elementIndices.push(i);
            seqRanges.push(seqStart, seqEnd);
        }
        const seqMap = { elementIndices, seqRanges: SortedRanges.ofSortedRanges(seqRanges) };
        seqMaps.set(key, seqMap);
    }

    const { findChainKey, findSequenceKey } = createLookUp(entities, chainMaps, seqMaps);

    const getEntityFromChain: CoarsedElementKeys['getEntityFromChain'] = c => {
        return chainToEntity[c];
    };

    return { chainKey, entityKey, findSequenceKey, findChainKey, getEntityFromChain };
}