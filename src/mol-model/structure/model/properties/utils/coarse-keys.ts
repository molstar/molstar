/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Entities } from '../common';
import { CoarseElementData, CoarsedElementKeys } from '../coarse';
import { ChainIndex, ElementIndex } from '../../indexing';

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

function createLookUp(entities: Entities, chain: Map<number, Map<string, number>>, seq: Map<number, Map<number, number>>) {
    const getEntKey = entities.getEntityIndex;
    const findChainKey: CoarsedElementKeys['findChainKey'] = (e, c) => {
        const eKey = getEntKey(e);
        if (eKey < 0) return -1 as ChainIndex;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1 as ChainIndex;
        return cm.get(c)! as ChainIndex;
    }
    // TODO consider implementing as binary search
    const findSequenceKey: CoarsedElementKeys['findSequenceKey'] = (e, c, s) => {
        const eKey = getEntKey(e);
        if (eKey < 0) return -1 as ElementIndex;
        const cm = chain.get(eKey);
        if (cm === undefined) return -1 as ElementIndex
        const cKey = cm.get(c)
        if (cKey === undefined) return -1 as ElementIndex
        const sm = seq.get(cKey)!
        if (!sm.has(s)) return -1 as ElementIndex;
        return sm.get(s)! as ElementIndex
    }
    return { findChainKey, findSequenceKey };
}

function missingEntity(k: string) {
    throw new Error(`Missing entity entry for entity id '${k}'.`);
}

export function getCoarseKeys(data: CoarseElementData, entities: Entities): CoarsedElementKeys {
    const { entity_id, asym_id, seq_id_begin, seq_id_end, count, chainSegments } = data;

    const seqMaps = new Map<number, Map<number, number>>();
    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };

    const chainKey = new Int32Array(count);
    const entityKey = new Int32Array(count);

    for (let i = 0; i < count; i++) {
        entityKey[i] = entities.getEntityIndex(entity_id.value(i));
        if (entityKey[i] < 0) missingEntity(entity_id.value(i));
    }

    for (let cI = 0; cI < chainSegments.count; cI++) {
        const start = chainSegments.offsets[cI], end = chainSegments.offsets[cI + 1];
        const map = getElementSubstructureKeyMap(chainMaps, entityKey[start]);
        const key = getElementKey(map, asym_id.value(start), chainCounter);
        for (let i = start; i < end; i++) chainKey[i] = key;

        // create seq_id map for the ranges defined by seq_id_begin and seq_id_end
        const seqMap: Map<number, number> = new Map()
        seqMaps.set(key, seqMap)
        for (let i = start; i < end; i++) {
            const seqStart = seq_id_begin.value(i)
            const seqEnd = seq_id_end.value(i)
            for (let j = seqStart; j <= seqEnd; j++) {
                seqMap.set(j, i)
            }
        }
    }



    const { findChainKey, findSequenceKey } = createLookUp(entities, chainMaps, seqMaps);

    return { chainKey, entityKey, findSequenceKey, findChainKey };
}