/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Entities } from '../common';
import { CoarseElementData, CoarsedElementKeys } from '../coarse';

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

function createLookUp(entities: Entities, chain: Map<number, Map<string, number>>) {
    const getEntKey = entities.getEntityIndex;
    const findChainKey: CoarsedElementKeys['findChainKey'] = (e, c) => {
        let eKey = getEntKey(e);
        if (eKey < 0) return -1;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1;
        return cm.get(c)!;
    }
    return { findChainKey };
}

function missingEntity(k: string) {
    throw new Error(`Missing entity entry for entity id '${k}'.`);
}

function missingModel(k: string) {
    throw new Error(`Missing entity entry for model id '${k}'.`);
}

export function getCoarseKeys(data: CoarseElementData, modelIndex: (id: number) => number, entities: Entities): CoarsedElementKeys {
    const { model_id, entity_id, asym_id, count, chainSegments } = data;

    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };
    const chainKey = new Int32Array(count);
    const entityKey = new Int32Array(count);
    const modelKey = new Int32Array(count);

    for (let i = 0; i < count; i++) {
        entityKey[i] = entities.getEntityIndex(entity_id.value(i));
        if (entityKey[i] < 0) missingEntity(entity_id.value(i));
        modelKey[i] = modelIndex(model_id.value(i));
        if (modelKey[i] < 0) missingModel('' + model_id.value(i));
    }

    for (let cI = 0; cI < chainSegments.count; cI++) {
        const start = chainSegments.segments[cI], end = chainSegments.segments[cI + 1];
        const map = getElementSubstructureKeyMap(chainMaps, entityKey[start]);
        const key = getElementKey(map, asym_id.value(start), chainCounter);
        for (let i = start; i < end; i++) chainKey[i] = key;
    }

    const { findChainKey } = createLookUp(entities, chainMaps);

    return {
        chainKey: chainKey,
        entityKey: entityKey,
        modelKey: modelKey,
        findChainKey
    };
}