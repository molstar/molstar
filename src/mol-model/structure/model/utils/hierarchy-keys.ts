/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { Data, Segments, Keys } from '../properties/hierarchy'
import { Interval, Segmentation } from 'mol-data/int'
import { Entities } from '../properties/common';

function getResidueId(comp_id: string, seq_id: number, ins_code: string) {
    return `${comp_id} ${seq_id} ${ins_code}`;
}

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

function createLookUp(entities: Entities, chain: Map<number, Map<string, number>>, residue: Map<number, Map<string, number>>) {
    const getEntKey = entities.getEntityIndex;
    const findChainKey: Keys['findChainKey'] = (e, c) => {
        let eKey = getEntKey(e);
        if (eKey < 0) return -1;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1;
        return cm.get(c)!;
    }
    const findResidueKey: Keys['findResidueKey'] = (e, c, name, seq, ins) => {
        let eKey = getEntKey(e);
        if (eKey < 0) return -1;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1;
        const rm = residue.get(cm.get(c)!)!
        const id = getResidueId(name, seq, ins);
        if (!rm.has(id)) return -1;
        return rm.get(id)!;
    }
    return { findChainKey, findResidueKey };
}

function checkMonotonous(xs: ArrayLike<number>) {
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i] < xs[i - 1]) {
            return false;
        }
    }
    return true;
}

function missingEntity(k: string) {
    throw new Error(`Missing entity entry for entity id '${k}'.`);
}

function create(data: Data, entities: Entities, segments: Segments): Keys {
    const { chains, residues } = data;

    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };
    const residueMaps = new Map<number, Map<string, number>>(), residueCounter = { index: 0 };

    const residueKey = new Int32Array(residues._rowCount);
    const chainKey = new Int32Array(chains._rowCount);
    const entityKey = new Int32Array(chains._rowCount);

    const { label_comp_id, auth_seq_id, pdbx_PDB_ins_code } = data.residues;
    const { label_entity_id, label_asym_id } = data.chains;

    const atomSet = Interval.ofBounds(0, data.atoms._rowCount);

    let isMonotonous = true;

    const chainsIt = Segmentation.transientSegments(segments.chainSegments, atomSet);
    while (chainsIt.hasNext) {
        const chainSegment = chainsIt.move();
        const cI = chainSegment.index;

        let eKey = entities.getEntityIndex(label_entity_id.value(cI));
        if (eKey < 0) missingEntity(label_entity_id.value(cI));
        const chainMap = getElementSubstructureKeyMap(chainMaps, eKey);
        const cKey = getElementKey(chainMap, label_asym_id.value(cI), chainCounter);

        chainKey[cI] = cKey;
        entityKey[cI] = eKey;

        const residueMap = getElementSubstructureKeyMap(residueMaps, cKey);
        const residuesIt = Segmentation.transientSegments(segments.residueSegments, atomSet, chainSegment);
        let last_seq_id = Number.NEGATIVE_INFINITY;
        while (residuesIt.hasNext) {
            const residueSegment = residuesIt.move();
            const rI = residueSegment.index;
            const seq_id = auth_seq_id.value(rI);
            if (seq_id < last_seq_id) isMonotonous = false;
            last_seq_id = seq_id;
            const residueId = getResidueId(label_comp_id.value(rI), auth_seq_id.value(rI), pdbx_PDB_ins_code.value(rI));
            residueKey[rI] = getElementKey(residueMap, residueId, residueCounter);
        }
    }

    const { findChainKey, findResidueKey } = createLookUp(entities, chainMaps, residueMaps);

    return {
        isMonotonous: isMonotonous && checkMonotonous(entityKey) && checkMonotonous(chainKey) && checkMonotonous(residueKey),
        residueKey: Column.ofIntArray(residueKey),
        chainKey: Column.ofIntArray(chainKey),
        entityKey: Column.ofIntArray(entityKey),
        findChainKey,
        findResidueKey
    };
}

export default create;