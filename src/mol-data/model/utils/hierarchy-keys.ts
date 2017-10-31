/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'
import { Data, Segments, Keys } from '../properties/hierarchy'
import Segmentation from '../../../mol-base/collections/integer/segmentation'
import Interval from '../../../mol-base/collections/integer/interval'

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

function createLookUp(entity: Map<string, number>, chain: Map<number, Map<string, number>>, residue: Map<number, Map<string, number>>) {
    const findEntityKey: Keys['findEntityKey'] = (id) => entity.has(id) ? entity.get(id)! : -1;
    const findChainKey: Keys['findChainKey'] = (e, c) => {
        if (!entity.has(e)) return -1;
        const cm = chain.get(entity.get(e)!)!;
        if (!cm.has(c)) return -1;
        return cm.get(c)!;
    }
    const findResidueKey: Keys['findResidueKey'] = (e, c, name, seq, ins) => {
        if (!entity.has(e)) return -1;
        const cm = chain.get(entity.get(e)!)!;
        if (!cm.has(c)) return -1;
        const rm = residue.get(cm.get(c)!)!
        const id = getResidueId(name, seq, ins);
        if (!rm.has(id)) return -1;
        return rm.get(id)!;
    }
    return { findEntityKey, findChainKey, findResidueKey };
}

function checkMonotonous(xs: ArrayLike<number>) {
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i] < xs[i - 1]) {
            return false;
        }
    }
    return true;
}

function create(data: Data, segments: Segments): Keys {
    const { chains, residues, entities } = data;

    const entityMap = Column.createFirstIndexMap(entities.id);
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

        const eKey = entityMap.get(label_entity_id.value(cI)) || 0;
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

    const { findEntityKey, findChainKey, findResidueKey } = createLookUp(entityMap, chainMaps, residueMaps);

    return {
        isMonotonous: isMonotonous && checkMonotonous(entityKey) && checkMonotonous(chainKey) && checkMonotonous(residueKey),
        residueKey: residueKey,
        chainKey: chainKey,
        entityKey: entityKey,
        findEntityKey,
        findChainKey,
        findResidueKey
    };
}

export default create;