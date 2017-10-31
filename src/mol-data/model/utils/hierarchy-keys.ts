/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'
import { HierarchyData, HierarchySegmentation, HierarchyKeys } from '../properties/hierarchy'
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
    const findEntity: HierarchyKeys['findEntity'] = (id) => entity.has(id) ? entity.get(id)! : -1;
    const findChain: HierarchyKeys['findChain'] = (e, c) => {
        if (!entity.has(e)) return -1;
        const cm = chain.get(entity.get(e)!)!;
        if (!cm.has(c)) return -1;
        return cm.get(c)!;
    }
    const findResidue: HierarchyKeys['findResidue'] = (e, c, name, seq, ins) => {
        if (!entity.has(e)) return -1;
        const cm = chain.get(entity.get(e)!)!;
        if (!cm.has(c)) return -1;
        const rm = residue.get(cm.get(c)!)!
        const id = getResidueId(name, seq, ins);
        if (!rm.has(id)) return -1;
        return rm.get(id)!;
    }
    return { findEntity, findChain, findResidue };
}

function isMonotonous(xs: ArrayLike<number>) {
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i] < xs[i - 1]) return false;
    }
    return true;
}

function create(data: HierarchyData, segments: HierarchySegmentation): HierarchyKeys {
    const { chains, residues, entities } = data;

    const entityMap = Column.createFirstIndexMap(entities.id);
    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };
    const residueMaps = new Map<number, Map<string, number>>(), residueCounter = { index: 0 };

    const residueKey = new Int32Array(residues._rowCount);
    const chainKey = new Int32Array(chains._rowCount);
    const entityKey = new Int32Array(chains._rowCount);

    const { auth_comp_id, auth_seq_id, pdbx_PDB_ins_code } = data.residues;
    const { label_entity_id, auth_asym_id } = data.chains;

    const chainsIt = Segmentation.transientSegments(segments.chains, Interval.ofBounds(0, data.atoms._rowCount));
    while (chainsIt.hasNext) {
        const chainSegment = chainsIt.move();
        const residuesIt = Segmentation.transientSegments(segments.residues, Interval.ofBounds(chainSegment.start, chainSegment.end));
        const cI = chainSegment.index;

        const eKey = entityMap.get(label_entity_id.value(cI)) || 0;
        const chainMap = getElementSubstructureKeyMap(chainMaps, eKey);
        const cKey = getElementKey(chainMap, auth_asym_id.value(cI), chainCounter);

        chainKey[cI] = cKey;
        entityKey[cI] = eKey;

        const residueMap = getElementSubstructureKeyMap(residueMaps, cKey);
        while (residuesIt.hasNext) {
            const rI = residuesIt.move().index;
            const residueId = getResidueId(auth_comp_id.value(rI), auth_seq_id.value(rI), pdbx_PDB_ins_code.value(rI));
            residueKey[rI] = getElementKey(residueMap, residueId, residueCounter);
        }
    }

    const { findEntity, findChain, findResidue } = createLookUp(entityMap, chainMaps, residueMaps);

    return {
        isMonotonous: isMonotonous(entityKey) && isMonotonous(chainKey) && isMonotonous(residueKey),
        residue: residueKey,
        chain: chainKey,
        entity: entityKey,
        findEntity,
        findChain,
        findResidue
    };
}

export default create;