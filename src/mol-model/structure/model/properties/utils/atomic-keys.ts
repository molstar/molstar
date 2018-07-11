/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AtomicData, AtomicSegments, AtomicKeys } from '../atomic'
import { Interval, Segmentation } from 'mol-data/int'
import { Entities } from '../common'

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
    const findChainKey: AtomicKeys['findChainKey'] = (e, c) => {
        let eKey = getEntKey(e);
        if (eKey < 0) return -1;
        const cm = chain.get(eKey)!;
        if (!cm.has(c)) return -1;
        return cm.get(c)!;
    }
    const findResidueKey: AtomicKeys['findResidueKey'] = (e, c, name, seq, ins) => {
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

function missingEntity(k: string) {
    throw new Error(`Missing entity entry for entity id '${k}'.`);
}

export function getAtomicKeys(data: AtomicData, entities: Entities, segments: AtomicSegments): AtomicKeys {
    const { chains, residues } = data;

    const chainMaps = new Map<number, Map<string, number>>(), chainCounter = { index: 0 };
    const residueMaps = new Map<number, Map<string, number>>(), residueCounter = { index: 0 };

    const residueKey = new Int32Array(residues._rowCount);
    const chainKey = new Int32Array(chains._rowCount);
    const entityKey = new Int32Array(chains._rowCount);

    const { label_comp_id, auth_seq_id, pdbx_PDB_ins_code } = data.residues;
    const { label_entity_id, label_asym_id } = data.chains;

    const atomSet = Interval.ofBounds(0, data.atoms._rowCount);

    const chainsIt = Segmentation.transientSegments(segments.chainAtomSegments, atomSet);
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
        const residuesIt = Segmentation.transientSegments(segments.residueAtomSegments, atomSet, chainSegment);
        while (residuesIt.hasNext) {
            const residueSegment = residuesIt.move();
            const rI = residueSegment.index;
            const residueId = getResidueId(label_comp_id.value(rI), auth_seq_id.value(rI), pdbx_PDB_ins_code.value(rI));
            residueKey[rI] = getElementKey(residueMap, residueId, residueCounter);
        }
    }

    const { findChainKey, findResidueKey } = createLookUp(entities, chainMaps, residueMaps);

    return { residueKey, chainKey, entityKey, findChainKey, findResidueKey };
}