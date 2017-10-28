/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RawData } from '../formats'
import mmCIF from '../../../mol-io/reader/cif/schema/mmcif'
import Model from '../../model'
import Interval from '../../../mol-base/collections/integer/interval'
import Segmentation from '../../../mol-base/collections/integer/segmentation'

function findModelBounds(data: mmCIF, startIndex: number) {
    const num = data.atom_site.pdbx_PDB_model_num;
    const atomCount = num.rowCount;
    if (!num.isDefined) return Interval.ofBounds(startIndex, atomCount);
    let endIndex = startIndex + 1;
    while (endIndex < atomCount && num.areValuesEqual(startIndex, endIndex)) endIndex++;
    return Interval.ofBounds(startIndex, endIndex);
}

function segment(data: mmCIF, bounds: Interval) {
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const residues = [start], chains = [start];

    const { label_entity_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code } = data.atom_site;

    for (let i = start + 1; i < end; i++) {
        const newEntity = !label_entity_id.areValuesEqual(i - 1, i);
        const newChain = newEntity || !auth_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i);

        if (newResidue) residues[residues.length] = i;
        if (newChain) chains[chains.length] = i;
    }

    residues[residues.length] = end;
    chains[chains.length] = end;

    return { residues: Segmentation.create(residues), chains: Segmentation.create(chains) };
}

function createModel(raw: RawData, data: mmCIF, bounds: Interval): Model {
    const segments = segment(data, bounds);
    return {
        data: raw,
        common: 0 as any,
        macromolecule: 0 as any,
        residues: segments.residues,
        chains: segments.chains
    };
}

function getModels(data: mmCIF): ArrayLike<Model> {
    const raw: RawData = { source: 'mmCIF', data };
    const models: Model[] = [];
    const atomCount = data.atom_site._rowCount;
    let modelStart = 0;
    while (modelStart < atomCount) {
        const bounds = findModelBounds(data, modelStart);
        const model = createModel(raw, data, bounds);
        models.push(model);
        modelStart = Interval.end(bounds);
    }
    return models;
}

export default getModels;

// function createStructure() {

// }

