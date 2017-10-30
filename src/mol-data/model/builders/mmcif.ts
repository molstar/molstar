/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RawData } from '../formats'
import { Frame as mmCIF } from '../../../mol-io/reader/cif/schema/mmcif'
import Model from '../../model'
import Interval from '../../../mol-base/collections/integer/interval'
import Segmentation from '../../../mol-base/collections/integer/segmentation'
import uuId from '../../../mol-base/utils/uuid'

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
    const residues = [0], chains = [0], entities = [0];

    const { label_entity_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, label_comp_id } = data.atom_site;

    let offset = 1;
    for (let i = start + 1; i < end; i++) {
        const newEntity = !label_entity_id.areValuesEqual(i - 1, i);
        const newChain = newEntity || !auth_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i)
            || !label_comp_id.areValuesEqual(i - 1, i);

        if (newEntity) entities[entities.length] = offset;
        if (newResidue) residues[residues.length] = offset;
        if (newChain) chains[chains.length] = offset;
        offset++;
    }

    residues[residues.length] = offset;
    chains[chains.length] = offset;
    entities[entities.length] = offset;

    return {
        residues: Segmentation.create(residues),
        chains: Segmentation.create(chains),
        entities: Segmentation.create(entities)
    };
}

function createModel(raw: RawData, data: mmCIF, bounds: Interval): Model {
    const segments = segment(data, bounds);
    return {
        id: uuId(),
        sourceData: raw,
        model_num: 0, // TODO: fix
        //common: 0 as any,
        macromolecule: 0 as any,
        conformation: 0 as any,
        version: { data: 0, conformation: 0 },
        atomCount: Interval.size(bounds),
        segments
    };
}

function buildModels(data: mmCIF): ArrayLike<Model> {
    const raw: RawData = { source: 'mmCIF', data };
    const models: Model[] = [];
    const atomCount = data.atom_site._rowCount;

    if (atomCount === 0) return models;

    let modelStart = 0;
    while (modelStart < atomCount) {
        const bounds = findModelBounds(data, modelStart);
        const model = createModel(raw, data, bounds);
        models.push(model);
        modelStart = Interval.end(bounds);
    }
    return models;
}

export default buildModels;