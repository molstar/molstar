/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RawData } from '../formats'
import { Frame as mmCIF } from '../../../mol-io/reader/cif/schema/mmcif'
import Model from '../../model'
//import Column from '../../../mol-base/collections/column'
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

function segmentOffsets(data: mmCIF, bounds: Interval) {
    const start = Interval.start(bounds), end = Interval.end(bounds);
    const residues = [start], chains = [start];

    const { label_entity_id, auth_asym_id, auth_seq_id, pdbx_PDB_ins_code, label_comp_id } = data.atom_site;

    for (let i = start + 1; i < end; i++) {
        const newChain = !label_entity_id.areValuesEqual(i - 1, i) || !auth_asym_id.areValuesEqual(i - 1, i);
        const newResidue = newChain
            || !auth_seq_id.areValuesEqual(i - 1, i)
            || !pdbx_PDB_ins_code.areValuesEqual(i - 1, i)
            || !label_comp_id.areValuesEqual(i - 1, i);

        if (newResidue) residues[residues.length] = i;
        if (newChain) chains[chains.length] = i;
    }

    return { residues, chains };
}

function createModel(raw: RawData, data: mmCIF, bounds: Interval): Model {
    const segments = segmentOffsets(data, bounds);
    return {
        id: uuId(),
        sourceData: raw,
        model_num: 0, // TODO: fix
        //common: 0 as any,
        macromolecule: 0 as any,
        conformation: 0 as any,
        version: { data: 0, conformation: 0 },
        atomCount: Interval.size(bounds),
        segments: {
            residues: Segmentation.ofOffsets(segments.residues, bounds),
            chains: Segmentation.ofOffsets(segments.chains, bounds),
        }
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