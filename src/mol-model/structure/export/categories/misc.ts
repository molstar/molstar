/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db';
import { mmCIF_Database, mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif';
import { CifWriter } from 'mol-io/writer/cif';
import { unionMany } from 'mol-util/set';
import { Model } from '../../model';
import { CifExportContext } from '../mmcif';
import CifCategory = CifWriter.Category
import { Structure } from '../../structure';

export const _chem_comp: CifCategory<CifExportContext> = {
    name: 'chem_comp',
    instance({ structures, cache }) {
        const chem_comp = getCifCategory(structures[0].model, 'chem_comp');
        if (!chem_comp) return CifCategory.Empty;
        const { id } = chem_comp;
        const names = cache.uniqueResidueNames || (cache.uniqueResidueNames = getUniqueResidueNames(structures));
        const indices = Column.indicesOf(id, id => names.has(id));
        return CifCategory.ofTable(chem_comp, indices);
    }
}

export const _pdbx_chem_comp_identifier: CifCategory<CifExportContext> = {
    name: 'pdbx_chem_comp_identifier',
    instance({ structures, cache }) {
        const pdbx_chem_comp_identifier = getCifCategory(structures[0].model, 'pdbx_chem_comp_identifier');
        if (!pdbx_chem_comp_identifier) return CifCategory.Empty;
        const { comp_id } = pdbx_chem_comp_identifier;
        const names = cache.uniqueResidueNames || (cache.uniqueResidueNames = getUniqueResidueNames(structures));
        const indices = Column.indicesOf(comp_id, id => names.has(id));
        return CifCategory.ofTable(pdbx_chem_comp_identifier, indices);
    }
}

function getCifCategory<K extends keyof mmCIF_Schema>(model: Model, name: K): mmCIF_Database[K] | undefined {
    if (model.sourceData.kind !== 'mmCIF') return;
    return model.sourceData.data[name];
}

function getUniqueResidueNames(structures: Structure[]) {
    return unionMany(structures.map(s => s.uniqueResidueNames));
}