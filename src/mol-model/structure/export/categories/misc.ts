/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from '../../../../mol-data/db';
import { CifWriter } from '../../../../mol-io/writer/cif';
import { CifExportContext } from '../mmcif';
import { getModelMmCifCategory, getUniqueResidueNamesFromStructures } from './utils';
import CifCategory = CifWriter.Category

export const _chem_comp: CifCategory<CifExportContext> = {
    name: 'chem_comp',
    instance({ firstModel, structures, cache }) {
        const chem_comp = getModelMmCifCategory(structures[0].model, 'chem_comp');
        if (!chem_comp) return CifCategory.Empty;
        const { id } = chem_comp;
        const names = cache.uniqueResidueNames || (cache.uniqueResidueNames = getUniqueResidueNamesFromStructures(structures));
        const indices = Column.indicesOf(id, id => names.has(id));
        return CifCategory.ofTable(chem_comp, indices);
    }
};

export const _pdbx_chem_comp_identifier: CifCategory<CifExportContext> = {
    name: 'pdbx_chem_comp_identifier',
    instance({ firstModel, structures, cache }) {
        const pdbx_chem_comp_identifier = getModelMmCifCategory(firstModel, 'pdbx_chem_comp_identifier');
        if (!pdbx_chem_comp_identifier) return CifCategory.Empty;
        const { comp_id } = pdbx_chem_comp_identifier;
        const names = cache.uniqueResidueNames || (cache.uniqueResidueNames = getUniqueResidueNamesFromStructures(structures));
        const indices = Column.indicesOf(comp_id, id => names.has(id));
        return CifCategory.ofTable(pdbx_chem_comp_identifier, indices);
    }
};

export const _pdbx_nonpoly_scheme: CifCategory<CifExportContext> = {
    name: 'pdbx_nonpoly_scheme',
    instance({ firstModel, structures, cache }) {
        const pdbx_nonpoly_scheme = getModelMmCifCategory(firstModel, 'pdbx_nonpoly_scheme');
        if (!pdbx_nonpoly_scheme) return CifCategory.Empty;
        // TODO: filter?
        return CifCategory.ofTable(pdbx_nonpoly_scheme);
    }
};