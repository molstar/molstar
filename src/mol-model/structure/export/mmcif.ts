/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import { Structure } from '../structure'
import { Model } from '../model'
import { _atom_site } from './categories/atom_site';

export interface CifExportContext {
    structure: Structure,
    model: Model
}

import CifCategory = CifWriter.Category

function copy_mmCif_category(name: keyof mmCIF_Schema) {
    return ({ model }: CifExportContext) => {
        if (model.sourceData.kind !== 'mmCIF') return CifCategory.Empty;
        const table = model.sourceData.data[name];
        if (!table || !table._rowCount) return CifCategory.Empty;
        return CifCategory.ofTable(name, table);
    };
}

function _entity({ model, structure }: CifExportContext): CifCategory {
    const keys = Structure.getEntityKeys(structure);
    return CifCategory.ofTable('entity', model.entities.data, keys);
}

const Categories = [
    copy_mmCif_category('entry'),
    copy_mmCif_category('exptl'),
    _entity,
    copy_mmCif_category('cell'),
    copy_mmCif_category('symmetry'),
    copy_mmCif_category('pdbx_struct_assembly'),
    copy_mmCif_category('pdbx_struct_assembly_gen'),
    copy_mmCif_category('pdbx_struct_oper_list'),
    // TODO: filter for actual present residues?
    copy_mmCif_category('chem_comp'),
    copy_mmCif_category('atom_sites'),
    _atom_site
];

namespace _Filters {
    export const AtomSitePositionsFieldNames = new Set<string>(<(keyof typeof mmCIF_Schema.atom_site)[]>['id', 'Cartn_x', 'Cartn_y', 'Cartn_z']);
}

export const mmCIF_Export_Filters = {
    onlyPositions: <CifWriter.Category.Filter>{
        includeCategory(name) { return name === 'atom_site'; },
        includeField(cat, field) { return _Filters.AtomSitePositionsFieldNames.has(field); }
    }
}

/** Doesn't start a data block */
export function encode_mmCIF_categories(encoder: CifWriter.Encoder, structure: Structure) {
    const models = Structure.getModels(structure);
    if (models.length !== 1) throw 'Can\'t export stucture composed from multiple models.';
    const model = models[0];

    const ctx: CifExportContext[] = [{ structure, model }];

    for (const cat of Categories) {
        encoder.writeCategory(cat, ctx);
    }
    for (const customProp of model.customProperties.all) {
        const cats = customProp.cifExport.categoryProvider(ctx[0]);
        for (const cat of cats) {
            encoder.writeCategory(cat, ctx);
        }
    }
}

function to_mmCIF(name: string, structure: Structure, asBinary = false) {
    const enc = CifWriter.createEncoder({ binary: asBinary });
    enc.startDataBlock(name);
    encode_mmCIF_categories(enc, structure);
    return enc.getData();
}

export default to_mmCIF