/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from 'mol-io/writer/cif'
import { mmCIF_Schema } from 'mol-io/reader/cif/schema/mmcif'
import { Structure } from '../structure'
import { _atom_site } from './categories/atom_site';
import CifCategory = CifWriter.Category
import { _struct_conf, _struct_sheet_range } from './categories/secondary-structure';
import { _pdbx_struct_mod_residue } from './categories/modified-residues';

export interface CifExportContext {
    structures: Structure[],
    cache: any
}

export namespace CifExportContext {
    export function create(structures: Structure | Structure[]): CifExportContext {
        return {
            structures: Array.isArray(structures) ? structures : [structures],
            cache: Object.create(null)
        };
    }
}

function copy_mmCif_category(name: keyof mmCIF_Schema): CifCategory<CifExportContext> {
    return {
        name,
        instance({ structures }) {
            const model = structures[0].model;
            if (model.sourceData.kind !== 'mmCIF') return CifCategory.Empty;
            const table = model.sourceData.data[name];
            if (!table || !table._rowCount) return CifCategory.Empty;
            return CifCategory.ofTable(table);
        }
    };
}

const _entity: CifCategory<CifExportContext> = {
    name: 'entity',
    instance({ structures }) {
        const indices = structures[0].entityIndices;
        return CifCategory.ofTable(structures[0].model.entities.data, indices);
    }
}

const Categories = [
    // Basics
    copy_mmCif_category('entry'),
    copy_mmCif_category('exptl'),
    _entity,

    // Symmetry
    copy_mmCif_category('cell'),
    copy_mmCif_category('symmetry'),

    // Assemblies
    copy_mmCif_category('pdbx_struct_assembly'),
    copy_mmCif_category('pdbx_struct_assembly_gen'),
    copy_mmCif_category('pdbx_struct_oper_list'),

    // Secondary structure
    _struct_conf,
    _struct_sheet_range,

    // Sequence
    copy_mmCif_category('struct_asym'), // TODO: filter only present chains?
    copy_mmCif_category('entity_poly'),
    copy_mmCif_category('entity_poly_seq'),

    // Branch
    copy_mmCif_category('pdbx_entity_branch'),
    copy_mmCif_category('pdbx_entity_branch_link'),
    copy_mmCif_category('pdbx_branch_scheme'),

    // Misc
    // TODO: filter for actual present residues?
    copy_mmCif_category('chem_comp'),
    copy_mmCif_category('pdbx_chem_comp_identifier'),
    copy_mmCif_category('atom_sites'),

    _pdbx_struct_mod_residue,

    // Atoms
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
export function encode_mmCIF_categories(encoder: CifWriter.Encoder, structures: Structure | Structure[], params?: { skipCategoryNames?: Set<string>, exportCtx?: CifExportContext }) {
    const first = Array.isArray(structures) ? structures[0] : (structures as Structure);
    const models = first.models;
    if (models.length !== 1) throw 'Can\'t export stucture composed from multiple models.';

    const _params = params || { };
    const ctx: CifExportContext = params && params.exportCtx ? params.exportCtx : CifExportContext.create(structures);

    for (const cat of Categories) {
        if (_params.skipCategoryNames && _params.skipCategoryNames.has(cat.name)) continue;
        encoder.writeCategory(cat, ctx);
    }

    for (const customProp of models[0].customProperties.all) {
        if (!customProp.cifExport || customProp.cifExport.categories.length === 0) continue;

        const prefix = customProp.cifExport.prefix;
        const cats = customProp.cifExport.categories;
        for (const cat of cats) {
            if (_params.skipCategoryNames && _params.skipCategoryNames.has(cat.name)) continue;
            if (cat.name.indexOf(prefix) !== 0) throw new Error(`Custom category '${cat.name}' name must start with prefix '${prefix}.'`);
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