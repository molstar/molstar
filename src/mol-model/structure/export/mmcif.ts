/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from '../../../mol-io/writer/cif'
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif'
import { Structure } from '../structure'
import { _atom_site } from './categories/atom_site';
import CifCategory = CifWriter.Category
import { _struct_conf, _struct_sheet_range } from './categories/secondary-structure';
import { _chem_comp, _pdbx_chem_comp_identifier, _pdbx_nonpoly_scheme } from './categories/misc';
import { Model } from '../model';
import { getUniqueEntityIndicesFromStructures, copy_mmCif_category } from './categories/utils';
import { _struct_asym, _entity_poly, _entity_poly_seq } from './categories/sequence';
import { CustomPropertyDescriptor } from '../common/custom-property';
import { atom_site_operator_mapping } from './categories/atom_site_operator_mapping';

export interface CifExportContext {
    structures: Structure[],
    firstModel: Model,
    cache: any
}

export namespace CifExportContext {
    export function create(structures: Structure | Structure[]): CifExportContext {
        const structureArray = Array.isArray(structures) ? structures : [structures];
        return {
            structures: structureArray,
            firstModel: structureArray[0].model,
            cache: Object.create(null)
        };
    }
}

const _entity: CifCategory<CifExportContext> = {
    name: 'entity',
    instance({ structures }) {
        const indices = getUniqueEntityIndicesFromStructures(structures);
        return CifCategory.ofTable(structures[0].model.entities.data, indices);
    }
}

function isWithoutSymmetry(structure: Structure) {
    return structure.units.every(u => u.conformation.operator.isIdentity)
}

const Categories = [
    // Basics
    copy_mmCif_category('entry'),
    copy_mmCif_category('exptl'),
    _entity,

    // Symmetry
    copy_mmCif_category('cell', isWithoutSymmetry),
    copy_mmCif_category('symmetry', isWithoutSymmetry),

    // Assemblies
    copy_mmCif_category('pdbx_struct_assembly', isWithoutSymmetry),
    copy_mmCif_category('pdbx_struct_assembly_gen', isWithoutSymmetry),
    copy_mmCif_category('pdbx_struct_oper_list', isWithoutSymmetry),

    // Secondary structure
    _struct_conf,
    _struct_sheet_range,

    // Sequence
    _struct_asym,
    _entity_poly,
    _entity_poly_seq,

    // Branch
    copy_mmCif_category('pdbx_entity_branch'),
    copy_mmCif_category('pdbx_entity_branch_link'),
    copy_mmCif_category('pdbx_branch_scheme'),

    // Misc
    // TODO: filter for actual present residues?
    _chem_comp,
    _pdbx_chem_comp_identifier,
    copy_mmCif_category('atom_sites'),

    _pdbx_nonpoly_scheme,

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

function encodeCustomProp(customProp: CustomPropertyDescriptor, ctx: CifExportContext, encoder: CifWriter.Encoder, params: encode_mmCIF_categories_Params) {
    if (!customProp.cifExport || customProp.cifExport.categories.length === 0) return;

    const prefix = customProp.cifExport.prefix;
    const cats = customProp.cifExport.categories;

    let propCtx = ctx;
    if (customProp.cifExport.context) {
        const propId = CustomPropertyDescriptor.getUUID(customProp);
        if (ctx.cache[propId + '__ctx']) propCtx = ctx.cache[propId + '__ctx'];
        else {
            propCtx = customProp.cifExport.context(ctx) || ctx;
            ctx.cache[propId + '__ctx'] = propCtx;
        }
    }
    for (const cat of cats) {
        if (params.skipCategoryNames && params.skipCategoryNames.has(cat.name)) continue;
        if (cat.name.indexOf(prefix) !== 0) throw new Error(`Custom category '${cat.name}' name must start with prefix '${prefix}.'`);
        encoder.writeCategory(cat, propCtx);
    }
}

type encode_mmCIF_categories_Params = { skipCategoryNames?: Set<string>, exportCtx?: CifExportContext }

/** Doesn't start a data block */
export function encode_mmCIF_categories(encoder: CifWriter.Encoder, structures: Structure | Structure[], params?: encode_mmCIF_categories_Params) {
    const first = Array.isArray(structures) ? structures[0] : (structures as Structure);
    const models = first.models;
    if (models.length !== 1) throw 'Can\'t export stucture composed from multiple models.';

    const _params = params || { };
    const ctx: CifExportContext = params && params.exportCtx ? params.exportCtx : CifExportContext.create(structures);

    for (const cat of Categories) {
        if (_params.skipCategoryNames && _params.skipCategoryNames.has(cat.name)) continue;
        encoder.writeCategory(cat, ctx);
    }

    if ((!_params.skipCategoryNames || !_params.skipCategoryNames.has('atom_site')) && encoder.isCategoryIncluded('atom_site')) {
        atom_site_operator_mapping(encoder, ctx);
    }

    for (const customProp of models[0].customProperties.all) {
        encodeCustomProp(customProp, ctx, encoder, _params);
    }

    const structureCustomProps = new Set<CustomPropertyDescriptor>();
    for (const s of ctx.structures) {
        if (!s.hasCustomProperties) continue;
        for (const p of s.customPropertyDescriptors.all) structureCustomProps.add(p);
    }
    structureCustomProps.forEach(customProp => encodeCustomProp(customProp, ctx, encoder, _params));
}

function to_mmCIF(name: string, structure: Structure, asBinary = false) {
    const enc = CifWriter.createEncoder({ binary: asBinary });
    enc.startDataBlock(name);
    encode_mmCIF_categories(enc, structure);
    return enc.getData();
}

export default to_mmCIF