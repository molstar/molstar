/**
 * Copyright (c) 2017-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CifWriter } from '../../../mol-io/writer/cif';
import { mmCIF_Schema } from '../../../mol-io/reader/cif/schema/mmcif';
import { Structure } from '../structure';
import { _atom_site } from './categories/atom_site';
import CifCategory = CifWriter.Category
import { _struct_conf, _struct_sheet_range } from './categories/secondary-structure';
import { _chem_comp, _chem_comp_bond, _pdbx_chem_comp_identifier, _pdbx_nonpoly_scheme } from './categories/misc';
import { Model } from '../model';
import { getUniqueEntityIndicesFromStructures, copy_mmCif_category, copy_source_mmCifCategory } from './categories/utils';
import { _struct_asym, _entity_poly, _entity_poly_seq } from './categories/sequence';
import { CustomPropertyDescriptor } from '../../custom-property';
import { atom_site_operator_mapping } from './categories/atom_site_operator_mapping';
import { MmcifFormat } from '../../../mol-model-formats/structure/mmcif';
import { molstar_bond_site } from './categories/molstar_bond_site';

export interface CifExportContext {
    structures: Structure[],
    firstModel: Model,
    cache: any
}

export type CifExportCategoryInfo =
    | [CifWriter.Category, any /** context */, CifWriter.Encoder.WriteCategoryOptions]
    | [CifWriter.Category, any /** context */]

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
};

function isWithoutSymmetry(structure: Structure) {
    return structure.units.every(u => u.conformation.operator.isIdentity);
}

function isWithoutOperator(structure: Structure) {
    return isWithoutSymmetry(structure) && structure.units.every(u => !u.conformation.operator.assembly && !u.conformation.operator.suffix);
}

const Categories = (options?: { keepAtomSiteId?: boolean }) => [
    // Basics
    copy_mmCif_category('entry'),
    copy_mmCif_category('exptl'),
    _entity,

    // Symmetry
    copy_mmCif_category('cell', isWithoutSymmetry),
    copy_mmCif_category('symmetry', isWithoutSymmetry),

    // Assemblies
    copy_mmCif_category('pdbx_struct_assembly', isWithoutOperator),
    copy_mmCif_category('pdbx_struct_assembly_gen', isWithoutOperator),
    copy_mmCif_category('pdbx_struct_oper_list', isWithoutOperator),

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

    // Struct conn
    copy_mmCif_category('struct_conn'),

    // Misc
    _chem_comp,
    _chem_comp_bond,
    _pdbx_chem_comp_identifier,
    copy_mmCif_category('atom_sites'),

    _pdbx_nonpoly_scheme,

    // Atoms
    _atom_site({ keepId: options?.keepAtomSiteId }),
];

namespace _Filters {
    export const AtomSitePositionsFieldNames = new Set<string>(<(keyof typeof mmCIF_Schema.atom_site)[]>['id', 'Cartn_x', 'Cartn_y', 'Cartn_z']);
}

export const mmCIF_Export_Filters = {
    onlyPositions: <CifWriter.Category.Filter>{
        includeCategory(name) { return name === 'atom_site'; },
        includeField(cat, field) { return _Filters.AtomSitePositionsFieldNames.has(field); }
    }
};

function getCustomPropCategories(customProp: CustomPropertyDescriptor, ctx: CifExportContext, params?: encode_mmCIF_categories_Params): CifExportCategoryInfo[] {
    if (!customProp.cifExport || customProp.cifExport.categories.length === 0) return [];

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

    const ret: CifExportCategoryInfo[] = [];
    for (const cat of cats) {
        if (params?.skipCategoryNames?.has(cat.name)) continue;
        if (params?.includedCategoryNames && !params.includedCategoryNames.has(cat.name)) continue;
        if (cat.name.indexOf(prefix) !== 0) throw new Error(`Custom category '${cat.name}' name must start with prefix '${prefix}.'`);
        ret.push([cat, propCtx]);
    }
    return ret;
}

type encode_mmCIF_categories_Params = {
    exportCtx?: CifExportContext,
    encoder?: CifWriter.Encoder,

    /** Skip provided categories */
    skipCategoryNames?: Set<string>,
    /** If defined, include only specified categories */
    includedCategoryNames?: Set<string>,
    /** If true, copy all categories present in an input mmCIF file */
    copyAllCategories?: boolean,
    /** If enabled, keep atom_site.id from the input data, otherwise use 1-based index of atom */
    doNotReindexAtomSiteId?: boolean,
    /** List of custom properties to include */
    customProperties?: CustomPropertyDescriptor[],

    extensions?: {
        /**
         * If specified, includes molstar_bond_site category with explicit bonds based on atom_site.id.
         * Forces doNotReindexAtomSiteId to true.
         * Note: This is not a standard mmCIF category and is currently specific to Mol*.
         */
        molstar_bond_site?: boolean,
    }
}

/** Doesn't start a data block */
export function encode_mmCIF_categories(encoder: CifWriter.Encoder, structures: Structure | Structure[], params?: encode_mmCIF_categories_Params) {
    const first = Array.isArray(structures) ? structures[0] : (structures as Structure);
    const models = first.models;
    if (models.length !== 1) throw new Error('Can\'t export stucture composed from multiple models.');

    const ctx: CifExportContext = params?.exportCtx || CifExportContext.create(structures);

    if (params?.copyAllCategories && MmcifFormat.is(models[0].sourceData)) {
        encode_mmCIF_categories_copyAll(encoder, ctx, params);
    } else {
        encode_mmCIF_categories_default(encoder, ctx, params);
    }
}

function encode_mmCIF_categories_default(encoder: CifWriter.Encoder, ctx: CifExportContext, params?: encode_mmCIF_categories_Params) {
    for (const cat of Categories({
        keepAtomSiteId: params?.extensions?.molstar_bond_site || params?.doNotReindexAtomSiteId
    })) {
        if (params?.skipCategoryNames && params?.skipCategoryNames.has(cat.name)) continue;
        if (params?.includedCategoryNames && !params.includedCategoryNames.has(cat.name)) continue;
        encoder.writeCategory(cat, ctx);
    }

    if (!params?.skipCategoryNames?.has('atom_site') && encoder.isCategoryIncluded('atom_site')) {
        const info = atom_site_operator_mapping(ctx);
        if (info) encoder.writeCategory(info[0], info[1], info[2]);
    }

    if (params?.extensions?.molstar_bond_site) {
        const info = molstar_bond_site(ctx);
        if (info) encoder.writeCategory(info[0], info[1], info[2]);
    }

    const _params = params || { };
    for (const customProp of ctx.firstModel.customProperties.all) {
        for (const [cat, propCtx] of getCustomPropCategories(customProp, ctx, _params)) {
            encoder.writeCategory(cat, propCtx);
        }
    }

    if (params?.customProperties) {
        for (const customProp of params?.customProperties) {
            for (const [cat, propCtx] of getCustomPropCategories(customProp, ctx, _params)) {
                encoder.writeCategory(cat, propCtx);
            }
        }
    }

    for (const s of ctx.structures) {
        if (!s.hasCustomProperties) continue;
        for (const customProp of s.customPropertyDescriptors.all) {
            for (const [cat, propCtx] of getCustomPropCategories(customProp, ctx, _params)) {
                encoder.writeCategory(cat, propCtx);
            }
        }
    }
}

function encode_mmCIF_categories_copyAll(encoder: CifWriter.Encoder, ctx: CifExportContext, params?: encode_mmCIF_categories_Params) {
    const providedCategories = new Map<string, CifExportCategoryInfo>();

    for (const cat of Categories({
        keepAtomSiteId: params?.extensions?.molstar_bond_site || params?.doNotReindexAtomSiteId
    })) {
        providedCategories.set(cat.name, [cat, ctx]);
    }

    const mapping = atom_site_operator_mapping(ctx);
    if (mapping) providedCategories.set(mapping[0].name, mapping);

    const _params = params || { };
    for (const customProp of ctx.firstModel.customProperties.all) {
        for (const info of getCustomPropCategories(customProp, ctx, _params)) {
            providedCategories.set(info[0].name, info);
        }
    }

    if (params?.customProperties) {
        for (const customProp of params?.customProperties) {
            for (const info of getCustomPropCategories(customProp, ctx, _params)) {
                providedCategories.set(info[0].name, info);
            }
        }
    }

    for (const s of ctx.structures) {
        if (!s.hasCustomProperties) continue;
        for (const customProp of s.customPropertyDescriptors.all) {
            for (const info of getCustomPropCategories(customProp, ctx)) {
                providedCategories.set(info[0].name, info);
            }
        }
    }

    const handled = new Set<string>();

    const data = (ctx.firstModel.sourceData as MmcifFormat).data;
    for (const catName of data.frame.categoryNames) {
        handled.add(catName);

        if (providedCategories.has(catName)) {
            const info = providedCategories.get(catName)!;
            encoder.writeCategory(info[0], info[1], info[2]);
        } else {
            if ((data.db as any)[catName]) {
                const cat = copy_mmCif_category(catName as any);
                encoder.writeCategory(cat, ctx);
            } else {
                const cat = copy_source_mmCifCategory(encoder, ctx, data.frame.categories[catName]);
                if (cat) encoder.writeCategory(cat);
            }
        }
    }

    providedCategories.forEach((info, name) => {
        if (!handled.has(name)) encoder.writeCategory(info[0], info[1], info[2]);
    });

    if (params?.extensions?.molstar_bond_site) {
        const info = molstar_bond_site(ctx);
        if (info) encoder.writeCategory(info[0], info[1], info[2]);
    }
}


function to_mmCIF(name: string, structure: Structure, asBinary = false, params?: encode_mmCIF_categories_Params) {
    const enc = params?.encoder ?? CifWriter.createEncoder({ binary: asBinary });
    enc.startDataBlock(name);
    encode_mmCIF_categories(enc, structure, params);
    return enc.getData();
}

export { to_mmCIF };