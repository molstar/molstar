/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Column, EnumCol, StrCol, IntCol, ListCol, FloatCol, CoordCol, MatrixCol, VectorCol } from './schema'
import * as Data from 'mol-io/reader/cif/data-model'
import { CifFrame } from 'mol-io/reader/cif/data-model';

export function getFieldType (type: string, description: string, values?: string[]): Column {
    switch (type) {
        case 'code':
        case 'ucode':
        case 'line':
        case 'uline':
        case 'text':
        case 'char':
        case 'uchar3':
        case 'uchar1':
        case 'boolean':
            return values && values.length ? EnumCol(values, 'str', description) : StrCol(description)
        case 'aliasname':
        case 'name':
        case 'idname':
        case 'any':
        case 'atcode':
        case 'fax':
        case 'phone':
        case 'email':
        case 'code30':
        case 'seq-one-letter-code':
        case 'author':
        case 'orcid_id':
        case 'sequence_dep':
        case 'pdb_id':
        case 'emd_id':
        // todo, consider adding specialised fields
        case 'yyyy-mm-dd':
        case 'yyyy-mm-dd:hh:mm':
        case 'yyyy-mm-dd:hh:mm-flex':
        case 'int-range':
        case 'float-range':
        case 'binary':
        case 'operation_expression':
        case 'point_symmetry':
        case '4x3_matrix':
        case '3x4_matrices':
        case 'point_group':
        case 'point_group_helical':
        case 'symmetry_operation':
        case 'date_dep':
        case 'url':
        case 'symop':
        case 'exp_data_doi':
        case 'asym_id':
            return StrCol(description)
        case 'int':
        case 'non_negative_int':
        case 'positive_int':
           return values && values.length ? EnumCol(values, 'int', description) : IntCol(description)
        case 'float':
            return FloatCol(description)
        case 'ec-type':
        case 'ucode-alphanum-csv':
        case 'id_list':
            return ListCol('str', ',', description)
        case 'id_list_spc':
            return ListCol('str', ' ', description)
    }
    console.log(`unknown type '${type}'`)
    return StrCol(description)
}

type FrameCategories = { [category: string]: Data.CifFrame }
type FrameLinks = { [k: string]: string }

interface FrameData {
    categories: FrameCategories
    links: FrameLinks
}

// get field from given or linked category
function getField (category: string, field: string, d: Data.CifFrame, ctx: FrameData): Data.CifField|undefined {
    const { categories, links } = ctx

    const cat = d.categories[category]
    if (cat) {
        return cat.getField(field)
    } else {
        if (d.header in links) {
            const linkName = links[d.header]
            if (linkName in categories) {
                return getField(category, field, categories[linkName], ctx)
            } else {
                console.log(`link '${linkName}' not found`)
            }
        } else {
            // console.log(`no links found for '${d.header}'`)
        }
    }
}

function getEnums (d: Data.CifFrame, ctx: FrameData) {
    const value = getField('item_enumeration', 'value', d, ctx)
    const enums: string[] = []
    if (value) {
        for (let i = 0; i < value.rowCount; ++i) {
            enums.push(value.str(i))
            // console.log(value.str(i))
        }
        return enums
    } else {
        // console.log(`item_enumeration.value not found for '${d.header}'`)
    }
}

function getCode (d: Data.CifFrame, ctx: FrameData): [string, string[]|undefined]|undefined {
    const code = getField('item_type', 'code', d, ctx)
    if (code) {
        return [ code.str(0), getEnums(d, ctx) ]
    } else {
        console.log(`item_type.code not found for '${d.header}'`)
    }
}

function getSubCategory (d: Data.CifFrame, ctx: FrameData): string|undefined {
    const value = getField('item_sub_category', 'id', d, ctx)
    if (value) {
        return value.str(0)
    }
}

function getDescription (d: Data.CifFrame, ctx: FrameData): string|undefined {
    const value = getField('item_description', 'description', d, ctx)
    if (value) {
        // trim (after newlines) and remove references to square brackets
        return value.str(0).trim()
            .replace(/(\r\n|\r|\n)([ \t]+)/g, '\n')
            .replace(/(\[[1-3]\])+ element/, 'elements')
            .replace(/(\[[1-3]\])+/, '')
    }
}

const reMatrixField = /\[[1-3]\]\[[1-3]\]/
const reVectorField = /\[[1-3]\]/

const FORCE_INT_FIELDS = [
    '_atom_site.id',
    '_atom_site.auth_seq_id',
    '_pdbx_struct_mod_residue.auth_seq_id',
    '_struct_conf.beg_auth_seq_id',
    '_struct_conf.end_auth_seq_id',
    '_struct_conn.ptnr1_auth_seq_id',
    '_struct_conn.ptnr2_auth_seq_id',
    '_struct_sheet_range.beg_auth_seq_id',
    '_struct_sheet_range.end_auth_seq_id',
];

const COMMA_SEPARATED_LIST_FIELDS = [
    '_atom_site.pdbx_struct_group_id',
    '_chem_comp.mon_nstd_parent_comp_id',
    '_diffrn_radiation.pdbx_wavelength_list',
    '_diffrn_source.pdbx_wavelength_list',
    '_em_diffraction.tilt_angle_list', // 20,40,50,55
    '_em_entity_assembly.entity_id_list',
    '_entity.pdbx_description', // Endolysin,Beta-2 adrenergic receptor
    '_entity.pdbx_ec',
    '_entity_poly.pdbx_strand_id', // A,B
    '_entity_src_gen.pdbx_gene_src_gene', // ADRB2, ADRB2R, B2AR
    '_pdbx_depui_entry_details.experimental_methods',
    '_pdbx_depui_entry_details.requested_accession_types',
    '_pdbx_soln_scatter_model.software_list', // INSIGHT II, HOMOLOGY, DISCOVERY, BIOPOLYMER, DELPHI
    '_pdbx_soln_scatter_model.software_author_list', // MSI
    '_pdbx_soln_scatter_model.entry_fitting_list', // Odd example: 'PDB CODE 1HFI, 1HCC, 1HFH, 1VCC'
    '_pdbx_struct_assembly_gen.entity_inst_id',
    '_pdbx_struct_assembly_gen.asym_id_list',
    '_pdbx_struct_assembly_gen.auth_asym_id_list',
    '_pdbx_struct_assembly_gen_depositor_info.asym_id_list',
    '_pdbx_struct_assembly_gen_depositor_info.chain_id_list',
    '_pdbx_struct_group_list.group_enumeration_type',
    '_reflns.pdbx_diffrn_id',
    '_refine.pdbx_diffrn_id',
    '_reflns_shell.pdbx_diffrn_id',
    '_struct_keywords.text',
];

const SPACE_SEPARATED_LIST_FIELDS = [
    '_chem_comp.pdbx_subcomponent_list', // TSM DPH HIS CHF EMR
    '_pdbx_soln_scatter.data_reduction_software_list', // OTOKO
    '_pdbx_soln_scatter.data_analysis_software_list', // SCTPL5 GNOM
];

const SEMICOLON_SEPARATED_LIST_FIELDS = [
    '_chem_comp.pdbx_synonyms' // GLYCERIN; PROPANE-1,2,3-TRIOL
]

/**
 * Useful when a dictionary extension will add enum values to an existing dictionary.
 * By adding them here, the dictionary extension can be tested before the added enum
 * values are available in the existing dictionary.
 */
const EXTRA_ENUM_VALUES: { [k: string]: string[] } = {

}

export function generateSchema (frames: CifFrame[]) {
    const schema: Database = {}

    const categories: FrameCategories = {}
    const links: FrameLinks = {}
    const ctx = { categories, links }

    // get category metadata
    frames.forEach(d => {
        if (d.header[0] === '_') return
        const categoryKeyNames = new Set<string>()
        const categoryKey = d.categories['category_key']
        if (categoryKey) {
            const categoryKey_names = categoryKey.getField('name')
            if (categoryKey_names) {
                for (let i = 0, il = categoryKey_names.rowCount; i < il; ++i) {
                    categoryKeyNames.add(categoryKey_names.str(i))
                }
            }
        }
        let description = ''
        const category = d.categories['category']
        if (category) {
            const category_description = category.getField('description')
            if (category_description) {
                description = category_description.str(0).trim()
                    .replace(/(\r\n|\r|\n)([ \t]+)/g, '\n') // remove padding after newlines
            } else {
                console.log(`no description given for category '${category}'`)
            }
        }
        if (categoryKeyNames.size === 0) {
            console.log(`no key given for category '${category}'`)
        }
        schema[d.header] = { description, key: categoryKeyNames, columns: {} }
        // console.log('++++++++++++++++++++++++++++++++++++++++++')
        // console.log('name', d.header)
        // console.log('desc', description)
        // console.log('key', categoryKeyNames)
    })

    // build list of links between categories
    frames.forEach(d => {
        if (d.header[0] !== '_') return
        categories[d.header] = d
        const item_linked = d.categories['item_linked']
        if (item_linked) {
            const child_name = item_linked.getField('child_name')
            const parent_name = item_linked.getField('parent_name')
            if (child_name && parent_name) {
                for (let i = 0; i < item_linked.rowCount; ++i) {
                    const childName = child_name.str(i)
                    const parentName = parent_name.str(i)
                    if (childName in links && links[childName] !== parentName) {
                        console.log(`${childName} linked to ${links[childName]}, ignoring link to ${parentName}`)
                    }
                    links[childName] = parentName
                }
            }
        }
    })

    // get field data
    Object.keys(categories).forEach(fullName => {
        const d = categories[fullName]
        if (!d) {
            console.log(`${fullName} not found, moving on`)
            return
        }
        const categoryName = d.header.substring(1, d.header.indexOf('.'))
        const itemName = d.header.substring(d.header.indexOf('.') + 1)
        let fields: { [k: string]: Column }
        if (categoryName in schema) {
            fields = schema[categoryName].columns
        } else {
            console.log(`category '${categoryName}' has no metadata`)
            fields = {}
            schema[categoryName] = {
                description: '',
                key: new Set(),
                columns: fields
            }
        }

        const description = getDescription(d, ctx) || ''

        // need to use regex to check for matrix or vector items
        // as sub_category assignment is missing for some entries
        const subCategory = getSubCategory(d, ctx)
        if (subCategory === 'cartesian_coordinate' || subCategory === 'fractional_coordinate') {
            fields[itemName] = CoordCol(description)
        } else if (FORCE_INT_FIELDS.includes(d.header)) {
            fields[itemName] = IntCol(description)
            console.log(`forcing int: ${d.header}`)
        } else if (subCategory === 'matrix') {
            fields[itemName.replace(reMatrixField, '')] = MatrixCol(3, 3, description)
        } else if (subCategory === 'vector') {
            fields[itemName.replace(reVectorField, '')] = VectorCol(3, description)
        } else {
            if (itemName.match(reMatrixField)) {
                fields[itemName.replace(reMatrixField, '')] = MatrixCol(3, 3, description)
                console.log(`${d.header} should have 'matrix' _item_sub_category.id`)
            } else if (itemName.match(reVectorField)) {
                fields[itemName.replace(reVectorField, '')] = VectorCol(3, description)
                console.log(`${d.header} should have 'vector' _item_sub_category.id`)
            } else {
                const code = getCode(d, ctx)
                if (code) {
                    let fieldType = getFieldType(code[0], description, code[1]);
                    if (fieldType.type === 'str') {
                        if (COMMA_SEPARATED_LIST_FIELDS.includes(d.header)) {
                            fieldType = ListCol('str', ',', description)
                            console.log(`forcing comma separated: ${d.header}`)
                        } else if (SPACE_SEPARATED_LIST_FIELDS.includes(d.header)) {
                            fieldType = ListCol('str', ' ', description)
                            console.log(`forcing space separated: ${d.header}`)
                        } else if (SEMICOLON_SEPARATED_LIST_FIELDS.includes(d.header)) {
                            fieldType = ListCol('str', ';', description)
                            console.log(`forcing space separated: ${d.header}`)
                        }
                    }
                    if (d.header in EXTRA_ENUM_VALUES) {
                        if (fieldType.type === 'enum') {
                            fieldType.values.push(...EXTRA_ENUM_VALUES[d.header])
                        } else {
                            console.warn(`expected enum: ${d.header}`)
                        }
                    }
                    fields[itemName] = fieldType
                } else {
                    console.log(`could not determine code for '${d.header}'`)
                }
            }
        }
    })

    return schema
}
