/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Column } from './json-schema'
import * as Data from 'mol-io/reader/cif/data-model'

export function getFieldType (type: string, values?: string[]): Column {
    switch (type) {
        case 'code':
        case 'ucode':
            if (values && values.length) {
                return { 'enum': values }
            } else {
                return 'str'
            }
        case 'line':
        case 'uline':
        case 'text':
        case 'char':
        case 'aliasname':
        case 'name':
        case 'idname':
        case 'any':
        case 'atcode':
        case 'fax':
        case 'phone':
        case 'email':
        case 'code30':
        case 'ec-type':
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
        case 'ucode-alphanum-csv':
        case 'point_symmetry':
        case 'id_list':
        case '4x3_matrix':
        case '3x4_matrices':
        case 'point_group':
        case 'point_group_helical':
        case 'boolean':
        case 'symmetry_operation':
        case 'date_dep':
        case 'uchar3':
        case 'uchar1':
        case 'url':
        case 'symop':
            return 'str'
        case 'int':
        case 'non_negative_int':
        case 'positive_int':
            return 'int'
        case 'float':
            return 'float'
    }
    console.log(`unknown type '${type}'`)
    return 'str'
}

type FrameCategories = { [category: string]: Data.Frame }
type FrameLinks = { [k: string]: string }

interface FrameData {
    categories: FrameCategories
    links: FrameLinks
}

// get field from given or linked category
function getField ( category: string, field: string, d: Data.Frame, ctx: FrameData): Data.Field|undefined {
    const { categories, links } = ctx

    const cat = d.categories[category]
    if (cat) {
        return cat.getField(field)
    } else {
        if (d.header in links) {
            return getField(category, field, categories[links[d.header]], ctx)
        } else {
            // console.log(`no links found for '${d.header}'`)
        }
    }
}

function getEnums (d: Data.Frame, ctx: FrameData): string[]|undefined {
    const value = getField('item_enumeration', 'value', d, ctx)
    if (value) {
        const enums: string[] = []
        for (let i = 0; i < value.rowCount; ++i) {
            enums.push(value.str(i))
            // console.log(value.str(i))
        }
        return enums
    } else {
        // console.log(`item_enumeration.value not found for '${d.header}'`)
    }
}

function getCode (d: Data.Frame, ctx: FrameData): [string, string[]]|undefined {
    const code = getField('item_type', 'code', d, ctx)
    if (code) {
        let c = code.str(0)
        let e = []
        if (c === 'ucode') {
            const enums = getEnums(d, ctx)
            if (enums) e.push(...enums)
        }
        return [c, e]
    } else {
        console.log(`item_type.code not found for '${d.header}'`)
    }
}

function getSubCategory (d: Data.Frame, ctx: FrameData): string|undefined {
    const value = getField('item_sub_category', 'id', d, ctx)
    if (value) {
        return value.str(0)
    }
}

const FORCE_INT_FIELDS = [
    '_struct_conf.beg_auth_seq_id',
    '_struct_conf.end_auth_seq_id',
    '_struct_sheet_range.beg_auth_seq_id',
    '_struct_sheet_range.end_auth_seq_id',
    '_struct_conn.ptnr1_auth_seq_id',
    '_struct_conn.ptnr2_auth_seq_id',
    '_pdbx_struct_mod_residue.auth_seq_id',
    '_atom_site.id',
    '_atom_site.auth_seq_id'
];

export function generateSchema (dic: Data.Block) {
    const schema: Database = {}

    const categories: FrameCategories = {}
    const links: FrameLinks = {}
    const ctx = { categories, links }

    dic.saveFrames.forEach(d => {
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

    Object.keys(categories).forEach(fullName => {
        const d = categories[fullName]
        const categoryName = d.header.substring(1, d.header.indexOf('.'))
        const itemName = d.header.substring(d.header.indexOf('.') + 1)
        let fields
        if (categoryName in schema) {
            fields = schema[categoryName]
        } else {
            fields = {}
            schema[categoryName] = fields
        }

        // need to use regex to check for matrix or vector items
        // as sub_category assignment is missing for some entries
        const subCategory = getSubCategory(d, ctx)
        if (subCategory === 'cartesian_coordinate' || subCategory === 'fractional_coordinate') {
            fields[itemName] = 'coord'
        } else if (FORCE_INT_FIELDS.includes(d.header)) {
            fields[itemName] = 'int'
        } else if (subCategory === 'matrix') {
            fields[itemName.replace(/\[[1-3]\]\[[1-3]\]/, '')] = { 'matrix': [ 3, 3 ] }
        } else if (subCategory === 'vector') {
            fields[itemName.replace(/\[[1-3]\]/, '')] = { 'vector': [ 3 ] }
        } else {
            if (itemName.match(/\[[1-3]\]\[[1-3]\]/)) {
                fields[itemName.replace(/\[[1-3]\]\[[1-3]\]/, '')] = { 'matrix': [ 3, 3 ] }
                // console.log(`${d.header} should have 'matrix' _item_sub_category.id`)
            } else if (itemName.match(/\[[1-3]\]/)) {
                fields[itemName.replace(/\[[1-3]\]/, '')] = { 'vector': [ 3 ] }
                // console.log(`${d.header} should have 'vector' _item_sub_category.id`)
            } else {
                const code = getCode(d, ctx)
                if (code) {
                    fields[itemName] = getFieldType(code[0], code[1])
                } else {
                    console.log(`could not determine code for '${d.header}'`)
                }
            }
        }
    })

    return schema
}
