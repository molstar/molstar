
// import dic from './dic'
import { Field, Category } from '../schema'
import * as Data from '../data-model'

const pooledStr = Field.pooledStr()
const str = Field.str()
const int = Field.int()
const float = Field.float()

function getFieldType (type: string) {
    switch (type) {
        case 'code':
        case 'ucode':
        case 'line':
        case 'uline':
        case 'text':
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
        case 'point_group':
        case 'point_group_helical':
        case 'boolean':
        case 'symmetry_operation':
        case 'date_dep':
            return str
        case 'uchar3':
        case 'uchar1':
        case 'symop':
            return pooledStr
        case 'int':
        case 'non_negative_int':
        case 'positive_int':
            return int
        case 'float':
            return float
    }
    console.log(`unknown type '${type}'`)
    return str
}

type SafeFrameCategories = { [category: string]: Data.SafeFrame }
type SafeFrameLinks = { [k: string]: string }

interface SafeFrameData {
    categories: SafeFrameCategories
    links: SafeFrameLinks
}

function getCode (d: Data.SafeFrame, ctx: SafeFrameData): string|undefined {
    const { categories, links } = ctx

    const item_type = d.categories['_item_type']
    if (item_type) {
        const code = item_type.getField('code')
        if (code) {
            return code.str(0)
        } else {
            console.log(`item_type.code not found for '${d.header}'`)
        }
    } else {
        if (d.header in links) {
            return getCode(categories[links[d.header]], ctx)
        } else {
            console.log(`no links found for '${d.header}'`)
        }
    }
}

export function getSchema (dic: Data.Block) {  // todo Block needs to be specialized with safe frames as well
    const schema: { [category: string]: Category.Schema } = {}

    const categories: SafeFrameCategories = {}
    const links: SafeFrameLinks = {}
    dic.saveFrames.forEach(d => {
        if (d.header[0] !== '_') return
        categories[d.header] = d
        const item_linked = d.categories['_item_linked']
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

    Object.keys(categories).forEach(categoryName => {
        const d = categories[categoryName]
        const itemName = d.header.substring(d.header.indexOf('.') + 1)
        let fields
        if (categoryName in schema) {
            fields = schema[categoryName]
        } else {
            fields = {}
            schema[categoryName] = fields
        }

        const code = getCode(d, { categories, links })
        if (code) {
            fields[itemName] = getFieldType(code)
        } else {
            console.log(`could not determine code for '${d.header}'`)
        }
    })

    return schema
}

// TODO
// support controlled vocabulary as a specialization string type field
// in the example below the string type would be `Y|N`
// _item_type.code               ucode
    // loop_
    // _item_enumeration.value
    // _item_enumeration.detail
          // Y  'Yes'
          // N  'No'
