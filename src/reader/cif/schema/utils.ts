
// import dic from './dic'
import { Field, Block, Category } from '../schema'
import * as Data from '../data-model'

const pooledStr = Field.pooledStr()
const str = Field.str()
const int = Field.int()
const float = Field.float()

export function getFieldType (type: string) {
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

// get field from given or linked category
function getField ( category: string, field: string, d: Data.SafeFrame, ctx: SafeFrameData): Data.Field|undefined {
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

// function getEnums (d: Data.SafeFrame, ctx: SafeFrameData): string[]|undefined {
//     const value = getField('_item_enumeration', 'value', d, ctx)
//     if (value) {
//         const enums: string[] = []
//         for (let i = 0; i < value.rowCount; ++i) {
//             enums.push(value.str(i))
//             // console.log(value.str(i))
//         }
//         return enums
//     } else {
//         // console.log(`item_enumeration.value not found for '${d.header}'`)
//     }
// }

function getCode (d: Data.SafeFrame, ctx: SafeFrameData): string|undefined {
    const code = getField('_item_type', 'code', d, ctx)
    if (code) {
        let c = code.str(0)
        // if (c === 'ucode') {
        //     const enums = getEnums(d, ctx)
        //     if (enums) c += `: ${enums.join('|')}`
        // }
        return c
    } else {
        console.log(`item_type.code not found for '${d.header}'`)
    }
}

export function getSchema (dic: Data.Block) {  // todo Block needs to be specialized with safe frames as well
    const schema: Block.Schema = {}  // { [category: string]: Category.Schema } = {}

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

        const code = getCode(d, { categories, links })
        if (code) {
            fields[itemName] = getFieldType(code)
        } else {
            console.log(`could not determine code for '${d.header}'`)
        }
    })

    return schema as Block.Instance<any>
}
