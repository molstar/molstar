
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

export function getSchema (dic: Data.Block) {  // todo Block needs to be specialized with safe frames as well
    const schema: { [category: string]: Category.Schema } = {}

    dic.saveFrames.forEach(d => {
        if (d.header[0] !== '_') {
            schema[d.header] = {}
        } else {
            const categoryName = d.header.substring(1, d.header.indexOf('.'))
            const itemName = d.header.substring(d.header.indexOf('.') + 1)
            let fields
            if (categoryName in schema) {
                fields = schema[categoryName]
            } else {
                fields = {}
                schema[categoryName] = fields
            }
            // console.log(util.inspect(d.categories, {showHidden: false, depth: 1}))
            const item_type = d.categories['_item_type']
            if (item_type) {
                const code = item_type.getField('code')
                if (code) {
                    fields[itemName] = getFieldType(code.str(0))
                } else {
                    console.log(`item_type.code not found for '${d.header}'`)
                }
            } else {
                // TODO check for _item_linked.parent_name and use its type
                console.log(`item_type not found for '${d.header}'`)
            }

        }
    })

    return schema
}
