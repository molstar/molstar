/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { validate } from './validate'
import { Database, getTypeAndArgs, Filter } from './json-schema'

function header (name: string, importDatabasePath = 'mol-data/db') {
    return `/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated '${name}' schema file
 *
 * @author mol-star package (src/apps/schema-generator/generate)
 */

import { Database, Column } from '${importDatabasePath}'

import Schema = Column.Schema

const str = Schema.str;
const int = Schema.int;
const float = Schema.float;
const coord = Schema.coord;

const Aliased = Schema.Aliased;
const Matrix = Schema.Matrix;
const Vector = Schema.Vector;
const List = Schema.List;`
}

function footer (name: string) {
    return `
export type ${name}_Schema = typeof ${name}_Schema;
export type ${name}_Database = Database<${name}_Schema>`
}

const value: { [k: string]: (...args: any[]) => string } = {
    enum: function (type: string, values: string[]) {
        return `Aliased<'${values.join(`' | '`)}'>(${type})`
    },
    matrix: function (rows: number, cols: number) {
        return `Matrix(${rows}, ${cols})`
    },
    vector: function (dim: number) {
        return `Vector(${dim})`
    },
    list: function (type: 'str'|'int'|'float', separator: string) {
        if (type === 'int') {
            return `List('${separator}', x => parseInt(x, 10))`
        } else if (type === 'float') {
            return `List('${separator}', x => parseFloat(x))`
        } else {
            return `List('${separator}', x => x)`
        }
    }
}

const reSafePropertyName = /^[a-zA-Z_$][0-9a-zA-Z_$]*$/
function safePropertyString(name: string) {
    return name.match(reSafePropertyName) ? name : `'${name}'`
}

export function generate (name: string, schema: Database, fields?: Filter, importDatabasePath?: string) {
    const validationResult = validate(schema)
    if (validationResult !== true) {
        throw validationResult
    }

    const codeLines: string[] = []

    codeLines.push(`export const ${name}_Schema = {`)
    Object.keys(schema).forEach(table => {
        if (fields && !fields[ table ]) return
        codeLines.push(`    ${safePropertyString(table)}: {`)
        const columns = schema[ table ]
        Object.keys(columns).forEach(columnName => {
            if (fields && !fields[ table ][ columnName ]) return
            let typeDef
            const fieldType = columns[ columnName ]
            if (typeof fieldType === 'object') {
                const { type, args } = getTypeAndArgs(fieldType)
                typeDef = value[ type ](...args)
            } else {
                typeDef = fieldType
            }
            codeLines.push(`        ${safePropertyString(columnName)}: ${typeDef},`)
        })
        codeLines.push('    },')
    })
    codeLines.push('}')

    return `${header(name, importDatabasePath)}\n\n${codeLines.join('\n')}\n${footer(name)}`
}
