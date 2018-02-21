/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { validate } from './validate'
import { Database, getTypeAndArgs, Filter } from './json-schema'

function header (name: string, importDatabasePath = 'mol-base/collections/database') {
    return `/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated '${name}' schema file
 *
 * @author mol-star package (src/scripts/schema-generation/generate)
 */

import { Database, Column } from '${importDatabasePath}'

import Schema = Column.Schema

const str = Schema.str;
const int = Schema.int;
const float = Schema.float;
const coord = Schema.coord;

const Aliased = Schema.Aliased;
const Matrix = Schema.Matrix;
const Vector = Schema.Vector;`
}

function footer (name: string) {
    return `
export type ${name}_Schema = typeof ${name}_Schema;
export interface ${name}_Database extends Database<${name}_Schema> { }`
}

const value: { [k: string]: (...args: any[]) => string } = {
    enum: function (...values: string[]) {
        return `Aliased<'${values.join(`' | '`)}'>(str)`
    },
    matrix: function (rows: number, cols: number) {
        return `Matrix(${rows}, ${cols})`
    },
    vector: function (dim: number) {
        return `Vector(${dim})`
    }
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
        codeLines.push(`\t'${table}': {`)
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
            codeLines.push(`\t\t'${columnName}': ${typeDef},`)
        })
        codeLines.push('\t},')
    })
    codeLines.push('}')

    return `${header(name, importDatabasePath)}\n\n${codeLines.join('\n')}\n${footer(name)}`
}
