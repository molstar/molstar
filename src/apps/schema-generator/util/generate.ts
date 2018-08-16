/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Filter, Column } from './schema'
import { indentString } from 'mol-util';

function header (name: string, info: string, importDatabasePath = 'mol-data/db') {
    return `/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated '${name}' schema file. ${info}
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
export interface ${name}_Database extends Database<${name}_Schema> {}`
}

function getTypeDef(c: Column): string {
    switch (c.type) {
        case 'str': return 'str'
        case 'int': return 'int'
        case 'float': return 'float'
        case 'coord': return 'coord'
        case 'enum':
            return `Aliased<'${c.values.map(v => v.replace(/'/g, '\\\'')).join(`' | '`)}'>(${c.subType})`
        case 'matrix':
            return `Matrix(${c.rows}, ${c.columns})`
        case 'vector':
            return `Vector(${c.length})`
        case 'list':
            if (c.subType === 'int') {
                return `List('${c.separator}', x => parseInt(x, 10))`
            } else if (c.subType === 'float' || c.subType === 'coord') {
                return `List('${c.separator}', x => parseFloat(x))`
            } else {
                return `List('${c.separator}', x => x)`
            }
    }
}

const reSafePropertyName = /^[a-zA-Z_$][0-9a-zA-Z_$]*$/
function safePropertyString(name: string) { return name.match(reSafePropertyName) ? name : `'${name}'` }

export function generate (name: string, info: string, schema: Database, fields?: Filter, importDatabasePath?: string) {
    const codeLines: string[] = []

    codeLines.push(`export const ${name}_Schema = {`)
    Object.keys(schema).forEach(table => {
        if (fields && !fields[table]) return
        codeLines.push(`    ${safePropertyString(table)}: {`)
        const columns = schema[table]
        Object.keys(columns).forEach(columnName => {
            if (fields && !fields[table][columnName]) return
            const typeDef = getTypeDef(columns[columnName])
            if (columns[columnName].description) {
                codeLines.push(`        /**`)
                codeLines.push(`${indentString(columns[columnName].description, 1, '         * ')}`)
                codeLines.push(`         */`)
            }
            codeLines.push(`        ${safePropertyString(columnName)}: ${typeDef},`)
        })
        codeLines.push('    },')
    })
    codeLines.push('}')

    return `${header(name, info, importDatabasePath)}\n\n${codeLines.join('\n')}\n${footer(name)}`
}
