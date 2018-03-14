/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Table, Column } from './json-schema'

const SimpleColumnTypes = [ 'str', 'int', 'float', 'coord' ]
const ComplexColumnTypes = [ 'enum', 'vector', 'matrix', 'list' ]

function allTrue<T> (list: T[], fn: (e: T) => boolean) {
    return list.reduce((a, v) => a && fn(v), true)
}

function allString (list: string[]) {
    return list.reduce((a, v) => a && typeof v === 'string', true)
}

function validateColumn (column: Column): true|Error {
    if (typeof column === 'string') {
        if (!SimpleColumnTypes.includes(column)) {
            return new Error(`simple column types must be one of '${SimpleColumnTypes.join(', ')}' not '${column}'`)
        }
        return true
    } else if (typeof column === 'object') {
        const keys = Object.keys(column)
        if (keys.length !== 1) {
            return new Error(`complex column object must have a single key`)
        }
        const type = keys[0]
        const args = column[ type ]
        if (!Array.isArray(args)) {
            return new Error(`complex column args must be an array`)
        }
        switch (type) {
            case 'enum':
                if (args.length !== 2 && (!allString(args[1]) && !allTrue(args[1], Number.isInteger))) {
                    return new Error(`enum column must have all string or all integer args ${args}`)
                }
                break;
            case 'vector':
                if (args.length !== 1 || !allTrue(args, Number.isInteger)) {
                    return new Error(`vector column must have one integer arg`)
                }
                break;
            case 'matrix':
                if (args.length !== 2 || !allTrue(args, Number.isInteger)) {
                    return new Error(`matrix column must have two integer args`)
                }
                break;
            case 'list':
                if (args.length !== 2 || !allString(args)) {
                    return new Error(`list column must have two string args`)
                }
                break;
            default:
                return new Error(`complex column types must be one of '${ComplexColumnTypes.join(', ')}' not '${type}'`)
        }
        return true
    }
    return new Error(`columns must be of type 'object' or 'string' not '${typeof column}'`)
}

function validateTable (table: Table): true|Error {
    if (typeof table !== 'object') {
        return new Error(`table must be of type 'object' not '${typeof table}'`)
    }
    for (const columnName in table) {
        // could check columnName with regex
        const r = validateColumn(table[columnName])
        if (r !== true) {
            return new Error(`[${columnName}] ${r.message}`)
        }
    }
    return true
}

function validateDatabase (database: Database): true|Error {
    if (typeof database !== 'object') {
        return new Error(`database must be of type 'object' not '${typeof database}'`)
    }
    for (const tableName in database) {
        // could check tableName with regex
        const r = validateTable(database[tableName])
        if (r !== true) {
            return new Error(`[${tableName}] ${r.message}`)
        }
    }
    return true
}

export function validate (schema: any): true|Error {
    return validateDatabase(schema)
}
