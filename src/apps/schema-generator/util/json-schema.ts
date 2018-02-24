/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export interface Database {
    [ tableName: string ]: Table
}

export interface Table {
    [ columnName: string ]: Column
}

export type Column = IntCol | StrCol | FloatCol | CoordCol | EnumCol | VectorCol | MatrixCol

type IntCol = 'int'
type StrCol = 'str'
type FloatCol = 'float'
type CoordCol = 'coord'

interface ComplexColumn {
    [ fieldType: string ]: any[]
}

interface EnumCol extends ComplexColumn {
    enum: string[]
}

interface VectorCol extends ComplexColumn {
    vector: [ number ]
}

interface MatrixCol extends ComplexColumn {
    matrix: [ number, number ]
}

export function getTypeAndArgs (column: ComplexColumn) {
    const type = Object.keys(column)[0] as string
    const args = column[ type ]
    return { type, args }
}

export type Filter = { [ table: string ]: { [ column: string ]: true } }

export function mergeFilters (...filters: Filter[]) {
    const n = filters.length
    const mergedFilter: Filter = {}
    const fields: Map<string, number> = new Map()
    filters.forEach(filter => {
        Object.keys(filter).forEach(category => {
            Object.keys(filter[ category ]).forEach(field => {
                const key = `${category}.${field}`
                const value = fields.get(key) || 0
                fields.set(key, value + 1)
            })
        })
    })
    fields.forEach((v, k) => {
        if (v !== n) return
        const [categoryName, fieldName] = k.split('.')
        if (categoryName in mergedFilter) {
            mergedFilter[categoryName][fieldName] = true
        } else {
            mergedFilter[categoryName] = { fieldName: true }
        }
    })
    return mergedFilter
}
