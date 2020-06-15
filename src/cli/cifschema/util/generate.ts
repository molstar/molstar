/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Database, Filter, Column } from './schema';
import { indentString } from '../../../mol-util/string';
import { FieldPath } from '../../../mol-io/reader/cif/schema';

function header (name: string, info: string, moldataImportPath: string) {
    return `/**
 * Copyright (c) 2017-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Code-generated '${name}' schema file. ${info}
 *
 * @author molstar/ciftools package
 */

import { Database, Column } from '${moldataImportPath}/db';

import Schema = Column.Schema;`;
}

function footer (name: string) {
    return `
export type ${name}_Schema = typeof ${name}_Schema;
export interface ${name}_Database extends Database<${name}_Schema> {};`;
}

function getTypeShorthands(schema: Database, fields?: Filter) {
    const types = new Set<string>();
    Object.keys(schema.tables).forEach(table => {
        if (fields && !fields[table]) return;
        const { columns } = schema.tables[table];
        Object.keys(columns).forEach(columnName => {
            if (fields && !fields[table][columnName]) return;
            types.add(schema.tables[table].columns[columnName].type);
        });
    });
    const shorthands: string[] = [];
    types.forEach(type => {
        switch (type) {
            case 'str': shorthands.push('const str = Schema.str;'); break;
            case 'int': shorthands.push('const int = Schema.int;'); break;
            case 'float': shorthands.push('const float = Schema.float;'); break;
            case 'coord': shorthands.push('const coord = Schema.coord;'); break;
            case 'enum': shorthands.push('const Aliased = Schema.Aliased;'); break;
            case 'matrix': shorthands.push('const Matrix = Schema.Matrix;'); break;
            case 'vector': shorthands.push('const Vector = Schema.Vector;'); break;
            case 'list': shorthands.push('const List = Schema.List;'); break;
        }
    });
    return shorthands.join('\n');
}

function getTypeDef(c: Column): string {
    switch (c.type) {
        case 'str': return 'str';
        case 'int': return 'int';
        case 'float': return 'float';
        case 'coord': return 'coord';
        case 'enum':
            return `Aliased<'${c.values.map(v => v.replace(/'/g, '\\\'')).join(`' | '`)}'>(${c.subType})`;
        case 'matrix':
            return `Matrix(${c.rows}, ${c.columns})`;
        case 'vector':
            return `Vector(${c.length})`;
        case 'list':
            if (c.subType === 'int') {
                return `List('${c.separator}', x => parseInt(x, 10))`;
            } else if (c.subType === 'float' || c.subType === 'coord') {
                return `List('${c.separator}', x => parseFloat(x))`;
            } else {
                return `List('${c.separator}', x => x)`;
            }
    }
}

const reSafePropertyName = /^[a-zA-Z_$][0-9a-zA-Z_$]*$/;
function safePropertyString(name: string) { return name.match(reSafePropertyName) ? name : `'${name}'`; }

function doc(description: string, spacesCount: number) {
    const spaces = ' '.repeat(spacesCount);
    return [
        `${spaces}/**`,
        `${indentString(description, 1, `${spaces} * `)}`.replace(/ +\n/g, '\n'),
        `${spaces} */`
    ].join('\n');
}

export function generate (name: string, info: string, schema: Database, fields: Filter | undefined, moldataImportPath: string, addAliases: boolean) {
    const codeLines: string[] = [];

    if (fields) {
        Object.keys(fields).forEach(table => {
            if (table in schema.tables) {
                const schemaTable = schema.tables[table];
                Object.keys(fields[table]).forEach(column => {
                    if (!(column in schemaTable.columns)) {
                        console.log(`filter field '${table}.${column}' not found in schema`);
                    }
                });
            } else {
                console.log(`filter category '${table}' not found in schema`);
            }
        });
    }

    codeLines.push(`export const ${name}_Schema = {`);
    Object.keys(schema.tables).forEach(table => {
        if (fields && !fields[table]) return;
        const { description, columns } = schema.tables[table];
        if (description) codeLines.push(doc(description, 4));
        codeLines.push(`    ${safePropertyString(table)}: {`);
        Object.keys(columns).forEach(columnName => {
            if (fields && !fields[table][columnName]) return;
            const c = columns[columnName];
            const typeDef = getTypeDef(c);
            if (c.description) codeLines.push(doc(c.description, 8));
            codeLines.push(`        ${safePropertyString(columnName)}: ${typeDef},`);
        });
        codeLines.push('    },');
    });
    codeLines.push('};');

    if (addAliases) {
        codeLines.push('');
        codeLines.push(`export const ${name}_Aliases = {`);
        Object.keys(schema.aliases).forEach(path => {
            const [ table, columnName ] = path.split('.');
            if (fields && !fields[table]) return;
            if (fields && !fields[table][columnName]) return;

            const filteredAliases = new Set<string>();
            schema.aliases[path].forEach(p => {
                if (!FieldPath.equal(p, path)) filteredAliases.add(FieldPath.canonical(p));
            });

            if (filteredAliases.size === 0) return;
            codeLines.push(`    ${safePropertyString(path)}: [`);
            filteredAliases.forEach(alias => {
                codeLines.push(`        '${alias}',`);
            });
            codeLines.push('    ],');
        });
        codeLines.push('};');
    }

    return `${header(name, info, moldataImportPath)}\n\n${getTypeShorthands(schema, fields)}\n\n${codeLines.join('\n')}\n${footer(name)}`;
}