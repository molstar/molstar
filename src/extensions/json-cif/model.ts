/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Table } from '../../mol-data/db';

export const JSONCifVERSION = '0.1.0';

export interface JSONCifFile {
    version: string;
    encoder: string;
    dataBlocks: JSONCifDataBlock[];
}

export interface JSONCifDataBlock {
    header: string,
    categoryNames: string[],
    categories: Record<string, JSONCifCategory>,
}

export interface JSONCifCategory<T extends Record<string, any> = Record<string, any>> {
    name: string,
    fieldNames: string[],
    rows: T[],
}

export function getJSONCifCategory<S extends Table.Schema>(block: JSONCifDataBlock, name: string): JSONCifCategory<Table.Row<S>> | undefined {
    return block.categories[name] as any;
}