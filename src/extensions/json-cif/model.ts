/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

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

export interface JSONCifCategory {
    name: string,
    fieldNames: string[],
    rows: Record<string, any>[],
}