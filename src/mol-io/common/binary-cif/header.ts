/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createMsgPackState, MsgPackState, parseMsgPackValue, readMsgPackArrayLength, readMsgPackMapLength, skipMsgPackValue } from '../msgpack/decode';

/**
 * Metadata of a BinaryCIF file without any column data.
 * Enough to content-sniff a file (e.g. in `DataFormatProvider.isApplicable`) at a small
 * fraction of the cost of decoding it with `decodeMsgPack`.
 */
export interface BinaryCifHeader {
    version: string,
    encoder: string,
    dataBlocks: BinaryCifBlockHeader[]
}

export interface BinaryCifBlockHeader {
    header: string,
    categories: BinaryCifCategoryHeader[]
}

export interface BinaryCifCategoryHeader {
    /** Category name including the leading underscore, e.g. `_atom_site`. */
    name: string,
    rowCount: number,
    columnNames: string[]
}

const HeaderCache = new WeakMap<Uint8Array, BinaryCifHeader>();

/** Memoized `decodeBinaryCifHeader`, keyed by the data buffer instance. */
export function getBinaryCifHeader(data: Uint8Array): BinaryCifHeader {
    const cached = HeaderCache.get(data);
    if (cached) return cached;
    const header = decodeBinaryCifHeader(data);
    HeaderCache.set(data, header);
    return header;
}

/** Throws if the data is not valid MessagePack, mirroring `decodeMsgPack`. */
export function decodeBinaryCifHeader(data: Uint8Array): BinaryCifHeader {
    const state = createMsgPackState(data);
    const result: BinaryCifHeader = { version: '', encoder: '', dataBlocks: [] };

    const count = readMsgPackMapLength(state);
    for (let i = 0; i < count; i++) {
        switch (parseMsgPackValue(state)) {
            case 'version': result.version = parseMsgPackValue(state); break;
            case 'encoder': result.encoder = parseMsgPackValue(state); break;
            case 'dataBlocks': result.dataBlocks = readDataBlocks(state); break;
            default: skipMsgPackValue(state);
        }
    }
    return result;
}

function readDataBlocks(state: MsgPackState) {
    const blocks: BinaryCifBlockHeader[] = [];
    const blockCount = readMsgPackArrayLength(state);
    for (let i = 0; i < blockCount; i++) {
        const block: BinaryCifBlockHeader = { header: '', categories: [] };
        const count = readMsgPackMapLength(state);
        for (let j = 0; j < count; j++) {
            switch (parseMsgPackValue(state)) {
                case 'header': block.header = parseMsgPackValue(state); break;
                case 'categories': block.categories = readCategories(state); break;
                default: skipMsgPackValue(state);
            }
        }
        blocks.push(block);
    }
    return blocks;
}

function readCategories(state: MsgPackState) {
    const categories: BinaryCifCategoryHeader[] = [];
    const categoryCount = readMsgPackArrayLength(state);
    for (let i = 0; i < categoryCount; i++) {
        const category: BinaryCifCategoryHeader = { name: '', rowCount: 0, columnNames: [] };
        const count = readMsgPackMapLength(state);
        for (let j = 0; j < count; j++) {
            switch (parseMsgPackValue(state)) {
                case 'name': category.name = parseMsgPackValue(state); break;
                case 'rowCount': category.rowCount = parseMsgPackValue(state); break;
                case 'columns': category.columnNames = readColumnNames(state); break;
                default: skipMsgPackValue(state);
            }
        }
        categories.push(category);
    }
    return categories;
}

function readColumnNames(state: MsgPackState) {
    const names: string[] = [];
    const columnCount = readMsgPackArrayLength(state);
    for (let i = 0; i < columnCount; i++) {
        const count = readMsgPackMapLength(state);
        for (let j = 0; j < count; j++) {
            // `data` and `mask` hold the encoded column payloads and are skipped without allocating.
            if (parseMsgPackValue(state) === 'name') names.push(parseMsgPackValue(state));
            else skipMsgPackValue(state);
        }
    }
    return names;
}

/** `categoryName` must include the leading underscore, e.g. `_atom_site`. */
export function binaryCifHasCategory(header: BinaryCifHeader, categoryName: string): boolean {
    for (const block of header.dataBlocks) {
        for (const category of block.categories) {
            if (category.name === categoryName) return true;
        }
    }
    return false;
}

/** `categoryName` must include the leading underscore, e.g. `_atom_site`. */
export function binaryCifHasColumn(header: BinaryCifHeader, categoryName: string, columnName: string): boolean {
    for (const block of header.dataBlocks) {
        for (const category of block.categories) {
            if (category.name !== categoryName) continue;
            for (const name of category.columnNames) {
                if (name === columnName) return true;
            }
        }
    }
    return false;
}
