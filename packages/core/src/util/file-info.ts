/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Russell Parker <russell@benchling.com>
 */

// TODO only support compressed files for which uncompression support is available???
// TODO store globally with decompression plugins?
const COMPRESSED_EXT_LIST = ['gz', 'zip'];

export interface FileNameInfo {
    path: string
    name: string
    ext: string
    base: string
    dir: string
    protocol: string
    query: string
}

export function getFileNameInfo(fileName: string): FileNameInfo {
    let path: string = fileName;
    let protocol = '';

    const queryIndex = path.lastIndexOf('?');
    const query = queryIndex !== -1 ? path.substring(queryIndex) : '';
    path = path.substring(0, queryIndex === -1 ? path.length : queryIndex);

    const name = path.replace(/^.*[\\/]/, '');
    let base = name.substring(0, name.lastIndexOf('.'));

    const nameSplit = name.split('.');
    let ext = nameSplit.length > 1 ? (nameSplit.pop() || '').toLowerCase() : '';

    const protocolMatch = path.match(/^(.+):\/\/(.+)$/);
    if (protocolMatch) {
        protocol = protocolMatch[1].toLowerCase();
        path = protocolMatch[2] || '';
    }

    const dir = path.substring(0, path.lastIndexOf('/') + 1);

    if (COMPRESSED_EXT_LIST.includes(ext)) {
        const n = path.length - ext.length - 1;
        ext = (path.substring(0, n).split('.').pop() || '').toLowerCase();
        const m = base.length - ext.length - 1;
        base = base.substring(0, m);
    }

    // Note: it appears that most of this data never gets used.
    return { path, name, ext, base, dir, protocol, query };
}