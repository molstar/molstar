/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** A File or Blob object or a URL string */
export type FileInput = File | Blob | string

// TODO only support compressed files for which uncompression support is available???
// TODO store globally with decompression plugins?
const compressedExtList = ['gz', 'zip'];

export interface FileInfo {
    path: string
    name: string
    ext: string
    base: string
    dir: string
    protocol: string
    query: string
    src: FileInput
}

export function getFileInfo(file: FileInput): FileInfo {
    let path: string;
    let protocol = '';

    if (file instanceof File) {
        path = file.name;
    } else if (file instanceof Blob) {
        path = '';
    } else {
        path = file;
    }
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

    if (compressedExtList.includes(ext)) {
        const n = path.length - ext.length - 1;
        ext = (path.substr(0, n).split('.').pop() || '').toLowerCase();
        const m = base.length - ext.length - 1;
        base = base.substr(0, m);
    }

    return { path, name, ext, base, dir, protocol, query, src: file };
}