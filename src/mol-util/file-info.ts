/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** A File or Blob object or a URL string */
export type FileInput = File | Blob | string

// TODO only support compressed files for which uncompression support is available???
// TODO store globally with decompression plugins?
const compressedExtList = [ 'gz', 'zip' ];

// TODO store globally with parser plugins?
const binaryExtList = [ 'bcif', 'ccp4', 'dcd' ];

export interface FileInfo {
    path: string
    name: string
    ext: string
    base: string
    dir: string
    compressed: string | boolean
    binary: boolean
    protocol: string
    query: string
    src: FileInput
}

export function getFileInfo (file: FileInput): FileInfo {
    let path: string;
    let compressed: string|false;
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
        protocol = protocolMatch[ 1 ].toLowerCase();
        path = protocolMatch[ 2 ] || '';
    }

    const dir = path.substring(0, path.lastIndexOf('/') + 1);

    if (compressedExtList.includes(ext)) {
        compressed = ext;
        const n = path.length - ext.length - 1;
        ext = (path.substr(0, n).split('.').pop() || '').toLowerCase();
        const m = base.length - ext.length - 1;
        base = base.substr(0, m);
    } else {
        compressed = false;
    }

    const binary = binaryExtList.includes(ext);

    return { path, name, ext, base, dir, compressed, binary, protocol, query, src: file };
}