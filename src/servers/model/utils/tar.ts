/**
 * Adapter from https://github.com/mafintosh/tar-stream
 * Copyright (c) 2014 Mathias Buus, MIT License (MIT)
 */

import { constants } from 'fs';

let alloc = Buffer.alloc;

let ZEROS = '0000000000000000000';
let SEVENS = '7777777777777777777';
let ZERO_OFFSET = '0'.charCodeAt(0);
let USTAR_MAGIC = Buffer.from('ustar\x00', 'binary');
let USTAR_VER = Buffer.from('00', 'binary');
let MASK = parseInt('7777', 8);
let MAGIC_OFFSET = 257;
let VERSION_OFFSET = 263;

let toTypeflag = function (flag: string) {
    switch (flag) {
        case 'file':
            return 0;
        case 'link':
            return 1;
        case 'symlink':
            return 2;
        case 'character-device':
            return 3;
        case 'block-device':
            return 4;
        case 'directory':
            return 5;
        case 'fifo':
            return 6;
        case 'contiguous-file':
            return 7;
        case 'pax-header':
            return 72;
    }

    return 0;
};

let indexOf = function (block: any, num: any, offset: any, end: any) {
    for (; offset < end; offset++) {
        if (block[offset] === num) return offset;
    }
    return end;
};

let cksum = function (block: any) {
    let sum = 8 * 32;
    for (let i = 0; i < 148; i++) sum += block[i];
    for (let j = 156; j < 512; j++) sum += block[j];
    return sum;
};

let encodeOct = function (val: any, n: any) {
    val = val.toString(8);
    if (val.length > n) return SEVENS.slice(0, n) + ' ';
    else return ZEROS.slice(0, n - val.length) + val + ' ';
};

let decodeStr = function (val: any, offset: any, length: any, encoding?: any) {
    return val.slice(offset, indexOf(val, 0, offset, offset + length)).toString(encoding);
};

let addLength = function (str: any) {
    let len = Buffer.byteLength(str);
    let digits = Math.floor(Math.log(len) / Math.log(10)) + 1;
    if (len + digits >= Math.pow(10, digits)) digits++;

    return (len + digits) + str;
};

exports.decodeLongPath = function (buf: any, encoding: any) {
    return decodeStr(buf, 0, buf.length, encoding);
};

exports.encodePax = function (opts: any) {
    let result = '';
    if (opts.name) result += addLength(' path=' + opts.name + '\n');
    if (opts.linkname) result += addLength(' linkpath=' + opts.linkname + '\n');
    let pax = opts.pax;
    if (pax) {
        for (let key in pax) {
            result += addLength(' ' + key + '=' + pax[key] + '\n');
        }
    }
    return Buffer.from(result);
};

exports.decodePax = function (buf: any) {
    let result: any = {};

    while (buf.length) {
        let i = 0;
        while (i < buf.length && buf[i] !== 32) i++;
        let len = parseInt(buf.slice(0, i).toString(), 10);
        if (!len) return result;

        let b = buf.slice(i + 1, len - 1).toString();
        let keyIndex = b.indexOf('=');
        if (keyIndex === -1) return result;
        result[b.slice(0, keyIndex)] = b.slice(keyIndex + 1);

        buf = buf.slice(len);
    }

    return result;
};

export interface Headers {
    name: string;
    mode?: number;
    uid?: number;
    gid?: number;
    size?: number;
    mtime?: Date;
    linkname?: string | null;
    type?:
    | 'file'
    | 'link'
    | 'symlink'
    | 'character-device'
    | 'block-device'
    | 'directory'
    | 'fifo'
    | 'contiguous-file'
    | 'pax-header'
    | 'pax-global-header'
    | 'gnu-long-link-path'
    | 'gnu-long-path'
    | null;
    uname?: string;
    gname?: string;
    devmajor?: number;
    devminor?: number;
    typeflag?: number
}

function modeToType(mode: number) {
    switch (mode & constants.S_IFMT) {
        case constants.S_IFBLK: return 'block-device';
        case constants.S_IFCHR: return 'character-device';
        case constants.S_IFDIR: return 'directory';
        case constants.S_IFIFO: return 'fifo';
        case constants.S_IFLNK: return 'symlink';
    }

    return 'file';
}

let DMODE = parseInt('755', 8);
let FMODE = parseInt('644', 8);

function normalizeHeader(header: Headers) {
    if (!header.size || header.type === 'symlink') header.size = 0;
    if (!header.type) header.type = modeToType(header.mode || 0);
    if (!header.mode) header.mode = header.type === 'directory' ? DMODE : FMODE;
    if (!header.uid) header.uid = 0;
    if (!header.gid) header.gid = 0;
    if (!header.mtime) header.mtime = new Date();
}

export const END_OF_TAR = alloc(1024);

export function encodeTarHeader(opts: Headers) {
    normalizeHeader(opts);

    let buf = alloc(512);
    let name = opts.name;
    let prefix = '';

    if (opts.typeflag === 5 && name[name.length - 1] !== '/') name += '/';
    if (Buffer.byteLength(name) !== name.length) return null; // utf-8

    while (Buffer.byteLength(name) > 100) {
        let i = name.indexOf('/');
        if (i === -1) return null;
        prefix += prefix ? '/' + name.slice(0, i) : name.slice(0, i);
        name = name.slice(i + 1);
    }

    if (Buffer.byteLength(name) > 100 || Buffer.byteLength(prefix) > 155) return null;
    if (opts.linkname && Buffer.byteLength(opts.linkname) > 100) return null;

    buf.write(name);
    buf.write(encodeOct(opts.mode! & MASK, 6), 100);
    buf.write(encodeOct(opts.uid, 6), 108);
    buf.write(encodeOct(opts.gid, 6), 116);
    buf.write(encodeOct(opts.size, 11), 124);
    buf.write(encodeOct((opts.mtime?.getTime()! / 1000) | 0, 11), 136);

    buf[156] = ZERO_OFFSET + toTypeflag(opts.type!);

    if (opts.linkname) buf.write(opts.linkname, 157);

    USTAR_MAGIC.copy(buf, MAGIC_OFFSET);
    USTAR_VER.copy(buf, VERSION_OFFSET);
    if (opts.uname) buf.write(opts.uname, 265);
    if (opts.gname) buf.write(opts.gname, 297);
    buf.write(encodeOct(opts.devmajor || 0, 6), 329);
    buf.write(encodeOct(opts.devminor || 0, 6), 337);

    if (prefix) buf.write(prefix, 345);

    buf.write(encodeOct(cksum(buf), 6), 148);

    return buf;
}