/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 *
 * - added `ungzip`
 */

import { writeUint, writeUshort, sizeUTF8, writeUTF8, readUshort, readUint, readUTF8, toInt32 } from './bin';
import { crc, adler } from './checksum';
import { _inflate } from './inflate';
import { _deflateRaw } from './deflate';
import { RuntimeContext } from '../../mol-task';

export async function unzip(runtime: RuntimeContext, buf: ArrayBuffer, onlyNames = false) {
    const out: { [k: string]: Uint8Array | { size: number, csize: number } } = Object.create(null);
    const data = new Uint8Array(buf);
    let eocd = data.length - 4;

    while(readUint(data, eocd) !== 0x06054b50) eocd--;

    let o = eocd;
    o += 4;	// sign  = 0x06054b50
    o += 4; // disks = 0;
    const cnu = readUshort(data, o);
    o += 2;
    // const cnt = readUshort(data, o);
    o += 2;

    // const csize = readUint(data, o);
    o += 4;
    const coffs = readUint(data, o);  o += 4;

    o = coffs;
    for(let i = 0; i < cnu; i++) {
        // const sign = readUint(data, o);
        o += 4;
        o += 4;  // versions;
        o += 4;  // flag + compr
        o += 4;  // time

        // const crc32 = readUint(data, o);
        o += 4;
        const csize = readUint(data, o);
        o += 4;
        const usize = readUint(data, o);
        o += 4;

        const nl = readUshort(data, o);
        const el = readUshort(data, o + 2);
        const cl = readUshort(data, o + 4);
        o += 6;  // name, extra, comment
        o += 8;  // disk, attribs

        const roff = readUint(data, o);  o += 4;
        o += nl + el + cl;

        await _readLocal(runtime, data, roff, out, csize, usize, onlyNames);
    }
    // console.log(out);
    return out;
}

async function _readLocal(runtime: RuntimeContext, data: Uint8Array, o: number, out: { [k: string]: Uint8Array | { size: number, csize: number } }, csize: number, usize: number, onlyNames: boolean) {
    // const sign  = readUint(data, o);
    o += 4;
    // const ver   = readUshort(data, o);
    o += 2;
    // const gpflg = readUshort(data, o);
    o += 2;
    // if((gpflg&8)!=0) throw "unknown sizes";
    const cmpr  = readUshort(data, o);
    o += 2;

    // const time  = readUint(data, o);
    o += 4;

    // const crc32 = readUint(data, o);
    o += 4;
    // var csize = rUi(data, o);  o+=4;
    // var usize = rUi(data, o);  o+=4;
    o += 8;

    const nlen = readUshort(data, o);
    o += 2;
    const elen = readUshort(data, o);
    o += 2;

    const name = readUTF8(data, o, nlen);
    o += nlen;  // console.log(name);
    o += elen;

    if(onlyNames) {
        out[name] = { size: usize, csize };
        return;
    }

    const file = new Uint8Array(data.buffer, o);
    if(cmpr === 0) {
        out[name] = new Uint8Array(file.buffer.slice(o, o + csize));
    } else if(cmpr === 8) {
        const buf = new Uint8Array(usize);
        await inflateRaw(runtime, file, buf);
        out[name] = buf;
    } else {
        throw `unknown compression method: ${cmpr}`;
    }
}

export async function inflateRaw(runtime: RuntimeContext, file: Uint8Array, buf?: Uint8Array) {
    return _inflate(runtime, file, buf);
}

export function inflate(runtime: RuntimeContext, file: Uint8Array, buf?: Uint8Array) {
    // const CMF = file[0]
    // const FLG = file[1]
    // const CM = (CMF&15)
    // const CINFO = (CMF>>>4);
    // console.log(CM, CINFO,CMF,FLG);
    return inflateRaw(runtime, new Uint8Array(file.buffer, file.byteOffset + 2, file.length - 6), buf);
}

// https://tools.ietf.org/html/rfc1952
export async function ungzip(runtime: RuntimeContext, file: Uint8Array, buf?: Uint8Array) {
    // const id1 = file[0]
    // const id2 = file[1]
    // const cm = file[2]
    const flg = file[3];
    // const mtime = readUint(file, 4)
    // const xfl = file[8]
    // const os = file[9]

    let o = 10;
    if (flg & 4) { // FEXTRA
        const xlen = readUshort(file, o);
        // console.log('FEXTRA', xlen)
        o += xlen;
    }
    if (flg & 8) { // FNAME
        let zero = o;
        while(file[zero] !== 0) ++zero;
        // const name = readUTF8(file, o, zero - o)
        // console.log('FNAME', name, zero - o)
        o = zero + 1;
    }
    if (flg & 16) { // FCOMMENT
        let zero = o;
        while(file[zero] !== 0) ++zero;
        // const comment = readUTF8(file, o, zero - o)
        // console.log('FCOMMENT', comment)
        o = zero + 1;
    }

    if (flg & 1) { // FHCRC
        // const hcrc = readUshort(file, o)
        // console.log('FHCRC', hcrc)
        o += 2;
    }

    const crc32 = toInt32(readUint(file, file.length - 8));
    const isize = readUint(file, file.length - 4);
    if (buf === undefined) buf = new Uint8Array(isize);

    const blocks = new Uint8Array(file.buffer, file.byteOffset + o, file.length - o - 8);
    const inflated = await inflateRaw(runtime, blocks, buf);
    const crcValue = crc(inflated, 0, inflated.length);
    if (crc32 !== crcValue) {
        console.error("ungzip: checksums don't match");
    }

    return inflated;
}

export function deflate(data: Uint8Array, opts?: { level: number }/* , buf, off*/) {
    if(opts === undefined) opts = { level: 6 };
    let off = 0;
    const buf = new Uint8Array(50 + Math.floor(data.length * 1.1));
    buf[off] = 120;  buf[off + 1] = 156;  off += 2;
    off = _deflateRaw(data, buf, off, opts.level);
    const crcValue = adler(data, 0, data.length);
    buf[off + 0] = ((crcValue >>> 24) & 255);
    buf[off + 1] = ((crcValue >>> 16) & 255);
    buf[off + 2] = ((crcValue >>> 8) & 255);
    buf[off + 3] = ((crcValue >>> 0) & 255);
    return new Uint8Array(buf.buffer, 0, off + 4);
}

function deflateRaw(data: Uint8Array, opts?: { level: number }) {
    if(opts === undefined) opts = { level: 6 };
    const buf = new Uint8Array(50 + Math.floor(data.length * 1.1));
    const off = _deflateRaw(data, buf, 0, opts.level);
    return new Uint8Array(buf.buffer, 0, off);
}

export function zip(obj: { [k: string]: Uint8Array }, noCmpr = false) {
    let tot = 0;
    const zpd: { [k: string]: { cpr: boolean, usize: number, crc: number, file: Uint8Array } } = {};
    for(const p in obj) {
        const cpr = !_noNeed(p) && !noCmpr, buf = obj[p];
        const crcValue = crc(buf, 0, buf.length);
        zpd[p] = {
            cpr,
            usize: buf.length,
            crc: crcValue,
            file: (cpr ? deflateRaw(buf) : buf)
        };
    }

    for(const p in zpd) tot += zpd[p].file.length + 30 + 46 + 2 * sizeUTF8(p);
    tot +=  22;

    const data = new Uint8Array(tot);
    let o = 0;
    const fof = [];

    for(const p in zpd) {
        const file = zpd[p];  fof.push(o);
        o = _writeHeader(data, o, p, file, 0);
    }
    let i = 0, ioff = o;
    for(const p in zpd) {
        const file = zpd[p];
        fof.push(o);
        o = _writeHeader(data, o, p, file, 1, fof[i++]);
    }
    const csize = o - ioff;

    writeUint(data, o, 0x06054b50);  o += 4;
    o += 4;  // disks
    writeUshort(data, o, i);  o += 2;
    writeUshort(data, o, i);  o += 2;	// number of c d records
    writeUint(data, o, csize);  o += 4;
    writeUint(data, o, ioff );  o += 4;
    o += 2;
    return data.buffer;
}

// no need to compress .PNG, .ZIP, .JPEG ....
function _noNeed(fn: string) {
    const ext = fn.split('.').pop()!.toLowerCase();
    return 'png,jpg,jpeg,zip'.indexOf(ext) !== -1;
}

function _writeHeader(data: Uint8Array, o: number, p: string, obj: { cpr: boolean, usize: number, crc: number, file: Uint8Array }, t: number, roff = 0) {
    const file = obj.file;

    writeUint(data, o, t === 0 ? 0x04034b50 : 0x02014b50); o += 4; // sign
    if(t === 1) o += 2;  // ver made by
    writeUshort(data, o, 20);  o += 2;	// ver
    writeUshort(data, o,  0);  o += 2;    // gflip
    writeUshort(data, o,  obj.cpr ? 8 : 0);  o += 2;	// cmpr

    writeUint(data, o,  0);  o += 4;	// time
    writeUint(data, o, obj.crc);  o += 4;	// crc32
    writeUint(data, o, file.length);  o += 4;	// csize
    writeUint(data, o, obj.usize);  o += 4;	// usize

    writeUshort(data, o, sizeUTF8(p));  o += 2;	// nlen
    writeUshort(data, o, 0);  o += 2;	// elen

    if(t === 1) {
        o += 2;  // comment length
        o += 2;  // disk number
        o += 6;  // attributes
        writeUint(data, o, roff);  o += 4;	// usize
    }
    const nlen = writeUTF8(data, o, p);  o += nlen;
    if(t === 0) {
        data.set(file, o);
        o += file.length;
    }
    return o;
}