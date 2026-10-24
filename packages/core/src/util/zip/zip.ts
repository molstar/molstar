/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 *
 * - added `ungzip`
 * - added `unzip`'s `only` (one entry read without inflating the rest)
 */

import { writeUint, writeUshort, sizeUTF8, writeUTF8, readUshort, readUint, readUTF8, toInt32 } from './bin.js';
import { crc, adler } from './checksum.js';
import { _inflate } from './inflate.js';
import { _deflateRaw } from './deflate.js';
import { type RuntimeContext, Task } from '@molstar/core/task';

export function Unzip(buf: ArrayBuffer, onlyNames = false) {
  return Task.create('Unzip', (ctx) => unzip(ctx, buf, onlyNames));
}

/**
 * What an archive may hold, checked from its central directory before anything is inflated: at most `maxEntries`
 * entries, their declared sizes summing to at most `maxBytes`, each local entry listed once.
 */
export interface UnzipLimits {
  readonly maxEntries: number;
  readonly maxBytes: number;
}

/** An entry as the central directory declares it (its name, sizes, CRC-32 and the local header's offset). */
interface CentralEntry {
  readonly name: string;
  readonly crc: number;
  readonly csize: number;
  readonly usize: number;
  readonly roff: number;
}

/**
 * The archive's entries by name: their bytes, or with `onlyNames` their sizes. `only`, when given, keeps the entries
 * whose names it accepts and reads nothing of the rest (the limits still cover the whole archive).
 */
export async function unzip(
  runtime: RuntimeContext,
  buf: ArrayBuffer,
  onlyNames = false,
  limits?: UnzipLimits,
  only?: (name: string) => boolean,
) {
  const out: { [k: string]: Uint8Array<ArrayBuffer> | { size: number; csize: number } } = Object.create(null);
  const data = new Uint8Array(buf);
  if (readUshort(data, 0) !== 0x4b50) {
    throw new Error('Invalid ZIP file. A valid ZIP file must start with two magic bytes \\x50\\x4b ("PK" in ASCII).');
  }
  // The end record is in the last 22 bytes, before a comment of at most 65,535; a truncated archive has none
  let eocd = data.length - 22;
  const first = Math.max(0, eocd - 65535);
  while (eocd >= first && readUint(data, eocd) !== 0x06054b50) eocd--;
  if (eocd < first) throw new Error('no end of central directory (the archive is cut short or not a zip)');

  let o = eocd;
  o += 4; // sign  = 0x06054b50
  o += 4; // disks = 0;
  const cnu = readUshort(data, o);
  o += 2;
  // const cnt = readUshort(data, o);
  o += 2;

  // const csize = readUint(data, o);
  o += 4;
  const coffs = readUint(data, o);
  o += 4;

  // The central directory is read whole first and its sizes and CRC are the ones used (the local header's are zeros
  // in an archive written with a data descriptor), and the limits are checked before anything is inflated.
  if (limits && cnu > limits.maxEntries) throw new Error(`${cnu} entries, over the bound of ${limits.maxEntries}`);
  const entries: CentralEntry[] = [];
  const offsets = new Set<number>();
  let total = 0;
  o = coffs;
  for (let i = 0; i < cnu; i++) {
    if (o + 46 > eocd || readUint(data, o) !== 0x02014b50) throw new Error(`central directory entry ${i} is not one`);
    o += 4;
    o += 4; // versions;
    o += 4; // flag + compr
    o += 4; // time

    const crc32 = readUint(data, o);
    o += 4;
    const csize = readUint(data, o);
    o += 4;
    const usize = readUint(data, o);
    o += 4;

    const nl = readUshort(data, o);
    const el = readUshort(data, o + 2);
    const cl = readUshort(data, o + 4);
    o += 6; // name, extra, comment
    o += 8; // disk, attribs

    const roff = readUint(data, o);
    o += 4;
    const name = readUTF8(data, o, nl);
    o += nl + el + cl;

    if (offsets.has(roff)) throw new Error(`"${name}" lists a local entry another lists too`);
    offsets.add(roff);
    total += usize;
    if (limits && total > limits.maxBytes) {
      throw new Error(`the entries inflate to over the bound of ${limits.maxBytes} bytes`);
    }
    entries.push({ name, crc: crc32, csize, usize, roff });
  }
  for (const e of entries) {
    if (!only || only(e.name)) await _readLocal(runtime, data, e, out, onlyNames);
  }
  // console.log(out);
  return out;
}

async function _readLocal(
  runtime: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  entry: CentralEntry,
  out: { [k: string]: Uint8Array | { size: number; csize: number } },
  onlyNames: boolean,
) {
  // The sizes, name and CRC are the central directory's; the local header gives the method and the lengths of
  // its name and extra field.
  const { name, csize, usize } = entry;
  let o = entry.roff;
  if (o + 30 > data.length || readUint(data, o) !== 0x04034b50) throw new Error(`"${name}": no local header`);
  // const sign  = readUint(data, o);
  o += 4;
  // const ver   = readUshort(data, o);
  o += 2;
  // const gpflg = readUshort(data, o);
  o += 2;
  const cmpr = readUshort(data, o);
  o += 2;

  // const time  = readUint(data, o);
  o += 4;

  o += 12; // crc32, csize, usize: the central directory's are used

  const nlen = readUshort(data, o);
  o += 2;
  const elen = readUshort(data, o);
  o += 2;

  o += nlen;
  o += elen;

  if (onlyNames) {
    out[name] = { size: usize, csize };
    return;
  }

  if (o + csize > data.length) throw new Error(`"${name}" runs past the end of the archive`);
  const file = new Uint8Array(data.buffer, o, csize);
  let got: Uint8Array<ArrayBuffer>;
  if (cmpr === 0) {
    got = file;
  } else if (cmpr === 8) {
    got = await inflateRaw(runtime, file, new Uint8Array(usize));
  } else {
    throw new Error(`unknown compression method: ${cmpr}`);
  }
  // The bytes are what was declared and stored: a short stream would leave zeros and a flipped byte read as data
  if (got.length !== usize) throw new Error(`"${name}" holds ${got.length} bytes, declared ${usize}`);
  if (crc(got, 0, got.length) >>> 0 !== entry.crc) {
    throw new Error(`"${name}" fails its CRC-32: the archive is damaged`);
  }
  out[name] = got;
}

export async function inflateRaw(
  runtime: RuntimeContext,
  file: Uint8Array<ArrayBuffer>,
  buf?: Uint8Array<ArrayBuffer>,
) {
  return _inflate(runtime, file, buf);
}

export function inflate(runtime: RuntimeContext, file: Uint8Array<ArrayBuffer>, buf?: Uint8Array<ArrayBuffer>) {
  // const CMF = file[0]
  // const FLG = file[1]
  // const CM = (CMF&15)
  // const CINFO = (CMF>>>4);
  // console.log(CM, CINFO,CMF,FLG);
  return inflateRaw(runtime, new Uint8Array(file.buffer, file.byteOffset + 2, file.length - 6), buf);
}

/** Whether `data` is a gzip stream (RFC 1952: ID1 ID2 and the deflate method, then room for a header and a trailer). */
export function isGzip(data: Uint8Array): boolean {
  return data.length >= 18 && data[0] === 0x1f && data[1] === 0x8b && data[2] === 8;
}

/** The inflated size a gzip stream declares (its trailer's ISIZE: the size modulo 2^32). */
export function gzipDeclaredSize(data: Uint8Array): number {
  return readUint(data, data.length - 4);
}

// https://tools.ietf.org/html/rfc1952
export async function ungzip(runtime: RuntimeContext, file: Uint8Array<ArrayBuffer>, buf?: Uint8Array<ArrayBuffer>) {
  // const id1 = file[0]
  // const id2 = file[1]
  // const cm = file[2]
  const flg = file[3];
  // const mtime = readUint(file, 4)
  // const xfl = file[8]
  // const os = file[9]

  let o = 10;
  if (flg & 4) {
    // FEXTRA
    const xlen = readUshort(file, o);
    // console.log('FEXTRA', xlen)
    o += 2 + xlen; // XLEN itself, then its bytes (RFC 1952 §2.3.1.1)
  }
  if (flg & 8) {
    // FNAME
    let zero = o;
    while (file[zero] !== 0) ++zero;
    // const name = readUTF8(file, o, zero - o)
    // console.log('FNAME', name, zero - o)
    o = zero + 1;
  }
  if (flg & 16) {
    // FCOMMENT
    let zero = o;
    while (file[zero] !== 0) ++zero;
    // const comment = readUTF8(file, o, zero - o)
    // console.log('FCOMMENT', comment)
    o = zero + 1;
  }

  if (flg & 1) {
    // FHCRC
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

export async function deflate(
  runtime: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  opts?: { level: number } /* , buf, off*/,
) {
  if (opts === undefined) opts = { level: 6 };
  let off = 0;
  const buf = new Uint8Array(50 + Math.floor(data.length * 1.1));
  buf[off] = 120;
  buf[off + 1] = 156;
  off += 2;
  off = await _deflateRaw(runtime, data, buf, off, opts.level);
  const crcValue = adler(data, 0, data.length);
  buf[off + 0] = (crcValue >>> 24) & 255;
  buf[off + 1] = (crcValue >>> 16) & 255;
  buf[off + 2] = (crcValue >>> 8) & 255;
  buf[off + 3] = (crcValue >>> 0) & 255;
  return new Uint8Array(buf.buffer, 0, off + 4);
}

async function deflateRaw(runtime: RuntimeContext, data: Uint8Array<ArrayBuffer>, opts?: { level: number }) {
  if (opts === undefined) opts = { level: 6 };
  const buf = new Uint8Array(50 + Math.floor(data.length * 1.1));
  const off = await _deflateRaw(runtime, data, buf, 0, opts.level);
  return new Uint8Array(buf.buffer, 0, off);
}

export function Zip(obj: { [k: string]: Uint8Array<ArrayBuffer> }, noCmpr = false) {
  return Task.create('Zip', (ctx) => zip(ctx, obj, noCmpr));
}

export async function zip(runtime: RuntimeContext, obj: { [k: string]: Uint8Array<ArrayBuffer> }, noCmpr = false) {
  let tot = 0;
  const zpd: { [k: string]: { cpr: boolean; usize: number; crc: number; file: Uint8Array } } = {};
  for (const p in obj) {
    const cpr = !_noNeed(p) && !noCmpr,
      buf = obj[p];
    const crcValue = crc(buf, 0, buf.length);
    zpd[p] = {
      cpr,
      usize: buf.length,
      crc: crcValue,
      file: cpr ? await deflateRaw(runtime, buf) : buf,
    };
  }

  for (const p in zpd) tot += zpd[p].file.length + 30 + 46 + 2 * sizeUTF8(p);
  tot += 22;

  const data = new Uint8Array(tot);
  let o = 0;
  const fof = [];

  for (const p in zpd) {
    const file = zpd[p];
    fof.push(o);
    o = _writeHeader(data, o, p, file, 0);
  }
  let i = 0;
  const ioff = o;
  for (const p in zpd) {
    const file = zpd[p];
    fof.push(o);
    o = _writeHeader(data, o, p, file, 1, fof[i++]);
  }
  const csize = o - ioff;

  writeUint(data, o, 0x06054b50);
  o += 4;
  o += 4; // disks
  writeUshort(data, o, i);
  o += 2;
  writeUshort(data, o, i);
  o += 2; // number of c d records
  writeUint(data, o, csize);
  o += 4;
  writeUint(data, o, ioff);
  o += 4;
  o += 2;
  return data.buffer;
}

// no need to compress .PNG, .ZIP, .JPEG ....
function _noNeed(fn: string) {
  const ext = fn.split('.').pop()!.toLowerCase();
  return 'png,jpg,jpeg,zip'.indexOf(ext) !== -1;
}

function _writeHeader(
  data: Uint8Array,
  o: number,
  p: string,
  obj: { cpr: boolean; usize: number; crc: number; file: Uint8Array },
  t: number,
  roff = 0,
) {
  const file = obj.file;

  writeUint(data, o, t === 0 ? 0x04034b50 : 0x02014b50);
  o += 4; // sign
  if (t === 1) o += 2; // ver made by
  writeUshort(data, o, 20);
  o += 2; // ver
  writeUshort(data, o, 0);
  o += 2; // gflip
  writeUshort(data, o, obj.cpr ? 8 : 0);
  o += 2; // cmpr

  writeUint(data, o, 0);
  o += 4; // time
  writeUint(data, o, obj.crc);
  o += 4; // crc32
  writeUint(data, o, file.length);
  o += 4; // csize
  writeUint(data, o, obj.usize);
  o += 4; // usize

  writeUshort(data, o, sizeUTF8(p));
  o += 2; // nlen
  writeUshort(data, o, 0);
  o += 2; // elen

  if (t === 1) {
    o += 2; // comment length
    o += 2; // disk number
    o += 6; // attributes
    writeUint(data, o, roff);
    o += 4; // usize
  }
  const nlen = writeUTF8(data, o, p);
  o += nlen;
  if (t === 0) {
    data.set(file, o);
    o += file.length;
  }
  return o;
}
