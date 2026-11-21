/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 *
 * - added `ungzip` (every member of a multi-member stream)
 * - added `unzip`'s `only` (one entry read without inflating the rest)
 */

import { writeUint, writeUshort, sizeUTF8, writeUTF8, readUshort, readUint, readUTF8 } from './bin.js';
import { crc, adler } from './checksum.js';
import { _inflate, _inflateMember, inflateStream } from './inflate.js';
import { checkCompressionStreamSupport } from './util.js';
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
/**
 * The inflated bytes of a gzip stream: every member of it, joined (several members: BGZF's blocks, `cat a.gz b.gz`).
 * Into `buf` when given, never past it. A member cut short, or bytes after a member that are not one, throw.
 */
export async function ungzip(runtime: RuntimeContext, file: Uint8Array<ArrayBuffer>, buf?: Uint8Array<ArrayBuffer>) {
  // The platform's gzip decoder checks the CRC-32 and the size itself, where the path below runs a CRC-32 in JS over
  // the output (perf-audit-2 01/P1). It refuses several members (trailing junk), so a BGZF file, told by its first
  // block's size, goes to the members' path at once; a file it refuses otherwise (plain members, a wrong CRC,
  // truncated) takes that path too.
  const first = gzipMemberHeader(file, 0);
  const several = first?.bsize !== undefined && first.bsize + 1 < file.length;
  if (!several && checkCompressionStreamSupport('gzip')) {
    try {
      return await inflateStream(runtime, file, 'gzip', buf ?? new Uint8Array(gzipDeclaredSize(file)));
    } catch (e) {
      if (Task.isAbort(e)) throw e;
    }
  }
  return ungzipMembers(runtime, file, buf);
}

/**
 * Deflate shrinks data at most about 1032 times: a member declaring more is corrupt or crafted, refused before its
 * buffer is allocated (as `util/_spec/assets.md` §5 does for a whole stream).
 */
const MaxDeflateRatio = 1032;

interface GzipMemberHeader {
  /** Where the deflate stream starts. */
  readonly data: number;
  /** BGZF's `BC` subfield: the member's size minus 1, when present. */
  readonly bsize?: number;
}

/** A gzip member's header at `o` (RFC 1952 §2.3), bounds-checked; undefined when there is none. */
function gzipMemberHeader(file: Uint8Array, o: number): GzipMemberHeader | undefined {
  if (o + 10 > file.length || file[o] !== 0x1f || file[o + 1] !== 0x8b || file[o + 2] !== 8) return;
  const flg = file[o + 3];
  let at = o + 10;
  let bsize: number | undefined;
  if (flg & 4) {
    // FEXTRA: XLEN, then subfields (SI1 SI2 LEN data); BGZF's is 'BC' with a 2-byte BSIZE
    if (at + 2 > file.length) return;
    const xlen = readUshort(file, at);
    const end = at + 2 + xlen;
    if (end > file.length) return;
    for (let s = at + 2; s + 4 <= end; ) {
      const len = readUshort(file, s + 2);
      if (file[s] === 66 && file[s + 1] === 67 && len === 2 && s + 6 <= end) bsize = readUshort(file, s + 4);
      s += 4 + len;
    }
    at = end;
  }
  for (const bit of [8, 16]) {
    // FNAME, FCOMMENT: zero-terminated
    if (!(flg & bit)) continue;
    while (at < file.length && file[at] !== 0) at++;
    if (at >= file.length) return;
    at++;
  }
  if (flg & 2) at += 2; // FHCRC
  if (at > file.length) return;
  // A BSIZE too small for the header and the trailer is not BGZF's: the member is found by its deflate stream instead
  if (bsize !== undefined && o + bsize + 1 < at + 8) bsize = undefined;
  return { data: at, bsize };
}

/**
 * Each member found (by BGZF's BSIZE when present, else where the JS inflater's deflate stream of the member ends),
 * inflated as one member with its own ISIZE, the outputs joined within the caller's buffer, never by growing one.
 */
async function ungzipMembers(runtime: RuntimeContext, file: Uint8Array<ArrayBuffer>, buf?: Uint8Array<ArrayBuffer>) {
  // When every member has its BSIZE (BGZF), the whole size is the sum of their ISIZEs, known before inflating
  if (!buf) {
    let total = 0;
    let o = 0;
    while (o < file.length) {
      const h = gzipMemberHeader(file, o);
      if (h?.bsize === undefined || o + h.bsize + 1 > file.length) break;
      const end = o + h.bsize + 1;
      total += readUint(file, end - 4);
      o = end;
    }
    if (o === file.length && o > 0) {
      if (total > file.length * MaxDeflateRatio + 1024) {
        throw new Error(`gzip members declare ${total} bytes inflated, more than ${file.length} bytes can hold`);
      }
      buf = new Uint8Array(total);
    }
  }
  const parts: Uint8Array<ArrayBuffer>[] = [];
  let total = 0;
  let o = 0;
  const only = (end: number) => o === 0 && end === file.length;
  while (o < file.length) {
    const h = gzipMemberHeader(file, o);
    if (!h) throw new Error(`no gzip member at byte ${o}`);
    const into = (size?: number) =>
      buf
        ? buf.subarray(total, size === undefined ? buf.length : total + size)
        : size === undefined
          ? undefined
          : new Uint8Array(size);
    let bytes: Uint8Array<ArrayBuffer>;
    let end: number;
    let crc32: number;
    if (h.bsize !== undefined) {
      end = o + h.bsize + 1;
      if (end > file.length || end - 8 < h.data) throw new Error(`the gzip member at byte ${o} is cut short`);
      const isize = readUint(file, end - 4);
      if (isize > (end - o) * MaxDeflateRatio + 1024) {
        throw new Error(`the gzip member at byte ${o} declares ${isize} bytes, more than it can hold`);
      }
      if (buf && total + isize > buf.length) throw new Error(`inflates past its declared ${buf.length} bytes`);
      crc32 = readUint(file, end - 8);
      if (checkCompressionStreamSupport('gzip')) {
        // the platform checks the member's CRC-32 and ISIZE
        bytes = await inflateStream(runtime, file.subarray(o, end), 'gzip', into(isize));
        if (bytes.length !== isize) throw new Error(`the gzip member at byte ${o} is cut short`);
        parts.push(bytes);
        total += bytes.length;
        o = end;
        continue;
      }
      bytes = (await _inflateMember(runtime, file.subarray(h.data, end - 8), into(isize))).bytes;
    } else {
      const r = await _inflateMember(runtime, file.subarray(h.data), into());
      bytes = r.bytes;
      end = h.data + r.used + 8;
      if (end > file.length) throw new Error(`the gzip member at byte ${o} is cut short`);
      crc32 = readUint(file, end - 8);
    }
    const isize = readUint(file, end - 4);
    const crcValue = crc(bytes, 0, bytes.length) >>> 0;
    if (crcValue !== crc32 || bytes.length % 2 ** 32 !== isize) {
      // One member read as it always was: a wrong CRC-32 logged, its bytes given
      if (!only(end)) throw new Error(`the gzip member at byte ${o} fails its CRC-32 or size: the file is damaged`);
      console.error("ungzip: checksums don't match");
    }
    parts.push(bytes);
    total += bytes.length;
    o = end;
  }
  if (buf) return total < buf.length ? buf.subarray(0, total) : buf;
  if (parts.length === 1) return parts[0];
  const out = new Uint8Array(total);
  for (let i = 0, at = 0; i < parts.length; at += parts[i].length, i++) out.set(parts[i], at);
  return out;
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

/**
 * `noCmpr`: no entry compressed (true), or the entries named (stored as they are: data that compresses little, whose
 * deflate would cost a buffer of its size and time for nothing)
 */
export function Zip(obj: { [k: string]: Uint8Array<ArrayBuffer> }, noCmpr: boolean | ReadonlySet<string> = false) {
  return Task.create('Zip', (ctx) => zip(ctx, obj, noCmpr));
}

export async function zip(
  runtime: RuntimeContext,
  obj: { [k: string]: Uint8Array<ArrayBuffer> },
  noCmpr: boolean | ReadonlySet<string> = false,
) {
  let tot = 0;
  const zpd: { [k: string]: { cpr: boolean; usize: number; crc: number; file: Uint8Array } } = {};
  for (const p in obj) {
    const cpr = !_noNeed(p) && (typeof noCmpr === 'boolean' ? !noCmpr : !noCmpr.has(p)),
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
