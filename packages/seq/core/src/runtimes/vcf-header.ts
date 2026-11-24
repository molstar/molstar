/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A VCF's header (_spec/README.md "Data", VCF): `##fileformat`, `##contig` (names and lengths, EC-188), `##INFO` and
// `##FORMAT` (`Number`, `Type`) and the sample names of the `#CHROM` line, read within `maxHeaderBytes` and
// `maxSamples` before anything is allocated for them. A header without its `#CHROM` line is malformed.

/** A declaration's `Type` (VCF 4.x). */
export type VcfType = 'Integer' | 'Float' | 'Flag' | 'Character' | 'String';

/** An `##INFO` or `##FORMAT` line: its key, its `Number` as written (`1`, `A`, `R`, `G`, `.`, `2` …) and its `Type`. */
export interface VcfDeclaration {
  readonly id: string;
  readonly number: string;
  readonly type: VcfType;
  readonly description?: string;
}

/** A `##contig` line: its name, and its length when it says one (a whole number). */
export interface VcfContig {
  readonly name: string;
  readonly length?: number;
}

export interface VcfHeader {
  /** `##fileformat`'s value (`VCFv4.2`), when the header has one. */
  readonly fileformat?: string;
  readonly contigs: readonly VcfContig[];
  readonly info: readonly VcfDeclaration[];
  readonly format: readonly VcfDeclaration[];
  /** The `#CHROM` line's columns past `FORMAT`, in order. */
  readonly samples: readonly string[];
  /** The header's lines and bytes (its newlines counted). */
  readonly lines: number;
  readonly bytes: number;
  /** `##` lines that name a structure they do not hold (a declaration without `ID`, a `Type` VCF has not): not read. */
  readonly unread: number;
}

export interface VcfHeaderLimits {
  /** The header's bytes (16 MiB: gnomAD's is a few hundred kB). */
  readonly maxHeaderBytes: number;
  /** Sample names on the `#CHROM` line (100,000). */
  readonly maxSamples: number;
}

export const DefaultVcfHeaderLimits: VcfHeaderLimits = { maxHeaderBytes: 16 * 1024 * 1024, maxSamples: 100_000 };

export type VcfHeaderResult =
  | { readonly ok: true; readonly header: VcfHeader }
  | { readonly ok: false; readonly kind: 'malformed'; readonly message: string };

const Types: ReadonlySet<string> = new Set(['Integer', 'Float', 'Flag', 'Character', 'String']);

/**
 * A structured meta line's fields (`<ID=AF,Number=A,Type=Float,Description="a, \"b\"">`): `key=value` pairs split at
 * commas outside quotes, a quoted value unescaped (`\"`, `\\`); undefined when the text is not one.
 */
export function metaFields(text: string): Map<string, string> | undefined {
  if (!text.startsWith('<') || !text.endsWith('>')) return undefined;
  const out = new Map<string, string>();
  const end = text.length - 1;
  let i = 1;
  while (i < end) {
    const eq = text.indexOf('=', i);
    if (eq < 0 || eq >= end) return undefined;
    const key = text.slice(i, eq).trim();
    let value = '';
    let j = eq + 1;
    if (text[j] === '"') {
      j++;
      while (j < end && text[j] !== '"') {
        if (text[j] === '\\' && j + 1 < end) j++;
        value += text[j++];
      }
      if (text[j] !== '"') return undefined;
      j++;
      while (j < end && text[j] !== ',') j++;
    } else {
      const comma = text.indexOf(',', j);
      const stop = comma < 0 || comma > end ? end : comma;
      value = text.slice(j, stop);
      j = stop;
    }
    if (!out.has(key)) out.set(key, value);
    i = j + 1;
  }
  return out;
}

/** Counts `byte` in `bytes[s, e)`. */
function count(bytes: Uint8Array, s: number, e: number, byte: number): number {
  let n = 0;
  for (let i = bytes.indexOf(byte, s); i >= 0 && i < e; i = bytes.indexOf(byte, i + 1)) n++;
  return n;
}

const decoder = new TextDecoder('utf-8');

/**
 * Reads a VCF's header from its text (the header alone, or the file's head: it stops at the `#CHROM` line), within
 * `limits`: past `maxHeaderBytes` before the `#CHROM` line, or more than `maxSamples` names on it, is `malformed`
 * (counted before the names are split), as is a text without a `#CHROM` line. `##contig`, `##INFO` and `##FORMAT` lines
 * without an `ID`, or with a `Type` VCF has not, are counted (`unread`), the first of an `ID` kept.
 */
export function readVcfHeader(text: Uint8Array | string, limits: Partial<VcfHeaderLimits> = {}): VcfHeaderResult {
  const { maxHeaderBytes, maxSamples } = { ...DefaultVcfHeaderLimits, ...limits };
  const bytes = typeof text === 'string' ? new TextEncoder().encode(text) : text;
  const contigs: VcfContig[] = [];
  const info: VcfDeclaration[] = [];
  const format: VcfDeclaration[] = [];
  const seen = { contig: new Set<string>(), info: new Set<string>(), format: new Set<string>() };
  let fileformat: string | undefined;
  let unread = 0;
  let lines = 0;
  for (let s = 0; s < bytes.length; ) {
    let e = bytes.indexOf(10, s);
    if (e < 0) e = bytes.length;
    const next = e + 1;
    if (e > s && bytes[e - 1] === 13) e--;
    if (e === s) {
      s = next;
      continue;
    }
    if (bytes[s] !== 35) break;
    if (e > maxHeaderBytes)
      return { ok: false, kind: 'malformed', message: `The header is over ${maxHeaderBytes} bytes` };
    lines++;
    if (bytes[s + 1] !== 35) {
      // the #CHROM line: eight fixed columns, then FORMAT and the samples
      const tabs = count(bytes, s, e, 9);
      const samples = Math.max(0, tabs - 8);
      if (samples > maxSamples)
        return { ok: false, kind: 'malformed', message: `${samples} samples, over the bound of ${maxSamples}` };
      const cells = decoder.decode(bytes.subarray(s, e)).split('\t');
      if (cells[0] !== '#CHROM' || cells.length < 8)
        return { ok: false, kind: 'malformed', message: 'The header’s last line is not `#CHROM POS ID REF ALT …`' };
      return {
        ok: true,
        header: {
          ...(fileformat !== undefined && { fileformat }),
          contigs,
          info,
          format,
          samples: cells.slice(9),
          lines,
          bytes: Math.min(next, bytes.length),
          unread,
        },
      };
    }
    const line = decoder.decode(bytes.subarray(s + 2, e));
    const eq = line.indexOf('=');
    const key = eq < 0 ? line : line.slice(0, eq);
    const value = eq < 0 ? '' : line.slice(eq + 1);
    if (key === 'fileformat') fileformat ??= value;
    else if (key === 'contig' || key === 'INFO' || key === 'FORMAT') {
      const f = metaFields(value);
      const id = f?.get('ID');
      if (!f || !id) unread++;
      else if (key === 'contig') {
        if (!seen.contig.has(id)) {
          seen.contig.add(id);
          const length = f.get('length');
          contigs.push(
            length !== undefined && /^\d+$/.test(length) ? { name: id, length: Number(length) } : { name: id },
          );
        }
      } else {
        const type = f.get('Type') ?? '';
        const number = f.get('Number') ?? '';
        const set = key === 'INFO' ? seen.info : seen.format;
        if (!Types.has(type) || !number) unread++;
        else if (!set.has(id)) {
          set.add(id);
          const description = f.get('Description');
          (key === 'INFO' ? info : format).push({
            id,
            number,
            type: type as VcfType,
            ...(description !== undefined && { description }),
          });
        }
      }
    }
    s = next;
  }
  return { ok: false, kind: 'malformed', message: 'No `#CHROM` line: the header is not a VCF’s' };
}
