/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ungzip } from '@molstar/core/util/zip/zip';

// A region's lines of a bgzipped file read through its tabix index, as a genome browser reads it (the TBI's bins and
// linear index give the chunks; only their BGZF blocks are inflated): the dataset's test checks its files by region
// this way, never by inflating them whole. A test's reader, independent of the runtimes it checks.

const gunzip = async (b: Uint8Array) => ungzip(SyncRuntimeContext, b as Uint8Array<ArrayBuffer>);

/** The bins that may hold `[beg, end)` (htslib's `reg2bins`, 14 bits and 5 levels). */
function regionBins(beg: number, end: number): number[] {
  const out = [0];
  const e = end - 1;
  for (const [first, shift] of [
    [1, 26],
    [9, 23],
    [73, 20],
    [585, 17],
    [4681, 14],
  ] as const)
    for (let k = first + (beg >> shift); k <= first + (e >> shift); k++) out.push(k);
  return out;
}

/** A TBI's reference names and, per reference, its bins' chunks (virtual offsets) and linear index. */
interface TabixIndex {
  readonly names: string[];
  readonly refs: { bins: Map<number, [bigint, bigint][]>; linear: bigint[] }[];
}

async function readIndex(tbi: Uint8Array): Promise<TabixIndex> {
  const b = await gunzip(tbi);
  const d = new DataView(b.buffer, b.byteOffset, b.byteLength);
  if (String.fromCharCode(b[0], b[1], b[2], b[3]) !== 'TBI\u0001') throw new Error('not a TBI');
  const nRef = d.getInt32(4, true);
  const lNm = d.getInt32(32, true);
  const names = new TextDecoder()
    .decode(b.subarray(36, 36 + lNm))
    .split('\u0000')
    .filter(Boolean);
  let o = 36 + lNm;
  const refs: TabixIndex['refs'] = [];
  for (let r = 0; r < nRef; r++) {
    const bins = new Map<number, [bigint, bigint][]>();
    const nBin = d.getInt32(o, true);
    o += 4;
    for (let i = 0; i < nBin; i++) {
      const bin = d.getUint32(o, true);
      const nChunk = d.getInt32(o + 4, true);
      o += 8;
      const chunks: [bigint, bigint][] = [];
      for (let c = 0; c < nChunk; c++, o += 16) chunks.push([d.getBigUint64(o, true), d.getBigUint64(o + 8, true)]);
      bins.set(bin, chunks);
    }
    const nIntv = d.getInt32(o, true);
    o += 4;
    const linear: bigint[] = [];
    for (let i = 0; i < nIntv; i++, o += 8) linear.push(d.getBigUint64(o, true));
    refs.push({ bins, linear });
  }
  return { names, refs };
}

/** The text from virtual offset `from` to `to`: the BGZF blocks between them inflated, one at a time. */
async function textBetween(file: Uint8Array, from: bigint, to: bigint): Promise<string> {
  const c0 = Number(from >> 16n);
  const c1 = Number(to >> 16n);
  const parts: Uint8Array[] = [];
  for (let c = c0; c <= c1 && c < file.length; ) {
    // BSIZE (the block's size less one) in the `BC` extra field at byte 16
    const size = (file[c + 16] | (file[c + 17] << 8)) + 1;
    const block = await gunzip(file.subarray(c, c + size));
    const u0 = c === c0 ? Number(from & 0xffffn) : 0;
    const u1 = c === c1 ? Number(to & 0xffffn) : block.length;
    parts.push(block.subarray(u0, u1));
    c += size;
  }
  return parts.map((p) => new TextDecoder().decode(p)).join('');
}

/**
 * A region's lines (each whose extent meets `[start, end)`, `extent` giving it from a line's columns) of a bgzipped
 * `file` through its TBI `index`; the blocks inflated, for the test to say it read the file by region.
 */
export async function tabixRegion(
  file: Uint8Array,
  index: Uint8Array,
  name: string,
  start: number,
  end: number,
  extent: (columns: string[]) => [number, number],
): Promise<{ lines: string[]; inflated: number }> {
  const tbi = await readIndex(index);
  const ref = tbi.refs[tbi.names.indexOf(name)];
  if (!ref) return { lines: [], inflated: 0 };
  const min = ref.linear[Math.min(start >> 14, ref.linear.length - 1)] ?? 0n;
  const chunks = regionBins(start, end)
    .flatMap((bin) => ref.bins.get(bin) ?? [])
    .filter(([, e]) => e > min)
    .sort((a, b) => (a[0] < b[0] ? -1 : 1));
  const seen = new Set<string>();
  const lines: string[] = [];
  let inflated = 0;
  for (const [s, e] of chunks) {
    const text = await textBetween(file, s > min ? s : min, e);
    inflated += text.length;
    for (const l of text.split('\n')) {
      const c = l.split('\t');
      if (!l || l.startsWith('#') || c[0] !== name || seen.has(l)) continue;
      const [b, x] = extent(c);
      if (b < end && x > start) {
        seen.add(l);
        lines.push(l);
      }
    }
  }
  return { lines, inflated };
}
