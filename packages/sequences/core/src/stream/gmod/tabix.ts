/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tabix-indexed VCF: items (README, "Tabix and BGZF"). A BGZF file is gzip members that each hold at most 64 kb, so a
// region is read by the index's chunks (virtual offsets: block file offset << 16 | offset in its inflated bytes).

import type { StreamProvider } from '@molstar/seq-core/format';
import { ITEM_LEVELS, type ItemsRead, itemsProvider } from './items.js';
import { inflate, once, type RangeReader, text, view } from './ranges.js';

const MAX_BLOCK = 65536;
const TYPES = ['SNV', 'MNV', 'insertion', 'deletion', 'structural'];
const SHIFT = 65536;

/** The size of the BGZF block at `at`: its `BC` extra subfield's BSIZE + 1; 0 where fewer than its header is held. */
function blockSize(bytes: Uint8Array, at: number): number {
  if (at + 18 > bytes.length) return 0;
  const dv = view(bytes);
  for (let x = at + 12, end = x + dv.getUint16(at + 10, true); x < end; x += 4 + dv.getUint16(x + 2, true))
    if (bytes[x] === 66 && bytes[x + 1] === 67) return dv.getUint16(x + 4, true) + 1;
  throw new Error('not a BGZF file: a block has no BC field');
}
/** What a whole index may inflate to (a large VCF's `.tbi` is tens of MB). */
const MAX_INDEX = 256 << 20;

/**
 * The complete blocks of a BGZF run inflated, each with the offset it started at (a read may cut the last). Each
 * block holds at most 64 KB and exactly its member's ISIZE, and the run at most `max` bytes in all (GI-07).
 */
async function bgzf(bytes: Uint8Array, max = MAX_INDEX) {
  const out: { at: number; size: number; data: Uint8Array }[] = [];
  const dv = view(bytes);
  let total = 0;
  for (
    let at = 0, size = blockSize(bytes, 0);
    size > 0 && at + size <= bytes.length;
    at += size, size = blockSize(bytes, at)
  ) {
    const data = await inflate(bytes.subarray(at, at + size), 'gzip', MAX_BLOCK);
    // A second line: the gzip stream checks the trailer itself (LA-34).
    if (data.length !== dv.getUint32(at + size - 4, true))
      throw new Error(`BGZF: the block at ${at} inflates to ${data.length} bytes, not its ISIZE`);
    total += data.length;
    if (total > max) throw new Error(`BGZF: the run inflates past ${max} bytes`);
    out.push({ at, size, data });
  }
  return out;
}
const joined = (parts: Uint8Array[]) => {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) {
    out.set(p, at);
    at += p.length;
  }
  return out;
};

interface Ref {
  readonly name: string;
  readonly bins: Map<number, number[]>; // bin -> chunk (begin, end) virtual offsets, flat
  readonly linear: number[];
}
/** The `.tbi` (BGZF, then its own layout; all little-endian): a VCF's `names` and per reference the bins and linear index. */
async function parseTbi(file: Uint8Array): Promise<Ref[]> {
  const b = joined((await bgzf(file)).map((x) => x.data));
  const dv = view(b);
  if (dv.getUint32(0, true) !== 0x01494254) throw new Error('not a tabix index: no TBI\\1 magic');
  if ((dv.getInt32(8, true) & 0xffff) !== 2) throw new Error('tabix: only VCF-format indexes are read');
  const names = text(b.subarray(36, 36 + dv.getInt32(32, true)))
    .split('\0')
    .filter(Boolean);
  const voff = (at: number) => dv.getUint32(at + 4, true) * 2 ** 32 + dv.getUint32(at, true);
  let at = 36 + dv.getInt32(32, true);
  return names.map((name) => {
    const bins = new Map<number, number[]>();
    const count = dv.getInt32(at, true);
    at += 4;
    for (let i = 0; i < count; i++) {
      const [bin, chunks] = [dv.getUint32(at, true), dv.getInt32(at + 4, true)];
      at += 8;
      bins.set(
        bin,
        Array.from({ length: 2 * chunks }, (_, k) => voff(at + 8 * k)),
      );
      at += 16 * chunks;
    }
    const linear = Array.from({ length: dv.getInt32(at, true) }, (_, k) => voff(at + 4 + 8 * k));
    at += 4 + 8 * linear.length;
    return { name, bins, linear };
  });
}

/** The bins that may hold a record meeting [beg, end) (the SAM specification's `reg2bins`). */
function bins(beg: number, end: number): number[] {
  const out = [0];
  for (const [first, shift] of [
    [1, 26],
    [9, 23],
    [73, 20],
    [585, 17],
    [4681, 14],
  ])
    for (let k = first + (beg >> shift); k <= first + ((end - 1) >> shift); k++) out.push(k);
  return out;
}

/**
 * Compressed bytes a record takes in its BGZF blocks until a file's own are measured (`measured`): ClinVar's, about
 * (measured on the recorded slices; ST4b).
 */
const BYTES_PER_ITEM = 28;
/** Records a measure needs before it replaces `BYTES_PER_ITEM`. */
const MEASURE_MIN = 50;

/**
 * Compressed bytes per record over the blocks that hold records only (no header line): a block's size over its lines,
 * which a block cut mid-line miscounts by one; 0: too few records to say.
 */
function measured(blocks: { size: number; data: Uint8Array }[]): number {
  let [bytes, records] = [0, 0];
  for (const b of blocks) {
    const t = text(b.data);
    if (/(^|\n)#/.test(t)) continue;
    bytes += b.size;
    for (let at = t.indexOf('\n'); at >= 0; at = t.indexOf('\n', at + 1)) records++;
  }
  return records >= MEASURE_MIN ? bytes / records : 0;
}

/** The chunks of the bins meeting [start, stop), less those before the linear index's smallest offset, merged. */
function chunksOf(ref: Ref, start: number, stop: number): [number, number][] {
  const floor = ref.linear[Math.min(start >> 14, ref.linear.length - 1)] ?? 0;
  const chunks: [number, number][] = [];
  for (const bin of bins(start, Math.max(stop, start + 1))) {
    const c = ref.bins.get(bin) ?? [];
    for (let i = 0; i < c.length; i += 2) if (c[i + 1] > floor) chunks.push([Math.max(c[i], floor), c[i + 1]]);
  }
  chunks.sort((a, b) => a[0] - b[0]);
  return chunks.reduce<[number, number][]>((m, c) => {
    if (m.length && c[0] <= m[m.length - 1][1]) m[m.length - 1][1] = Math.max(m[m.length - 1][1], c[1]);
    else m.push([...c]);
    return m;
  }, []);
}

/**
 * Compressed bytes per bin of `w` positions over [start, stop), from the linear index: a 16 kb window's bytes are the
 * distance between its first record's block and the next window's, spread evenly over the window. No data is read.
 */
function windowBytes(ref: Ref, start: number, stop: number, w: number): Float64Array {
  const out = new Float64Array(Math.max(1, Math.ceil((stop - start) / w)));
  const L = ref.linear;
  // The last window ends where the reference's last chunk does.
  let end = 0;
  for (const c of ref.bins.values()) for (let i = 1; i < c.length; i += 2) end = Math.max(end, c[i]);
  for (let k = start >> 14; k <= (stop - 1) >> 14 && k < L.length; k++) {
    const next = k + 1 < L.length ? Math.max(L[k + 1], L[k]) : end;
    const bytes = Math.max(0, Math.floor(next / SHIFT) - Math.floor(L[k] / SHIFT));
    const [a, e] = [Math.max(start, k << 14), Math.min(stop, (k + 1) << 14)];
    for (let j = Math.floor((a - start) / w); j < out.length && start + j * w < e; j++) {
      const over = Math.min(e, start + (j + 1) * w) - Math.max(a, start + j * w);
      if (over > 0) out[j] += (bytes * over) / 16384;
    }
  }
  return out;
}

/**
 * A VCF's `.tbi`: its bytes, or a function that fetches them under the open's signal (aborted when every request
 * waiting on the open has gone), called again when an open fails.
 */
export type VcfIndex = Uint8Array | Promise<Uint8Array> | ((signal?: AbortSignal) => Promise<Uint8Array>);

/** A tabix-indexed VCF as an `items` provider (`index`: the `.tbi` file, whole); `info` names INFO keys kept as columns. */
export function vcfProvider(
  data: RangeReader,
  index: VcfIndex,
  o: { info?: string[]; maxItems?: number; lengths?: Record<string, number> } = {},
): StreamProvider {
  // The estimates' compressed bytes a record: the fixed figure until a data block is read (the header's, or a region's).
  let perItem = BYTES_PER_ITEM;
  let learnt = false;
  const learn = (blocks: { size: number; data: Uint8Array }[]) => {
    const m = learnt ? 0 : measured(blocks);
    if (m > 0) [perItem, learnt] = [m, true];
  };
  const open = once(async (signal: AbortSignal) => {
    const refs = new Map(
      (await parseTbi(await (typeof index === 'function' ? index(signal) : index))).map((r) => [r.name, r]),
    );
    // Lengths are not in the index: they are the header's `##contig` lines, in the first blocks.
    let header = '';
    let blocks: Awaited<ReturnType<typeof bgzf>> = [];
    for (let n = 16384; !/\n#CHROM[^\n]*\n/.test(header) && n <= 64 * MAX_BLOCK; n *= 4) {
      blocks = await bgzf(await data(0, n, signal));
      header = text(joined(blocks.map((x) => x.data)));
    }
    // The first blocks past the header, read already, measure the file's bytes a record.
    learn(blocks);
    const lengths = new Map<string, number>();
    for (const m of header.matchAll(/##contig=<ID=([^,>]+),length=(\d+)/g)) lengths.set(m[1], Number(m[2]));
    return { refs, lengths };
  });
  return itemsProvider(
    {
      vocabulary: { category: TYPES },
      estimated: true,
      async systems(signal) {
        const { refs, lengths } = await open(signal);
        // The host's `lengths`, else the header's `##contig` lines, else where the linear index (16 kb windows) ends.
        return [...refs.values()].map((r) => ({
          key: r.name,
          extent: o.lengths?.[r.name] ?? lengths.get(r.name) ?? r.linear.length * 16384,
        }));
      },
      // ST4b: what a region would cost, from the index alone, as JBrowse 2 does before it fetches: its chunks'
      // compressed bytes, and its items from the linear index's 16 kb windows (compressed bytes over the file's bytes a record).
      async estimate(system, start, stop, signal) {
        const ref = (await open(signal)).refs.get(system);
        if (!ref) return null;
        let bytes = 0;
        for (const [from, to] of chunksOf(ref, start, stop)) bytes += Math.floor(to / SHIFT) - Math.floor(from / SHIFT);
        let items = 0;
        for (const x of windowBytes(ref, start, stop, stop - start)) items += x / perItem;
        return { items: Math.round(items), bytes };
      },
      async bins(system, start, stop, level, signal) {
        const ref = (await open(signal)).refs.get(system);
        if (!ref) return null;
        const w = ITEM_LEVELS[level];
        return windowBytes(ref, start, stop, w).map((x) => Math.round(x / perItem));
      },
      async read(system, start, stop, signal): Promise<ItemsRead> {
        const ref = (await open(signal)).refs.get(system);
        if (!ref) throw new Error(`VCF: no chromosome ${system} in the index`);
        const merged = chunksOf(ref, start, stop);
        const t: Record<string, (string | number)[]> = {};
        for (const k of ['start', 'stop', 'id', 'label', 'ref', 'alt', 'category', 'name', ...(o.info ?? [])])
          t[k] = [];
        const [starts, stops]: number[][] = [[], []];
        for (const [from, to] of merged) {
          const [cbeg, ubeg, cend, uend] = [Math.floor(from / SHIFT), from % SHIFT, Math.floor(to / SHIFT), to % SHIFT];
          // The chunk ends `uend` bytes into the block at `cend`: that block is read too (its size is in its header).
          const last = uend > 0 ? blockSize(await data(cend, 18, signal), 0) : 0;
          const read = await bgzf(await data(cbeg, cend - cbeg + last, signal));
          learn(read);
          const blocks = read.map((x) => x.data);
          const bytes = joined(blocks);
          const lines = text(
            bytes.subarray(ubeg, last > 0 ? bytes.length - blocks[blocks.length - 1].length + uend : bytes.length),
          ).split('\n');
          for (const line of lines) {
            if (!line || line[0] === '#') continue;
            const f = line.split('\t');
            const pos = Number(f[1]) - 1;
            if (f[0] !== system) continue;
            const end = Number(/(?:^|;)END=(\d+)/.exec(f[7])?.[1] ?? pos + f[3].length);
            if (end <= start || pos >= stop) continue;
            for (const alt of f[4].split(',')) {
              const type =
                alt[0] === '<'
                  ? 4
                  : alt.length === f[3].length
                    ? alt.length === 1
                      ? 0
                      : 1
                    : alt.length > f[3].length
                      ? 2
                      : 3;
              const row = {
                start: pos,
                stop: end,
                id: `${pos + 1}:${f[3]}>${alt}`,
                label: f[2] === '.' ? `${f[3]}>${alt}` : f[2],
                ref: f[3],
                alt,
                category: TYPES[type],
                name: f[2] === '.' ? '' : f[2],
              };
              for (const k in row) t[k].push(row[k as keyof typeof row]);
              for (const k of o.info ?? []) t[k].push(new RegExp(`(?:^|;)${k}=([^;]*)`).exec(f[7])?.[1] ?? '');
              starts.push(pos);
              stops.push(end);
            }
          }
        }
        return { table: t, start: starts, stop: stops };
      },
    },
    o.maxItems,
  );
}
