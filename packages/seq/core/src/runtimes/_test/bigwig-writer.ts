/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';

// A bigWig writer for tests (never shipped): chromosomes of per-position values (NaN: no data) written as UCSC's
// format lays them out (the header, zoom headers, a total summary, the chromosome B+ tree, the data sections and
// their R-tree, each zoom level's records and R-tree), so the runtime reads what a real file holds. Sections of a
// chromosome are bedGraph, varStep or fixedStep; the zoom records are aligned to their reduction (`bins` of the
// runtime handles records that are not); blocks are zlib-compressed or not; little- or big-endian. Small trees
// (`blockSize`) give an index with inner nodes on a small file.

export type SectionType = 'bedGraph' | 'varStep' | 'fixedStep';

export interface WriterChrom {
  readonly name: string;
  readonly values: Float32Array;
  /** The sections' type (bedGraph by default). */
  readonly type?: SectionType;
}

export interface WriterOptions {
  /** The zoom levels' reductions (bases a record), ascending. */
  readonly zooms?: readonly number[];
  /** Items a data section (a block) holds at most: 64 by default. */
  readonly itemsPerSection?: number;
  /** Items an index node holds at most: 4 by default, so a small file has inner nodes. */
  readonly blockSize?: number;
  readonly compress?: boolean;
  readonly bigEndian?: boolean;
}

/** A signal of `length` positions from `seed`: values in [0, 1) by 1/256, runs of a value, gaps of no data. */
export function syntheticSignal(seed: number, length: number, gaps: readonly (readonly [number, number])[] = []) {
  const rng = new PCG(seed);
  const out = new Float32Array(length);
  for (let i = 0; i < length; ) {
    const run = 1 + (rng.int() % 6);
    const v = (rng.int() % 256) / 256;
    for (let j = 0; j < run && i < length; j++, i++) out[i] = v;
  }
  for (const [a, b] of gaps) out.fill(Number.NaN, a, b);
  return out;
}

async function deflate(bytes: Uint8Array): Promise<Uint8Array> {
  const stream = new Blob([bytes as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new CompressionStream('deflate'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** A growing little- or big-endian buffer. */
class Out {
  private buf = new Uint8Array(1024);
  private dv = new DataView(this.buf.buffer);
  length = 0;
  constructor(readonly le: boolean) {}
  private room(n: number) {
    if (this.length + n <= this.buf.length) return;
    const next = new Uint8Array(Math.max(this.buf.length * 2, this.length + n));
    next.set(this.buf);
    this.buf = next;
    this.dv = new DataView(next.buffer);
  }
  u8(x: number) {
    this.room(1);
    this.dv.setUint8(this.length, x);
    this.length += 1;
  }
  u16(x: number) {
    this.room(2);
    this.dv.setUint16(this.length, x, this.le);
    this.length += 2;
  }
  u32(x: number) {
    this.room(4);
    this.dv.setUint32(this.length, x, this.le);
    this.length += 4;
  }
  u64(x: number) {
    const lo = x % 0x100000000;
    const hi = Math.floor(x / 0x100000000);
    this.u32(this.le ? lo : hi);
    this.u32(this.le ? hi : lo);
  }
  f32(x: number) {
    this.room(4);
    this.dv.setFloat32(this.length, x, this.le);
    this.length += 4;
  }
  f64(x: number) {
    this.room(8);
    this.dv.setFloat64(this.length, x, this.le);
    this.length += 8;
  }
  bytes(b: Uint8Array) {
    this.room(b.length);
    this.buf.set(b, this.length);
    this.length += b.length;
  }
  /** Writes `n` zero bytes. */
  zeros(n: number) {
    this.room(n);
    this.buf.fill(0, this.length, this.length + n);
    this.length += n;
  }
  done() {
    return this.buf.slice(0, this.length);
  }
}

/** An R-tree leaf item: a block's extent and where it is. */
interface Item {
  chrom: number;
  start: number;
  end: number;
  offset: number;
  size: number;
}

/** Writes an R-tree over `items` (file order) at the end of `out`: header, then its levels from the root down. */
function writeRTree(out: Out, items: readonly Item[], blockSize: number, itemsPerSlot: number, endOfData: number) {
  // levels from the leaves up: each node an extent and its children
  type Node = { chrom: number; start: number; endChrom: number; end: number; children: (Node | Item)[] };
  const extent = (xs: (Node | Item)[]): Node => ({
    chrom: xs[0].chrom,
    start: xs[0].start,
    endChrom: 'endChrom' in xs[xs.length - 1] ? (xs[xs.length - 1] as Node).endChrom : xs[xs.length - 1].chrom,
    end: xs[xs.length - 1].end,
    children: xs,
  });
  let level: Node[] = [];
  for (let i = 0; i < items.length; i += blockSize) level.push(extent(items.slice(i, i + blockSize)));
  if (!level.length) level = [{ chrom: 0, start: 0, endChrom: 0, end: 0, children: [] }];
  const levels = [level];
  while (level.length > 1) {
    const up: Node[] = [];
    for (let i = 0; i < level.length; i += blockSize) up.push(extent(level.slice(i, i + blockSize)));
    levels.unshift(up);
    level = up;
  }
  const root = levels[0][0];
  out.u32(0x2468ace0);
  out.u32(blockSize);
  out.u64(items.length);
  out.u32(root.chrom);
  out.u32(root.start);
  out.u32(root.endChrom);
  out.u32(root.end);
  out.u64(endOfData);
  out.u32(itemsPerSlot);
  out.u32(0);
  // offsets: the levels one after another, each node its items
  const leafDepth = levels.length - 1;
  const offsets = new Map<Node, number>();
  let at = out.length;
  levels.forEach((nodes, d) => {
    for (const n of nodes) {
      offsets.set(n, at);
      at += 4 + n.children.length * (d === leafDepth ? 32 : 24);
    }
  });
  levels.forEach((nodes, d) => {
    for (const n of nodes) {
      const leaf = d === leafDepth;
      out.u8(leaf ? 1 : 0);
      out.u8(0);
      out.u16(n.children.length);
      for (const c of n.children) {
        const endChrom = 'endChrom' in c ? c.endChrom : c.chrom;
        out.u32(c.chrom);
        out.u32(c.start);
        out.u32(endChrom);
        out.u32(c.end);
        if (leaf) {
          out.u64((c as Item).offset);
          out.u64((c as Item).size);
        } else out.u64(offsets.get(c as Node)!);
      }
    }
  });
}

/** The data sections of a chromosome: runs of data cut into sections of at most `per` items. */
function sections(id: number, c: WriterChrom, per: number, le: boolean): Uint8Array[] {
  const v = c.values;
  const type = c.type ?? 'bedGraph';
  // items: [start, end, value], a bedGraph item a run of one value, the others one position each
  const items: [number, number, number][] = [];
  for (let i = 0; i < v.length; ) {
    if (Number.isNaN(v[i])) {
      i++;
      continue;
    }
    let j = i + 1;
    if (type === 'bedGraph') while (j < v.length && v[j] === v[i]) j++;
    items.push([i, j, v[i]]);
    i = j;
  }
  const out: Uint8Array[] = [];
  for (let i = 0; i < items.length; ) {
    // a fixedStep section is a run without gaps
    let j = i + 1;
    while (j < items.length && j - i < per && (type !== 'fixedStep' || items[j][0] === items[j - 1][1])) j++;
    const part = items.slice(i, j);
    const s = new Out(le);
    s.u32(id);
    s.u32(part[0][0]);
    s.u32(part[part.length - 1][1]);
    s.u32(type === 'fixedStep' ? 1 : 0);
    s.u32(type === 'bedGraph' ? 0 : 1);
    s.u8(type === 'bedGraph' ? 1 : type === 'varStep' ? 2 : 3);
    s.u8(0);
    s.u16(part.length);
    for (const [a, b, x] of part) {
      if (type === 'bedGraph') {
        s.u32(a);
        s.u32(b);
      } else if (type === 'varStep') s.u32(a);
      s.f32(x);
    }
    out.push(s.done());
    i = j;
  }
  return out;
}

/** The zoom records of a chromosome at `reduction`: one a bin with data, aligned to multiples of it. */
function zoomRecords(id: number, values: Float32Array, reduction: number) {
  const out: { start: number; end: number; count: number; min: number; max: number; sum: number; sq: number }[] = [];
  for (let a = 0; a < values.length; a += reduction) {
    const b = Math.min(values.length, a + reduction);
    let count = 0;
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let sum = 0;
    let sq = 0;
    for (let p = a; p < b; p++) {
      const x = values[p];
      if (Number.isNaN(x)) continue;
      count++;
      min = Math.min(min, x);
      max = Math.max(max, x);
      sum += x;
      sq += x * x;
    }
    if (count) out.push({ start: a, end: b, count, min, max, sum, sq });
  }
  return out;
}

/** The bytes of a bigWig of `chroms` (in the order given: their ids). */
export async function writeBigWig(chroms: readonly WriterChrom[], options: WriterOptions = {}): Promise<Uint8Array> {
  const le = !options.bigEndian;
  const zooms = options.zooms ?? [];
  const per = options.itemsPerSection ?? 64;
  const blockSize = Math.max(2, options.blockSize ?? 4);
  const compress = options.compress ?? true;
  const out = new Out(le);
  let uncompress = 0;
  const block = async (raw: Uint8Array) => {
    uncompress = Math.max(uncompress, raw.length);
    return compress ? deflate(raw) : raw;
  };

  // header, patched at the end
  out.zeros(64 + 24 * zooms.length);
  const totalSummary = out.length;
  let covered = 0;
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;
  let sum = 0;
  let sq = 0;
  for (const c of chroms)
    for (const x of c.values) {
      if (Number.isNaN(x)) continue;
      covered++;
      min = Math.min(min, x);
      max = Math.max(max, x);
      sum += x;
      sq += x * x;
    }
  out.u64(covered);
  out.f64(covered ? min : 0);
  out.f64(covered ? max : 0);
  out.f64(sum);
  out.f64(sq);

  // the chromosome tree: one leaf of every chromosome, names zero-padded to the longest
  const chromTree = out.length;
  const keySize = Math.max(1, ...chroms.map((c) => c.name.length));
  out.u32(0x78ca8c91);
  out.u32(Math.max(1, chroms.length));
  out.u32(keySize);
  out.u32(8);
  out.u64(chroms.length);
  out.u64(0);
  out.u8(1);
  out.u8(0);
  out.u16(chroms.length);
  for (let i = 0; i < chroms.length; i++) {
    const key = new Uint8Array(keySize);
    for (let k = 0; k < chroms[i].name.length; k++) key[k] = chroms[i].name.charCodeAt(k);
    out.bytes(key);
    out.u32(i);
    out.u32(chroms[i].values.length);
  }

  // the full data: a count of sections, then the sections, then their index
  const fullData = out.length;
  const items: Item[] = [];
  const all = chroms.flatMap((c, i) => sections(i, c, per, le).map((s) => [i, s] as const));
  out.u64(all.length);
  for (const [id, s] of all) {
    const v = new DataView(s.buffer, s.byteOffset, s.byteLength);
    const bytes = await block(s);
    items.push({
      chrom: id,
      start: v.getUint32(4, le),
      end: v.getUint32(8, le),
      offset: out.length,
      size: bytes.length,
    });
    out.bytes(bytes);
  }
  const fullIndex = out.length;
  writeRTree(out, items, blockSize, per, fullIndex);

  // each zoom level: a count of records, the records in blocks, their index
  const zoomHeads: { reduction: number; data: number; index: number }[] = [];
  for (const reduction of zooms) {
    const data = out.length;
    const records = chroms.flatMap((c, i) => zoomRecords(i, c.values, reduction).map((r) => ({ id: i, ...r })));
    out.u32(records.length);
    const zoomItems: Item[] = [];
    for (let i = 0; i < records.length; ) {
      let j = i + 1;
      while (j < records.length && j - i < per && records[j].id === records[i].id) j++;
      const s = new Out(le);
      for (const r of records.slice(i, j)) {
        s.u32(r.id);
        s.u32(r.start);
        s.u32(r.end);
        s.u32(r.count);
        s.f32(r.min);
        s.f32(r.max);
        s.f32(r.sum);
        s.f32(r.sq);
      }
      const bytes = await block(s.done());
      zoomItems.push({
        chrom: records[i].id,
        start: records[i].start,
        end: records[j - 1].end,
        offset: out.length,
        size: bytes.length,
      });
      out.bytes(bytes);
      i = j;
    }
    const index = out.length;
    writeRTree(out, zoomItems, blockSize, per, index);
    zoomHeads.push({ reduction, data, index });
  }

  // the header
  const h = new Out(le);
  h.u32(0x888ffc26);
  h.u16(4);
  h.u16(zooms.length);
  h.u64(chromTree);
  h.u64(fullData);
  h.u64(fullIndex);
  h.u16(0);
  h.u16(0);
  h.u64(0);
  h.u64(totalSummary);
  h.u32(compress ? uncompress : 0);
  h.u64(0);
  for (const z of zoomHeads) {
    h.u32(z.reduction);
    h.u32(0);
    h.u64(z.data);
    h.u64(z.index);
  }
  const bytes = out.done();
  bytes.set(h.done(), 0);
  return bytes;
}
