/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The container BigWig and BigBed share (README, "BigWig and BigBed"): the header, the chromosome tree, the R-tree that
// finds the data blocks of a region, and the blocks' inflate. What a block holds is the format's (bigwig.ts, bigbed.ts).

import { cstring, inflate, type RangeReader, u64, view } from './ranges.js';

const MAGIC = { bigwig: 0x888ffc26, bigbed: 0x8789f2eb } as const;
const CHROM_TREE = 0x78ca8c91;
const HEADER = 64;
const ZOOM_HEADER = 24;
const MAX_ZOOMS = 16;
/** How deep the chromosome tree and the R-tree may be: a genome's are a few levels (GI-06). */
const MAX_TREE_DEPTH = 32;
/** Header-derived read lengths (GI-11): a million contigs' tree is about 40 MB, an autoSql a few KB. */
const MAX_CHROM_TREE = 64 << 20;
const MAX_AUTOSQL = 1 << 20;
const MAX_NODE_ITEMS = 65536;

/** A length the header gives, checked before it is read (GI-11: a damaged one asked for a terabyte). */
function bounded(what: string, length: number, max: number): number {
  if (!(length >= 0 && length <= max))
    throw new Error(`the header's ${what} is ${length} bytes: past its bound of ${max}`);
  return length;
}

export interface Chrom {
  readonly name: string;
  readonly id: number;
  readonly size: number;
}
export interface Zoom {
  /** Positions per summary of this level, as the file says (`levels`, in the manifest). */
  readonly reduction: number;
  readonly indexOffset: number;
}
export interface Bbi {
  readonly kind: keyof typeof MAGIC;
  readonly le: boolean;
  readonly chroms: Map<string, Chrom>;
  /** Finest first. */
  readonly zooms: Zoom[];
  readonly fieldCount: number;
  readonly definedFieldCount: number;
  readonly autoSql: string | null;
  /** The whole file's value range (a BigWig's `domain`); null: no summary. */
  readonly range: [number, number] | null;
  /** The decoded blocks of the index at `indexOffset` that may hold [start, stop) of `chrom`: raw bytes, inflated. */
  blocks(indexOffset: number, chrom: Chrom, start: number, stop: number, signal?: AbortSignal): Promise<Uint8Array[]>;
  readonly dataIndexOffset: number;
}

export async function openBbi(read: RangeReader, signal?: AbortSignal): Promise<Bbi> {
  const head = await read(0, HEADER + MAX_ZOOMS * ZOOM_HEADER, signal);
  const dv = view(head);
  const kind = (Object.keys(MAGIC) as (keyof typeof MAGIC)[]).find(
    (k) => dv.getUint32(0, true) === MAGIC[k] || dv.getUint32(0, false) === MAGIC[k],
  );
  if (!kind) throw new Error('not a BigWig or BigBed file: no magic number');
  const le = dv.getUint32(0, true) === MAGIC[kind]; // the file's own byte order, whatever the machine that wrote it
  const [zoomCount, treeAt, dataAt, indexAt] = [dv.getUint16(6, le), u64(dv, 8, le), u64(dv, 16, le), u64(dv, 24, le)];
  const [fieldCount, definedFieldCount, sqlAt, summaryAt] = [
    dv.getUint16(32, le),
    dv.getUint16(34, le),
    u64(dv, 36, le),
    u64(dv, 44, le),
  ];
  const uncompressed = dv.getUint32(52, le);
  if (zoomCount > MAX_ZOOMS) throw new Error(`${zoomCount} zoom levels: at most ${MAX_ZOOMS} are read`);
  const zooms = Array.from({ length: zoomCount }, (_, k) => ({
    reduction: dv.getUint32(HEADER + k * ZOOM_HEADER, le),
    indexOffset: u64(dv, HEADER + k * ZOOM_HEADER + 16, le),
  })).sort((a, b) => a.reduction - b.reduction);

  // The chromosome tree is one run of nodes between its offset and the data's: read it whole and walk it in memory.
  const tree = await read(treeAt, bounded('chromosome tree', dataAt - treeAt, MAX_CHROM_TREE), signal);
  const tv = view(tree);
  if (tv.getUint32(0, le) !== CHROM_TREE) throw new Error('bad chromosome tree magic');
  const [keySize, valSize, items] = [tv.getUint32(8, le), tv.getUint32(12, le), u64(tv, 16, le)];
  const chroms = new Map<string, Chrom>();
  // Iterative, each node once, at most MAX_TREE_DEPTH deep and `items` leaves (GI-06: a node pointing twice at the
  // next one doubled the walk per level, 2^30 leaves from a 1 KB file).
  const seen = new Set<number>();
  const stack: [number, number][] = [[32, 0]];
  while (stack.length > 0) {
    const [at, depth] = stack.pop() as [number, number];
    if (seen.has(at)) continue;
    seen.add(at);
    if (depth > MAX_TREE_DEPTH) throw new Error(`the chromosome tree is deeper than ${MAX_TREE_DEPTH}`);
    const [leaf, n] = [tree[at], tv.getUint16(at + 2, le)];
    const step = keySize + (leaf ? valSize : 8);
    // A node's children are pushed last first, so they pop in the tree's (sorted) order: the manifest's `systems`.
    if (!leaf)
      for (let i = n - 1; i >= 0; i--) stack.push([u64(tv, at + 4 + i * step + keySize, le) - treeAt, depth + 1]);
    else
      for (let i = 0, e = at + 4; i < n; i++, e += step) {
        const name = cstring(tree.subarray(e, e + keySize));
        chroms.set(name, { name, id: tv.getUint32(e + keySize, le), size: tv.getUint32(e + keySize + 4, le) });
        if (chroms.size > items) throw new Error(`the chromosome tree holds more than its ${items} items`);
      }
  }

  const sqlLength = sqlAt > 0 ? bounded('autoSql', (summaryAt || treeAt) - sqlAt, MAX_AUTOSQL) : 0;
  const sql = sqlAt > 0 ? cstring(await read(sqlAt, sqlLength, signal)) : null;
  let range: [number, number] | null = null;
  if (summaryAt > 0) {
    const s = view(await read(summaryAt, 40, signal));
    if (s.getBigUint64(0, le) > 0n) range = [s.getFloat64(8, le), s.getFloat64(16, le)];
  }

  /** The leaf items of the R-tree at `at` whose chromosome range meets (id, start, stop): [offset, size] of a block. */
  async function search(
    at: number,
    id: number,
    start: number,
    stop: number,
    signal?: AbortSignal,
  ): Promise<[number, number][]> {
    const root = view(await read(at, 48, signal));
    const [blockSize, rootAt] = [root.getUint32(4, le), at + 48];
    if (blockSize > MAX_NODE_ITEMS)
      throw new Error(`the R-tree's block size ${blockSize}: past ${MAX_NODE_ITEMS} items`);
    const meets = (v: DataView, o: number) => {
      const [sc, sb, ec, eb] = [
        v.getUint32(o, le),
        v.getUint32(o + 4, le),
        v.getUint32(o + 8, le),
        v.getUint32(o + 12, le),
      ];
      return (sc < id || (sc === id && sb < stop)) && (ec > id || (ec === id && eb > start));
    };
    // Each node once and at most MAX_TREE_DEPTH deep (GI-06: a root pointing at itself twice doubled the visits).
    const seen = new Set<number>();
    const visit = async (off: number, depth = 0): Promise<[number, number][]> => {
      if (seen.has(off)) return [];
      seen.add(off);
      if (depth > MAX_TREE_DEPTH) throw new Error(`the R-tree is deeper than ${MAX_TREE_DEPTH}`);
      const node = await read(off, 4 + blockSize * 32, signal);
      const nv = view(node);
      const [leaf, n] = [node[0], nv.getUint16(2, le)];
      const found: Promise<[number, number][]>[] = [];
      for (let i = 0, o = 4; i < n; i++, o += leaf ? 32 : 24) {
        if (!meets(nv, o)) continue;
        if (leaf) found.push(Promise.resolve([[u64(nv, o + 16, le), u64(nv, o + 24, le)]]));
        else found.push(visit(u64(nv, o + 16, le), depth + 1));
      }
      return (await Promise.all(found)).flat();
    };
    return visit(rootAt);
  }
  return {
    kind,
    le,
    chroms,
    zooms,
    fieldCount,
    definedFieldCount,
    autoSql: sql,
    range,
    dataIndexOffset: indexAt,
    async blocks(indexOffset, chrom, start, stop, signal) {
      const found = await search(indexOffset, chrom.id, start, stop, signal);
      return Promise.all(
        found.map(async ([offset, size]) => {
          const raw = await read(offset, size, signal);
          return uncompressed > 0 ? inflate(raw, 'deflate', uncompressed) : raw;
        }),
      );
    },
  };
}
