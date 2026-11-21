/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { columnValues } from '@molstar/core/repr/tree';
import {
  bgzip,
  type IndexedLine,
  type IndexOptions,
  memoryRanges,
  tabixIndex,
} from '@molstar/test-harness/genome/bgzf-writer';
import type { AssetReader, RegionAnswer } from '../../data/protocol.js';

// The BED family's generated files (seeded, on core's PCG) and what the tests share over them: BED with BED12 blocks,
// names that repeat, strands and scores; bedGraph; BEDPE; each bgzipped and tabix-indexed by the shared test-only
// writer, and an asset port over files in memory.

export interface Generated {
  readonly text: string;
  readonly lines: readonly IndexedLine[];
  readonly names: readonly string[];
}

const enc = new TextEncoder();

/** Sorted rows as text lines, with their places for the index. */
function textOf(
  rows: { seq: number; begin: number; end: number; text: string }[],
  names: readonly string[],
  header: string[] = [],
): Generated {
  rows.sort((a, b) => a.seq - b.seq || a.begin - b.begin);
  let at = 0;
  const parts: string[] = [];
  for (const h of header) {
    parts.push(`${h}\n`);
    at += enc.encode(`${h}\n`).length;
  }
  const lines: IndexedLine[] = [];
  for (const r of rows) {
    const t = `${r.text}\n`;
    const n = enc.encode(t).length;
    lines.push({ seq: r.seq, begin: r.begin, end: Math.max(r.end, r.begin + 1), start: at, stop: at + n });
    parts.push(t);
    at += n;
  }
  return { text: parts.join(''), lines, names };
}

/**
 * `n` BED12 lines over `names`: a third with several blocks, names drawn from a small set (so `name@start-end` repeats
 * only for identical extents, and a few lines are written twice), strands and scores; positions up to `span`.
 */
export function generateBed12(n: number, names: readonly string[], seed = 1, span = 400_000): Generated {
  const rng = new PCG(seed);
  const rows: { seq: number; begin: number; end: number; text: string }[] = [];
  for (let i = 0; i < n; i++) {
    const seq = Math.floor(rng.float() * names.length);
    const begin = Math.floor(rng.float() * span);
    const length = 1 + Math.floor(rng.float() ** 3 * 20_000);
    const end = begin + length;
    const count = rng.float() < 0.35 ? 1 + Math.floor(rng.float() * 5) : 1;
    // blocks inside the item, sorted, the last ending at its end
    const starts: number[] = [];
    const sizes: number[] = [];
    if (count === 1 || length < 4 * count) {
      starts.push(0);
      sizes.push(length);
    } else {
      const step = Math.floor(length / count);
      for (let b = 0; b < count; b++) {
        const s = b * step;
        const z = b === count - 1 ? length - s : Math.max(1, Math.floor(step / 2));
        starts.push(s);
        sizes.push(z);
      }
    }
    const strand = ['+', '-', '.'][Math.floor(rng.float() * 3)];
    const name = `f${Math.floor(rng.float() * 40)}`;
    const score = Math.floor(rng.float() * 1000);
    const text = [
      names[seq],
      begin,
      end,
      name,
      score,
      strand,
      begin + starts[0],
      end,
      '0,0,0',
      starts.length,
      `${sizes.join(',')},`,
      `${starts.join(',')},`,
    ].join('\t');
    rows.push({ seq, begin, end, text });
    if (rng.float() < 0.02) rows.push({ seq, begin, end, text });
  }
  return textOf(rows, names, ['track name=generated', '# a comment']);
}

/** `n` bedGraph intervals over `names`, some overlapping. */
export function generateBedGraph(n: number, names: readonly string[], seed = 2, span = 200_000): Generated {
  const rng = new PCG(seed);
  const rows: { seq: number; begin: number; end: number; text: string }[] = [];
  for (let i = 0; i < n; i++) {
    const seq = Math.floor(rng.float() * names.length);
    const begin = Math.floor(rng.float() * span);
    const end = begin + 1 + Math.floor(rng.float() * 300);
    const value = Math.round(rng.float() * 1000) / 10;
    rows.push({ seq, begin, end, text: `${names[seq]}\t${begin}\t${end}\t${value}` });
  }
  return textOf(rows, names, ['track type=bedGraph']);
}

/** `n` BEDPE lines (10 columns) over `names`, ends anywhere; positions up to `span`. */
export function generateBedpe(n: number, names: readonly string[], seed = 3, span = 1_000_000, named = true): string {
  const rng = new PCG(seed);
  const lines: string[] = [];
  for (let i = 0; i < n; i++) {
    const a = names[Math.floor(rng.float() * names.length)];
    const b = names[Math.floor(rng.float() * names.length)];
    const s1 = Math.floor(rng.float() * span);
    const s2 = Math.floor(rng.float() * span);
    const e1 = s1 + 1 + Math.floor(rng.float() * 2000);
    const e2 = s2 + 1 + Math.floor(rng.float() * 2000);
    const name = named ? (rng.float() < 0.2 ? '.' : `link${i}`) : '.';
    const strands = `${rng.float() < 0.5 ? '+' : '-'}\t${rng.float() < 0.5 ? '+' : '-'}`;
    lines.push(`${a}\t${s1}\t${e1}\t${b}\t${s2}\t${e2}\t${name}\t${Math.floor(rng.float() * 100)}\t${strands}`);
  }
  return `${lines.join('\n')}\n`;
}

/** BED's index options (0-based, columns 1, 2, 3, `#` meta) over `names`. */
export const bedIndex = (names: readonly string[]): IndexOptions => ({
  kind: 'tbi',
  format: 0,
  zeroBased: true,
  columns: [1, 2, 3],
  meta: '#',
  skip: 0,
  names,
});

/** A generated file bgzipped (blocks of `blockSize`) with its tabix index. */
export function indexed(g: Generated, blockSize = 4000) {
  const file = bgzip(enc.encode(g.text), blockSize);
  return { file: file.bytes, index: tabixIndex(g.lines, file, bedIndex(g.names)) };
}

/** An asset port over files in memory, with range reads unless `ranges` is false. */
export function memoryAssets(files: Record<string, Uint8Array | string>, ranges = true): AssetReader {
  const bytes = (uri: string) => {
    const f = files[uri];
    if (f === undefined) throw new Error(`no file ${uri}`);
    return typeof f === 'string' ? enc.encode(f) : f;
  };
  return {
    readText: async (uri) => new TextDecoder().decode(bytes(uri)),
    readBytes: async (uri) => bytes(uri),
    ...(ranges ? { openRange: async (uri: string) => memoryRanges(bytes(uri)) } : {}),
  };
}

/** An answer's ranges list as plain rows (its columns' values by row), for comparing two reads. */
export function plainRows(a: RegionAnswer): unknown[] {
  if (!a.ok) return [a.kind];
  const t = a.tables.ranges ?? a.tables.values;
  if (!t) return [];
  if (t.form === 'positions') return [t.start, ...Array.from(t.value, (v) => (Number.isNaN(v) ? null : v))];
  if (t.form !== 'list') return [t.form];
  const cols: [string, unknown[]][] = [
    ['start', Array.from(t.start)],
    ['end', Array.from(t.end)],
    ['strand', t.strand ? Array.from(t.strand) : []],
    ['id', t.id ? columnValues(t.id) : []],
    ['itemStart', t.itemStart ? Array.from(t.itemStart) : []],
    ['itemEnd', t.itemEnd ? Array.from(t.itemEnd) : []],
    ...(t.fields ?? []).map(
      (f) => [f.name, columnValues(f.column).map((v) => (Number.isNaN(v) ? null : v))] as [string, unknown[]],
    ),
  ];
  return Array.from(t.start, (_, i) => Object.fromEntries(cols.map(([k, v]) => [k, v[i]])));
}
