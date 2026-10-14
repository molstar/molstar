/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The projection core's fixtures (VW5): views opened with a projection over them and readers of what it lights, case
// 14's gene on either strand (filtered, or placed on a chromosome), constructs (di-ubiquitin, a fusion) and rows bound
// to chains. Shared by project.test.ts and project.heap.ts.

import { createSeqViewSpec } from '@molstar/seq-core/builder';
import { type CoordinateSystem, Extents, Mappings, type SeqSet } from '@molstar/seq-core/model';
import { expect } from 'vitest';
import {
  compile,
  type Projection,
  projectionOf,
  type ResolvedView,
  resolveView,
  type SeqViewSpecDoc,
  supportedAll,
} from '../../src/format/index.js';
import { fixtureLoaders, registries } from './load-fixtures.js';

const reg = registries();
const all = supportedAll(reg);
export type Entry = [key: string, start: number, stop: number];

export async function open(doc: SeqViewSpecDoc, ref = 'v') {
  const { set, issues } = await compile(doc, reg, fixtureLoaders());
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const r = resolveView(doc, set, ref, reg, all);
  expect(r.issues.filter((i) => i.severity === 'error')).toEqual([]);
  return opened(set, r.view);
}
export function opened(set: SeqSet, view: ResolvedView) {
  const p = projectionOf(view, set);
  const system = (key: string): CoordinateSystem => {
    const s = { kind: -1, id: -1, row: -1 };
    if (!Mappings.systemOf(set, key, s)) throw new Error(`no system ${key}`);
    return s;
  };
  const extents = (entries: Entry[]) => {
    const e = Extents.create(1);
    for (const [key, s, t] of entries) Extents.add(e, system(key), s, t);
    return e;
  };
  /** A block by its label. */
  const block = (label: string) => {
    for (let b = 0; b < view.blocks.count; b++) if (set.strings.get(view.blocks.label[b]) === label) return b;
    throw new Error(`no block ${label}`);
  };
  /** The runs the entries light on a block (closed first when `closure` is containing), merged and sorted. */
  const lit = (label: string, entries: Entry[] | Extents, closure: 'element' | 'containing' = 'element') => {
    const e = Array.isArray(entries) ? extents(entries) : entries;
    if (closure === 'containing') p.close(e);
    return runs(p, e, block(label));
  };
  /** The extents' entries by key. */
  const keyed = (e: Extents) =>
    Array.from({ length: e.count }, (_, i) => {
      const k = e.kind[i];
      const key =
        k === 0 ? set.molecules[e.id[i]]?.key : k === 1 ? set.units.key[e.id[i]] : (set.alignments[e.id[i]]?.key ?? -1);
      const name = `${set.strings.get(key ?? -1)}${k === 2 ? `/${e.row[i]}` : k === 3 ? ' columns' : ''}`;
      return [name, e.start[i], e.stop[i]];
    });
  /** The entries after closure, sorted. */
  const closed = (entries: Entry[]) => {
    const e = extents(entries);
    p.close(e);
    return keyed(e).sort((a, b) => String(a).localeCompare(String(b)));
  };
  return { set, view, p, system, extents, block, lit, keyed, closed };
}

/** Merged sorted runs of a projection. */
export function runs(p: Projection, e: Extents, block: number): number[][] {
  const out: number[][] = [];
  p.project(e, block, { add: (s, t) => out.push([s, t]) });
  out.sort((a, b) => a[0] - b[0]);
  const merged: number[][] = [];
  for (const r of out) {
    const last = merged.at(-1);
    if (last && r[0] <= last[1]) last[1] = Math.max(last[1], r[1]);
    else merged.push([...r]);
  }
  return merged;
}

export const seq = (n: number, offset = 0) =>
  Array.from({ length: n }, (_, i) => 'ACDEFGHIKLMNPQRSTVWY'[(i + offset) % 20]).join('');
const GENE =
  'GAATTCGCCACC' +
  'ATGGCTTCCAAGGAACTGACCCCACGTGACGGCC' +
  'GTAAGTCTTTCCTAACTTCTGTTTCCACAG' +
  'ATTGGTACTTCAACAGCCAGGTGGAGTAA' +
  'GCGGCCGCATAAAGT';
const PROTEIN = 'MASKELTPRDGHWYFNSQVE';
const COMPLEMENT: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A' };
const reverseComplement = (s: string) =>
  Array.from(s, (c) => COMPLEMENT[c])
    .reverse()
    .join('');

/**
 * Case 14's gene (or its reverse complement on strand -) with the protein's letters over it, a protein section and,
 * optionally, a filter to the exons and a chromosome placing the gene at offset 100 in a section of its own.
 */
export function gene(
  o: { minus?: boolean; filter?: boolean; chromosome?: boolean; hideTrack?: boolean; noProteinSection?: boolean } = {},
) {
  const b = createSeqViewSpec();
  const g = b.data.sequence({ ref: 'gene', alphabet: 'dna', value: o.minus ? reverseComplement(GENE) : GENE });
  b.data.sequence({ ref: 'protein', alphabet: 'protein', value: PROTEIN });
  const cds = o.minus
    ? {
        cds: [
          [74, 108],
          [15, 44],
        ],
        strand: '-',
      }
    : {
        cds: [
          [12, 46],
          [76, 105],
        ],
      };
  b.data.mapping({ ref: 'm', from: 'protein', to: 'gene', method: 'codon', ...cds } as never);
  const exons = b.data.annotation({
    ref: 'exons',
    geometry: 'ranges',
    target: g,
    value: { start: [0, 76], stop: [46, 120] },
  });
  const v = b.view({ ref: 'v', label: 'V', ...(o.filter ? { positions: exons } : {}) });
  v.block({ sequence: g, label: 'gene' })
    .track({ label: 'Protein', hidden: o.hideTrack ?? false })
    .letters({ of: 'protein' });
  if (!o.noProteinSection)
    v.section({ ref: 'v.protein' }).block({ sequence: 'protein', label: 'protein' }).track({ label: 'P' }).letters();
  if (o.chromosome) {
    b.data.sequence({ ref: 'chr', alphabet: 'dna', value: `${'A'.repeat(100)}${GENE}${'A'.repeat(80)}` });
    b.data.mapping({ ref: 'placed', from: 'gene', to: 'chr', method: 'segments', segments: [[0, 100, 120]] });
    v.section({ ref: 'v.chr' }).block({ sequence: 'chr', label: 'chr' }).track({ label: 'DNA' }).letters();
  }
  return b.getDocument();
}

/** Ub, linear di-ubiquitin (one mapping, two segments over Ub) and a GFP-p53 fusion with a linker. */
export function constructs() {
  const b = createSeqViewSpec();
  for (const [ref, n] of [
    ['ub', 76],
    ['diub', 152],
    ['gfp', 238],
    ['p53', 393],
    ['fusion', 636],
  ] as const)
    b.data.sequence({ ref, value: seq(n, ref.length) });
  const ub = [
    [0, 0, 76],
    [76, 0, 76],
  ];
  b.data.mapping({ ref: 'diub-ub', from: 'diub', to: 'ub', method: 'segments', segments: ub } as never);
  b.data.mapping({ ref: 'f-gfp', from: 'fusion', to: 'gfp', method: 'segments', segments: [[0, 0, 238]] } as never);
  b.data.mapping({ ref: 'f-p53', from: 'fusion', to: 'p53', method: 'segments', segments: [[243, 0, 393]] } as never);
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: 'ub', label: 'ub' }).track({ label: 'L' }).letters();
  for (const s of ['diub', 'gfp', 'p53', 'fusion'])
    v.section({ ref: `v.${s}` })
      .block({ sequence: s, label: s })
      .track({ label: 'L' })
      .letters();
  return b.getDocument();
}

export function rows(n = 4) {
  const b = createSeqViewSpec();
  const value = ['MK-VL', 'MKIVL', 'MKIV-', ...Array.from({ length: n - 3 }, () => 'MKIVL')];
  const a = b.data.alignment({ ref: 'a', value });
  for (const ref of ['x', 'y', 'z']) b.data.sequence({ ref, value: 'GGMKVLGG' });
  b.data.mapping({ ref: 'r0', from: 'a/0', to: 'x', method: 'segments', segments: [[0, 2, 4]] } as never);
  for (const to of ['y', 'z'])
    b.data.mapping({ ref: `r1${to}`, from: 'a/1', to, method: 'segments', segments: [[0, 2, 2]] } as never);
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ alignment: a, label: 'a' }).track({ label: 'L' }).letters();
  const w = v.section({ ref: 'v.chains' });
  w.block({ sequence: 'x', label: 'x' });
  w.block({ sequence: 'y', label: 'y' });
  return b.getDocument();
}
