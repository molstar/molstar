/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

// What the rows and pairs layers name (_spec/README.md "Render", Accessibility; `LayerType.items`): a row by its label,
// a pair by its arcs and by its dots, on a track, a ring and a grid, as the layer returns them and as the exported
// SVG's accessible layer titles them.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type Frame, type ItemDescription, PairsType, renderSvg, RowsType } from '../index.js';
import { settle } from './helpers.js';
import { arcsSpec, gridSpec } from './pairs-fixtures.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLS--', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1)gc:1)rest:1)root;';
const Pairs: [number, number][] = [
  [2, 30],
  [4, 20],
  [6, 12],
];
const PairKeys = ['2-30', '4-20', '6-12'];

type Spec = ReturnType<typeof rowsSpec>;

function rowsSpec(params: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln, tree, ...params } as never)),
    ),
  );
  return b.build();
}

type Described = { name: string; items?: (lay: never) => readonly ItemDescription[] };

/** The items every layer of a frame of `type` names (the track's layers, then the section's own). */
function itemsOf(f: Frame, type: Described) {
  const out: ItemDescription[] = [];
  for (const s of f.layout.sections)
    for (const l of [...s.lays.flatMap((x) => x ?? []), ...s.sectionLays])
      if (l.type.name === type.name) out.push(...type.items!(l.lay as never));
  return out;
}

/** The titles of the exported SVG's item boxes. */
async function titles(spec: Spec, size = { width: 600, height: 300 }) {
  const r = await renderSvg(spec, { size });
  const doc = new DOMParser().parseFromString(r.svg, 'image/svg+xml');
  return [...doc.querySelectorAll('[data-seq-item] > title')].map((t) => t.textContent ?? '');
}

const pairTitles = (t: string[]) =>
  t
    .filter((x) => /^\d+-\d+, /.test(x))
    .map((x) => x.split(',')[0])
    .sort();

describe('the rows layer’s items', () => {
  it('names each row in view by its label, in the order, one item a row across the band', async () => {
    const s = await settle(rowsSpec());
    const f = await s.frame({
      size: { width: 600, height: 300 },
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6 } },
    });
    const items = itemsOf(f, RowsType);
    expect(items.map((i) => i.label)).toEqual(Labels);
    expect(items.map((i) => [i.start, i.end])).toEqual(Labels.map((_, k) => [k, k + 1]));
    for (const i of items) {
      expect(i.x1).toBeGreaterThan(i.x0);
      expect(i.height).toBeGreaterThan(0);
    }
    for (let k = 1; k < items.length; k++) expect(items[k].y!).toBeGreaterThan(items[k - 1].y!);
  });

  it('names the rows in the exported SVG: each label with its extent', async () => {
    const t = await titles(rowsSpec());
    expect(t.map((x) => x.split(',')[0])).toEqual(Labels);
  });

  it('leaves a folded clade’s summary row without an item, its rows named no more', async () => {
    const s = await settle(rowsSpec({ order: 'tree' }));
    const f = await s.frame({
      size: { width: 600, height: 300 },
      settled: true,
      viewports: {
        '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { rowScale: 20, folds: ['hb'] } as never } },
      },
    });
    // HBA and HBB fold into one summary row, which `rowItems` leaves out (rows.ts)
    expect(itemsOf(f, RowsType).map((i) => i.label)).toEqual(['MYG/1-6', 'GLB/3-8', 'CYG/1-6']);
  });
});

describe('the pairs layer’s items', () => {
  it('names each arc on a track by its key and its first end', async () => {
    const s = await settle(arcsSpec({ pairs: Pairs }));
    const f = await s.frame({ size: { width: 720, height: 300 }, settled: true });
    const items = itemsOf(f, PairsType);
    expect(items.map((i) => [i.label, i.start, i.end]).sort()).toEqual([
      ['2-30', 2, 3],
      ['4-20', 4, 5],
      ['6-12', 6, 7],
    ]);
    for (const i of items) expect(i.x1).toBeGreaterThan(i.x0);
  });

  // BUG (reported, not fixed): `accessibleOf` (accessibility.ts) reads `TrackLayout.layers` only, and the arcs are a
  // `section` layer (`SectionLayout.sectionLays`), so the exported SVG names none of them. Remove `.fails` with the fix.
  it.fails('names each arc in the exported SVG, on a track and on a ring', async () => {
    for (const layout of [undefined, 'rings' as const]) {
      const t = await titles(arcsSpec({ pairs: Pairs, ...(layout ? { layout } : {}) }), {
        width: 720,
        height: 400,
      });
      expect(pairTitles(t)).toEqual(PairKeys);
    }
  });

  it('names a pair by its label when it has one, else its key', async () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'rna', letters: 'ACGUACGUACGUACGUACGU' });
    const p = b.data.pairs('p', {
      a: s,
      rows: {
        a_start: [1, 3],
        a_end: [2, 4],
        b_start: [15, 10],
        b_end: [16, 11],
        id: ['k1', 'k2'],
        label: ['stem 1', ''],
      } as never,
    });
    b.views.view('v', (v) => v.section('s', { of: s }).track('pairs').layer('arcs', { type: 'pairs', data: p }));
    const st = await settle(b.build());
    const f = await st.frame({ size: { width: 720, height: 300 }, settled: true });
    expect(
      itemsOf(f, PairsType)
        .map((i) => i.label)
        .sort(),
    ).toEqual(['k2', 'stem 1']);
  });

  it('names each pair by its dot on a grid', async () => {
    const s = await settle(gridSpec({ n: 40, pairs: Pairs }));
    const f = await s.frame({ size: { width: 720, height: 900 }, settled: true });
    const items = itemsOf(f, PairsType);
    expect(items.map((i) => [i.label, i.start, i.end]).sort()).toEqual([
      ['2-30', 2, 3],
      ['4-20', 4, 5],
      ['6-12', 6, 7],
    ]);
    const t = await titles(gridSpec({ n: 40, pairs: Pairs }), { width: 720, height: 900 });
    expect(pairTitles(t)).toEqual(PairKeys);
  });
});
