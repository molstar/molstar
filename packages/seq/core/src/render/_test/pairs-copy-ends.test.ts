/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A pair's ends on copies (_spec/README.md "Render", The pairs layer, Ends on blocks), by the real example's shape
// (Globins › Hemoglobin, copies split and merged: 4HHB's two α and two β chains and four hemes as units of a complex):
// an end on a unit is on that unit's block; an end on a molecule or a sequence reaches its units (EC-48); a pair whose
// ends are on one block of one system is drawn again on each copy of that block, copies as the resolver groups them
// (EC-65: one copy group, or one system and one range); a pair between two systems with an end on split copies its
// data names no unit of is drawn on no copy, said once (`seq.render.pairs.copy`). Each arc's feet at its data, at the
// fit view, zoomed onto one copy and across a tile edge; hit as its pair, carets too; canvas and SVG alike.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { drawListsToSvg, type Frame, hitView, type PairsLay, PlacePoint } from '../index.js';
import { BigRuntime } from './fixtures.js';
import { type Settled, settle } from './helpers.js';

const Size = { width: 900, height: 300 };
const Hba = 'MVLSPADKTNVKAAWGKVGA';
const Hbb = 'MVHLTPEEKSAVTALWGKVN';
/** The links' contacts: each heme's unit to a residue of its own chain. */
const Contacts = { heme: ['ha', 'hb', 'hc', 'hd'], chain: ['A', 'B', 'C', 'D'], at: [5, 6, 7, 8] };
/** Pairs within α (both ends on one system): its copies each draw them. */
const Intra: [number, number][] = [
  [2, 10],
  [4, 12],
];
/** Contacts of the heme molecule to α's residues: the data names no unit. */
const Across = [5, 7];

interface Options {
  readonly copies?: 'each' | 'merge';
  /** The pairs layer's data: the complex's links, the pairs within α, the heme's contacts on the molecules. */
  readonly data: 'links' | 'intra' | 'across';
  /** The layer's data filter (a query after the address). */
  readonly filter?: string;
  /** Blocks of α as a sequence, no complex: its whole twice, or its whole, a zoom on 5–15 and its whole again. */
  readonly blocks?: 'twice' | 'zoom';
}

function spec(o: Options): Node {
  const b = createSeqSpec({ title: 'hemoglobin' });
  const hba = b.data.sequence('hba', { alphabet: 'protein', letters: Hba });
  const hbb = b.data.sequence('hbb', { alphabet: 'protein', letters: Hbb });
  const hem = b.data.molecule('hem', { class: 'ligand', residues: ['HEM'] });
  b.data.pairs('intra', {
    a: hba,
    b: hba,
    rows: {
      a_start: Intra.map((p) => p[0]),
      a_end: Intra.map((p) => p[0] + 1),
      b_start: Intra.map((p) => p[1]),
      b_end: Intra.map((p) => p[1] + 1),
      id: Intra.map((p) => `${p[0]}-${p[1]}`),
      type: ['salt', 'hbond'],
    },
  });
  b.data.pairs('across', {
    a: hem,
    b: hba,
    rows: {
      a_start: [0, 0],
      a_end: [1, 1],
      b_start: Across,
      b_end: Across.map((r) => r + 1),
      id: ['h5', 'h7'],
      type: ['near', 'far'],
    },
  });
  b.data.complex('hb', { name: 'Hemoglobin' }, (c) => {
    for (const id of ['A', 'B', 'C', 'D']) c.unit(id, { of: id === 'A' || id === 'C' ? hba : hbb, chain: id });
    for (const id of Contacts.heme) c.unit(id, { of: hem, chain: id.toUpperCase(), label: 'heme' });
    c.links('contacts', {
      rows: {
        a_unit: Contacts.heme,
        a_start: [0, 0, 0, 0],
        a_end: [1, 1, 1, 1],
        b_unit: Contacts.chain,
        b_start: Contacts.at,
        b_end: Contacts.at.map((r) => r + 1),
        type: ['contact', 'contact', 'contact', 'contact'],
      },
    });
  });
  const data = `${o.data === 'links' ? '/data/hb/contacts' : `/data/${o.data}`}${o.filter ?? ''}`;
  b.views.view('v', (v) => {
    const params = {
      ...(o.blocks ? { of: hba } : { of: '/data/hb' as never }),
      ...(o.copies ? { copies: o.copies } : {}),
    };
    v.section('s', params, (s) => {
      if (o.blocks) {
        s.block('a1', { of: hba });
        if (o.blocks === 'zoom') s.block('zoom', { of: hba, range: [5, 15] });
        s.block('a2', { of: hba });
      }
      s.track('letters').layer('letters', { type: 'letters' });
      s.track('pairs').layer('pairs', { type: 'pairs', data: data as never });
    });
  });
  return b.build();
}

/** Settled once a spec's options (fixtures once per file). */
const Settles = new Map<string, Promise<Settled>>();
const settled = (o: Options) => {
  const key = JSON.stringify(o);
  let s = Settles.get(key);
  if (!s) Settles.set(key, (s = settle(spec(o))));
  return s;
};

function pairsLay(f: Frame): PairsLay {
  for (const l of f.layout.sections[0].sectionLays) if (l.type.name === 'pairs') return l.lay as PairsLay;
  throw new Error('no pairs layer');
}

const Pt = new PlacePoint();

/** A display unit's page x. */
function xOf(f: Frame, u: number): number {
  expect(f.layout.sections[0].placement.place(u, Pt)).toBe(true);
  return Pt.x;
}

/** The display unit of position `p` on the block whose system is `system` (a chip: its middle). */
function unitOf(f: Frame, system: string, p: number): number {
  const blocks = f.layout.sections[0].section.blocks!.filter(
    (b) => b.axis.system === system || b.systems.includes(system),
  );
  const b = blocks[0];
  return b.form === 'whole' ? b.offset! + 0.5 : b.offset! + p - b.axis.ranges![0] + 0.5;
}

/** Each arc: its pair's key and its feet's page x, sorted. */
function feet(f: Frame): { key: string | number; x: [number, number] }[] {
  const lay = pairsLay(f);
  const d = lay.prep!.data!;
  const out = [];
  for (let i = 0; i < lay.n; i++) {
    const c = lay.curve;
    const x = [c[8 * i], c[8 * i + 6]].sort((p, q) => p - q) as [number, number];
    out.push({ key: d.keys[lay.pair[i]], x });
  }
  return out;
}

/** Feet expected from display units, sorted. */
const want = (f: Frame, key: string | number, ua: number, ub: number) => ({
  key,
  x: [xOf(f, ua), xOf(f, ub)].sort((p, q) => p - q),
});

/** Each arc hit at its middle as its pair; the SVG has each arc's curve (canvas and SVG alike). */
function hitAndSvg(f: Frame) {
  const lay = pairsLay(f);
  const curves = [...drawListsToSvg(f.lists).matchAll(/d="M ([^"C]*) C ([^"C]*)"/g)].map((m) =>
    `${m[1]} ${m[2]}`.split(' ').map(Number),
  );
  for (let i = 0; i < lay.n; i++) {
    const c = lay.curve.subarray(8 * i, 8 * i + 8);
    const x = (c[0] + 3 * c[2] + 3 * c[4] + c[6]) / 8;
    const y = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
    expect(hitView(f.layout, x, y)?.element?.key).toBe(lay.prep!.data!.keys[lay.pair[i]]);
    // the SVG's numbers rounded to hundredths (the draw lists' points are 32-bit)
    expect(curves.some((q) => q.every((v, j) => Math.abs(v - c[j]) < 0.02))).toBe(true);
  }
}

describe('a pair’s ends on copies (EC-48, EC-65)', () => {
  for (const copies of ['each', 'merge'] as const)
    it(`(1, 6) the links name their units: each heme's contact on its own chain (${copies})`, async () => {
      const f = await (await settled({ data: 'links', copies })).frame({ size: Size });
      expect(f.layout.issues).toEqual([]);
      // merged: a heme ×4 and the α and β ×2 blocks, each contact on the block its unit is merged into
      const keys = pairsLay(f).prep!.data!.keys;
      const expected = Contacts.heme.map((h, i) =>
        want(f, keys[i], unitOf(f, `/data/hb/${h}`, 0), unitOf(f, `/data/hb/${Contacts.chain[i]}`, Contacts.at[i])),
      );
      expect(feet(f)).toEqual(expected);
      hitAndSvg(f);
    });

  it('(1) contacts on the molecules naming no units, split: on no copy, never a guessed one, said once', async () => {
    const f = await (await settled({ data: 'across', copies: 'each' })).frame({ size: Size });
    expect(feet(f)).toEqual([]);
    expect(pairsLay(f).carets).toEqual([]);
    expect(f.layout.issues.map((i) => [i.code, i.severity, i.message])).toEqual([
      [
        'seq.render.pairs.copy',
        'warning',
        '2 pairs on split copies (hem ×4, hba ×2) not drawn: the data names no units',
      ],
    ]);
  });

  it('(1) contacts on the molecules, merged: on the heme ×4 and α ×2 blocks', async () => {
    const f = await (await settled({ data: 'across' })).frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    expect(feet(f)).toEqual(
      Across.map((r, i) => want(f, ['h5', 'h7'][i], unitOf(f, '/data/hb/ha', 0), unitOf(f, '/data/hb/A', r))),
    );
  });

  it('(2) a pair within one split system draws on each copy: α’s units A and C', async () => {
    const f = await (await settled({ data: 'intra', copies: 'each' })).frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    // a pair, then the pair again on each copy
    const expected = Intra.flatMap(([a, b]) =>
      ['/data/hb/A', '/data/hb/C'].map((unit) => want(f, `${a}-${b}`, unitOf(f, unit, a), unitOf(f, unit, b))),
    );
    expect(feet(f)).toEqual(expected);
    hitAndSvg(f);
  });

  it('(2) merged, a pair within one system is one arc a pair', async () => {
    const f = await (await settled({ data: 'intra' })).frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    expect(feet(f)).toEqual(
      Intra.map(([a, b]) => want(f, `${a}-${b}`, unitOf(f, '/data/hb/A', a), unitOf(f, '/data/hb/A', b))),
    );
  });

  it('(2) a sequence’s blocks split: each copy draws its pairs; a zoom block is no copy', async () => {
    for (const blocks of ['twice', 'zoom'] as const) {
      const f = await (await settled({ data: 'intra', copies: 'each', blocks })).frame({ size: Size });
      expect(f.layout.issues).toEqual([]);
      const sec = f.layout.sections[0].section.blocks!;
      const copies = sec.filter((b) => b.copy === sec[0].copy);
      expect(copies.length).toBe(2);
      const expected = Intra.flatMap(([p, q]) =>
        copies.map((b) => want(f, `${p}-${q}`, b.offset! + p + 0.5, b.offset! + q + 0.5)),
      );
      expect(feet(f)).toEqual(expected);
    }
  });

  it('(2) the default merge mode with a zoom block: each pair once, on the first block, no issue', async () => {
    const f = await (await settled({ data: 'intra', blocks: 'zoom' })).frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    expect(f.layout.sections[0].section.blocks!.length).toBe(2);
    expect(feet(f)).toEqual(Intra.map(([p, q]) => want(f, `${p}-${q}`, p + 0.5, q + 0.5)));
  });

  it('(3) the issue counts the pairs shown, not those the filter leaves out, and names both systems', async () => {
    const filter = '/(except :of (item :type far))';
    const f = await (await settled({ data: 'across', copies: 'each', filter })).frame({ size: Size });
    expect(f.layout.issues.map((i) => i.message)).toEqual([
      '1 pair on split copies (hem ×4, hba ×2) not drawn: the data names no units',
    ]);
    expect(pairsLay(f).prep!.filtered).toBe(1);
  });

  it('(3) a pair within one system the filter leaves out is drawn on no copy', async () => {
    const filter = '/(except :of (item :type salt))';
    const f = await (await settled({ data: 'intra', copies: 'each', filter })).frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    expect(feet(f).map((a) => a.key)).toEqual(['4-12', '4-12']);
    expect(pairsLay(f).prep!.filtered).toBe(1);
  });

  it('(4) an end out of its block’s range beside a copy end: its mark on each copy, or `pairs.systems`', async () => {
    // x's 0–8 twice, l (one position) twice: pairs within x and from l to x, an end at 12 out of every block's range
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'protein', letters: Hba });
    const l = b.data.sequence('l', { alphabet: 'protein', letters: 'A' });
    const within = b.data.pairs('within', {
      a: x,
      b: x,
      rows: { a_start: [2], a_end: [3], b_start: [12], b_end: [13], id: ['w'] },
    });
    const across = b.data.pairs('across', {
      a: l,
      b: x,
      rows: { a_start: [0, 0], a_end: [1, 1], b_start: [3, 12], b_end: [4, 13], id: ['in', 'out'] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: x, copies: 'each' }, (s) => {
        for (const id of ['x1', 'x2']) s.block(id, { of: x, range: [0, 8] });
        for (const id of ['l1', 'l2']) s.block(id, { of: l });
        s.track('within').layer('within', { type: 'pairs', data: within });
        s.track('across').layer('across', { type: 'pairs', data: across });
      }),
    );
    const f = await (await settle(b.build())).frame({ size: Size });
    const [w, a] = f.layout.sections[0].sectionLays.map((q) => q.lay as PairsLay);
    // within one system: the placed end's caret on each copy (a mark, no issue: its partner is not in view)
    expect(w.n).toBe(0);
    expect(w.carets.length / 4).toBe(2);
    expect(new Set([w.carets[0], w.carets[4]]).size).toBe(2);
    // across two systems: no copy guessed, no caret; the end out of range said by `pairs.systems`
    expect(a.n).toBe(0);
    expect(a.carets).toEqual([]);
    expect(f.layout.issues.map((i) => [i.address, i.code, i.message])).toEqual([
      [a.key, 'seq.render.pairs.copy', '2 pairs on split copies (l ×2, x ×2) not drawn: the data names no units'],
      [a.key, 'seq.render.pairs.systems', 'Some pairs’ ends are on no system the section shows'],
    ]);
  });

  it('across a tile edge: each copy’s arcs at their data, the letters read by region', async () => {
    const b = createSeqSpec();
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const p = b.data.pairs('p', {
      a: chr,
      b: chr,
      rows: { a_start: [1020, 1010], a_end: [1021, 1011], b_start: [1030, 1060], b_end: [1031, 1061], id: ['p', 'q'] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: chr, copies: 'each' }, (s) => {
        for (const id of ['c1', 'c2']) s.block(id, { of: chr, range: [1000, 1100] });
        s.track('letters').layer('letters', { type: 'letters' });
        s.track('pairs').layer('pairs', { type: 'pairs', data: p });
      }),
    );
    const s = await settle(b.build(), [BigRuntime]);
    // the second copy's 1018–1036: across the tile edge at 1024
    const f = await s.frame({ size: Size, viewport: { start: 118, end: 136 } });
    expect(f.layout.issues).toEqual([]);
    // p within the window on the second copy (the first's taken past the window's edge): its feet at its data
    expect(feet(f).filter((a) => a.key === 'p')).toEqual([want(f, 'p', 100 + 20.5, 100 + 30.5)]);
    hitAndSvg(f);
  });

  it('zoomed onto copy C: its arcs at their data, an end out of the window a caret, hit as its pair', async () => {
    const s = await settled({ data: 'intra', copies: 'each' });
    const at = (u: string) => f0.layout.sections[0].section.blocks!.find((b) => b.axis.system === u)!.offset!;
    const f0 = await s.frame({ size: Size });
    // C's 2 to 11: pair 2–10 inside, 4–12 with one end out of the window, A's ends out of it
    const c = at('/data/hb/C');
    const f = await s.frame({ size: Size, viewport: { start: c + 1, end: c + 11 } });
    expect(f.layout.issues).toEqual([]);
    const lay = pairsLay(f);
    expect(feet(f).map((a) => a.key)).toContain('2-10');
    const arc = feet(f).find((a) => a.key === '2-10')!;
    expect(arc).toEqual(want(f, '2-10', c + 2.5, c + 10.5));
    // an end past the window is taken `reach` past its edge, the rest a caret, hit as its pair
    for (let i = 0; i + 3 < lay.carets.length; i += 4) {
      const y = lay.carets[i + 1] + (lay.carets[i + 2] > 0 ? -3 : 3);
      expect(hitView(f.layout, lay.carets[i], y)?.element?.key).toBe(lay.prep!.data!.keys[lay.carets[i + 3]]);
    }
  });

  it('(5) a caret is hit as its pair, a tie to the later drawn; a second layout reads the ends kept by the section', async () => {
    // α's pairs with one end on no block (α's 0–8 twice): x and y at one place, z beside them
    const b = createSeqSpec();
    const hba = b.data.sequence('hba', { alphabet: 'protein', letters: Hba });
    const p = b.data.pairs('p', {
      a: hba,
      b: hba,
      rows: { a_start: [2, 2, 4], a_end: [3, 3, 5], b_start: [10, 12, 14], b_end: [11, 13, 15], id: ['x', 'y', 'z'] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: hba, copies: 'each' }, (s) => {
        s.block('a1', { of: hba, range: [0, 8] });
        s.block('a2', { of: hba, range: [0, 8] });
        s.track('pairs').layer('pairs', { type: 'pairs', data: p });
      }),
    );
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    const lay = pairsLay(f);
    expect(lay.n).toBe(0);
    // a caret a pair on each copy: x's and y's at one place hit as y (drawn later), z's as z
    expect(lay.carets.length / 4).toBe(6);
    const keys = [];
    for (let i = 0; i + 3 < lay.carets.length; i += 4) {
      const y = lay.carets[i + 1] + (lay.carets[i + 2] > 0 ? -3 : 3);
      keys.push(hitView(f.layout, lay.carets[i], y)?.element?.key);
    }
    expect(keys).toEqual(['y', 'y', 'y', 'y', 'z', 'z']);
    // the ends are kept by the section: a second layout reads them again
    const prep = lay.prep!;
    const ends = prep.units.get(f.layout.sections[0].section);
    await s.frame({ size: Size });
    expect(prep.units.get(f.layout.sections[0].section)).toBe(ends);
  });
});
