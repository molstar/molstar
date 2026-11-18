/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis drawn (_spec/README.md "Render", Tracks layout, Layer types `letters`, Lanes and Hits): a reversed axis's
// transform both ways, letters through a codon path and a one-to-one path, the complement, breaks and their hits, the
// direction mark; hit-testing the inverse of placement over random sections.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { AxisTransform, Cut } from '../axis.js';
import { rangesLocation } from '@molstar/seqx-core/query';
import { hitView, type LettersLay, LightTokens, type RangesLay, recordDrawLists, renderSvg } from '../index.js';
import { settle } from './helpers.js';

const Size = { width: 640, height: 200 };

// a protein of 4 residues on a 30-nucleotide gene, exons [2, 7) and [12, 19): the second codon (minus strand: the
// third) is split by the intron
const Gene = 'ACGTACGGATCCGTAGCATGCATCGATCGA';

function spec(gene: Record<string, unknown>, strand: '+' | '-' = '-', layer: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
  const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
  b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [7, 19] }, strand });
  b.views.view('v', (v) => {
    const s = v.section('gene', { of: g, ...(gene as object) });
    s.track('residues').layer('l', { type: 'letters', data: p, ...(layer as object) });
    s.track('bases').layer('l', { type: 'letters' });
  });
  return b.build();
}

describe('the axis transform', () => {
  it('a reversed piece shows `end − 1 − k` at offset `k`, spans mirrored with their cut ends on screen', () => {
    const a = new AxisTransform('/s', [12, 19, 2, 7], 0, 120, 1, undefined, { dir: -1 });
    expect(a.positionAt(0.5)).toBe(18);
    expect(a.positionAt(119.5)).toBe(2);
    expect(a.cell(18)).toBe(true);
    expect([a.sx0, a.sx1]).toEqual([0, 10]);
    expect(a.visible([])).toEqual([12, 19, 2, 7]);
    // a range over the intron: a piece each side, its start (5) on the right, cut on the left by the break
    expect(a.spans(5, 14, [])).toEqual([50, 70, Cut.right, 70, 90, Cut.left]);
    expect(a.breaks([])).toEqual([70, 7, 12]);
    // a range ending on a piece's edge at a break is cut there (a projected band, two ranges across the intron)
    expect(a.spans(12, 14, [])).toEqual([50, 70, Cut.right]);
    expect(a.spans(14, 16, [])).toEqual([30, 50, 0]);
  });

  it('hit-testing is the inverse of placement over random sections, both directions', () => {
    const rng = new PCG(11);
    for (let t = 0; t < 100; t++) {
      const pieces: number[] = [];
      let at = 0;
      for (let k = 0; k < 1 + Math.floor(rng.float() * 6); k++) {
        at += Math.floor(rng.float() * 5);
        const len = 1 + Math.floor(rng.float() * 8);
        pieces.push(at, at + len);
        at += len;
      }
      const dir = rng.float() < 0.5 ? 1 : -1;
      const display: number[] = [];
      if (dir === 1) display.push(...pieces);
      else for (let i = pieces.length - 2; i >= 0; i -= 2) display.push(pieces[i], pieces[i + 1]);
      const a = new AxisTransform('/s', display, 10, 300, 2, undefined, { dir });
      for (const [s, e] of chunks(a.visible([])))
        for (let p = s; p < e; p++) {
          expect(a.cell(p)).toBe(true);
          expect(a.positionAt((a.sx0 + a.sx1) / 2)).toBe(p);
        }
    }
  });
});

describe('letters through a codon path', () => {
  it('a residue a codon wide on the gene’s coding axis, along the protein, each part of a split codon a cell', async () => {
    const s = await settle(spec({ direction: { along: '/data/p' }, pieces: 'coding' }));
    expect(s.view.issues).toEqual([]);
    const f = await s.frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as LettersLay;
    // residues 0, 1, 2 and 3, left to right; residue 2 split by the intron: two cells (2 nucleotides and 1)
    const cells = Array.from({ length: lay.n }, (_, j) => [lay.pos[j], lay.x1[j] - lay.x0[j]]);
    expect(cells.map(([p]) => p)).toEqual([0, 1, 2, 2, 3]);
    expect(lay.x0[0]).toBeLessThan(lay.x0[4]);
    const w = cells[0][1];
    expect(cells[2][1] + cells[3][1]).toBeCloseTo(w, 0);
    // the letters: M, E, E, E (both parts with room at this zoom), P
    expect(Array.from(lay.code.subarray(0, lay.n), (c) => String.fromCharCode(c)).join('')).toBe('MEEEP');
  });

  it('a hit on a residue carries its extent on the protein and its codon as the hull', async () => {
    const s = await settle(spec({ direction: { along: '/data/p' }, pieces: 'coding' }));
    const f = await s.frame({ size: Size });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as LettersLay;
    const j = 2;
    const h = hitView(f.layout, (lay.x0[j] + lay.x1[j]) / 2, lay.box.y + 2)!;
    expect(h.system).toBe('/data/g');
    expect(h.extent?.parts.map((p) => [p.system, [...p.ranges]])).toEqual([['/data/p', [2, 3]]]);
    // residue 2 on the minus strand: nucleotides 12, 6, 5
    expect(h.hull).toEqual([5, 13]);
  });

  it('refuses the inverse placement with its reason; a many-to-one path too', async () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
    b.data.map('m', { from: p, to: g, method: 'codon', exons: { start: [2, 12], end: [7, 19] }, strand: '+' });
    b.views.view('v', (v) => void v.section('protein', { of: p }).track('t').layer('l', { type: 'letters', data: g }));
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    expect(f.layout.issues.map((i) => i.code)).toEqual(['seq.render.layer.letters']);
    expect(f.layout.issues[0].message).toContain('a third of a residue');
  });

  it('draws the conflict ink on a residue its codon does not translate to, once the gene’s letters are held', async () => {
    const s = await settle(spec({ pieces: 'coding' }, '+'));
    const f = await s.frame({ size: Size });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as LettersLay;
    const conflict = LightTokens.colors['seq-conflict'];
    // the toy gene translates to none of MEEP: every residue is in conflict
    for (let j = 0; j < lay.n; j++) expect(lay.fill[j]).toBe(conflict);
  });

  it('three-letter codes where the narrowest part in view has room (`code: auto`), or as written', async () => {
    const wide = await settle(spec({ pieces: 'coding' }, '+', { code: 'auto' }));
    const zoomed = await wide.frame({ size: Size, viewport: { start: 0, end: 6 } });
    expect((zoomed.layout.sections[0].tracks[0].layers[0].lay as LettersLay).three).toBe(true);
    expect(recordDrawLists(zoomed.lists, ['data'])).toContain('Met');
    const fit = await wide.frame({ size: { width: 200, height: 200 } });
    expect((fit.layout.sections[0].tracks[0].layers[0].lay as LettersLay).three).toBe(false);
    const one = await settle(spec({ pieces: 'coding' }, '+', { code: 'one' }));
    const z = await one.frame({ size: Size, viewport: { start: 0, end: 6 } });
    expect((z.layout.sections[0].tracks[0].layers[0].lay as LettersLay).three).toBe(false);
  });
});

describe('letters through a one-to-one path', () => {
  it('a letter its image’s cell; a position no letter lands on, a gap mark', async () => {
    const b = createSeqSpec();
    const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
    const h = b.data.sequence('h', { alphabet: 'dna', letters: 'ACGTACGTAC' });
    b.data.map('o', { from: h, to: g, method: 'offset', by: 3 });
    b.views.view('v', (v) => void v.section('g', { of: g }).track('t').layer('l', { type: 'letters', data: h }));
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    expect(f.layout.issues).toEqual([]);
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as LettersLay;
    expect(lay.n).toBe(10);
    const a = f.layout.sections[0].axis!;
    a.cell(3);
    expect(lay.x0[0]).toBeCloseTo(a.snap(a.sx0), 6);
    // gaps: [0, 3) and [13, 30)
    expect(lay.gaps.length).toBe(6);
  });
});

describe('the complement', () => {
  // the letters left to right on screen
  const bases = (lay: LettersLay) =>
    Array.from({ length: lay.n }, (_, j) => j)
      .sort((a, b) => lay.x0[a] - lay.x0[b])
      .map((j) => String.fromCharCode(lay.code[j]))
      .join('');

  it('a reversed DNA axis shows the other strand read left to right (`auto`); `never` and `always` as written', async () => {
    const fwd = await settle(spec({}));
    const rev = await settle(spec({ direction: 'reversed' }));
    const a = (await fwd.frame({ size: Size })).layout.sections[0].tracks[1].layers[0].lay as LettersLay;
    const b = (await rev.frame({ size: Size })).layout.sections[0].tracks[1].layers[0].lay as LettersLay;
    expect(bases(a)).toBe(Gene);
    // reversed: the positions from 29 down, each its complement: the reverse complement
    const rc = [...Gene]
      .reverse()
      .map((c) => ({ A: 'T', C: 'G', G: 'C', T: 'A' })[c])
      .join('');
    expect(bases(b)).toBe(rc);
  });

  it('`never` keeps the stored letters on a reversed axis; `always` complements a forward one', async () => {
    const b = createSeqSpec();
    const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
    b.views.view('v', (v) => {
      v.section('r', { of: g, direction: 'reversed' }).track('t').layer('l', { type: 'letters', complement: 'never' });
      v.section('f', { of: g }).track('t').layer('l', { type: 'letters', complement: 'always' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size });
    const r = f.layout.sections[0].tracks[0].layers[0].lay as LettersLay;
    const fw = f.layout.sections[1].tracks[0].layers[0].lay as LettersLay;
    expect(bases(r)).toBe([...Gene].reverse().join(''));
    expect(bases(fw)).toBe([...Gene].map((c) => ({ A: 'T', C: 'G', G: 'C', T: 'A' })[c]).join(''));
  });
});

describe('breaks and the direction mark', () => {
  it('a hairline a break across the band, a notch above it, hit within 3 px of the notch; cells beside it win', async () => {
    const s = await settle(spec({ follow: undefined, pieces: 'coding', direction: { along: '/data/p' } }));
    const f = await s.frame({ size: Size });
    const sec = f.layout.sections[0];
    expect(sec.breaks.length).toBe(3);
    const [x, hs, he] = sec.breaks;
    expect([hs, he]).toEqual([7, 12]);
    const band = sec.band!;
    const h = hitView(f.layout, x + 2, band.y - 2)!;
    expect(h.target).toBe('break');
    expect(h.extent?.parts.map((p) => [...p.ranges])).toEqual([[7, 12]]);
    expect(hitView(f.layout, x + 5, band.y - 2)).toBeUndefined();
    // on data beside it, the cells
    expect(hitView(f.layout, x + 1, sec.tracks[1].box.y + 2)?.target).toBeUndefined();
    const data = recordDrawLists(f.lists, ['data']);
    expect(data).toContain('#8a847a');
  });

  it('a reversed axis says so in its first label, `along` names what it reads along; `none` leaves it out', async () => {
    const along = await settle(spec({ pieces: 'coding', direction: { along: '/data/p' } }));
    expect(recordDrawLists((await along.frame({ size: Size })).lists, ['labels'])).toContain('along p');
    const rev = await settle(spec({ direction: 'reversed' }));
    expect(recordDrawLists((await rev.frame({ size: Size })).lists, ['labels'])).toContain('← reversed');
    const none = await settle(spec({ direction: 'reversed' }), [], { axis: { directionMark: 'none' } });
    expect(recordDrawLists((await none.frame({ size: Size })).lists, ['labels'])).not.toContain('reversed');
    const fwd = await settle(spec({}));
    expect(recordDrawLists((await fwd.frame({ size: Size })).lists, ['labels'])).not.toContain('←');
  });
});

describe('what a break hides is not dropped', () => {
  function hiddenSpec() {
    const b = createSeqSpec();
    const p = b.data.sequence('p', { alphabet: 'protein', letters: 'MEEP' });
    const g = b.data.sequence('g', { alphabet: 'dna', letters: Gene });
    const m = b.data.map('m', {
      from: p,
      to: g,
      method: 'codon',
      exons: { start: [2, 12], end: [7, 19] },
      strand: '+',
    });
    // one item in an exon, two wholly in the intron [7, 12), one across it
    const r = b.data.ranges('sites', {
      target: g,
      rows: { start: [3, 8, 9, 6], end: [5, 9, 11, 14], id: ['a', 'b', 'c', 'd'] },
    });
    b.views.view(
      'v',
      (v) => void v.section('g', { of: g, pieces: m }).track('t').layer('l', { type: 'ranges', data: r }),
    );
    return b.build();
  }

  it('items wholly in a hidden stretch draw one tick on its break a lane, hit as the break with their count', async () => {
    const s = await settle(hiddenSpec());
    const f = await s.frame({ size: Size });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    // `b` and `c` are hidden: in two lanes or one, a tick a lane; `a` and `d` are marks
    expect(lay.tn).toBeGreaterThanOrEqual(1);
    let count = 0;
    for (let t = 0; t < lay.tn; t++) count += lay.ticks[5 * t + 2];
    expect(count).toBe(2);
    const x = lay.ticks[0];
    const y = lay.box.y + (lay.ticks[1] + 0.5) * lay.laneHeight;
    const h = hitView(f.layout, x + 1, y)!;
    expect(h.target).toBe('break');
    expect(h.extent?.parts.map((p) => [...p.ranges])).toEqual([[7, 12]]);
    expect(h.count).toBeGreaterThanOrEqual(1);
  });

  it('a Location wholly hidden (a selection) lights its break', async () => {
    const s = await settle(hiddenSpec());
    const f = await s.frame({
      size: Size,
      emphasis: { selection: rangesLocation('/data/g', [8, 10]) },
    });
    const x = f.layout.sections[0].breaks[0];
    const rects = recordDrawLists(f.lists, ['selection'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => Number(l.split(' ')[1]));
    expect(rects.some((r) => Math.abs(r - x) <= 1)).toBe(true);
  });
});

describe('the export', () => {
  it('says each section’s direction and pieces; a viewport in positions needs no knowledge of the axis', async () => {
    const s = spec({ pieces: 'coding', direction: { along: '/data/p' } });
    const r = await renderSvg(s, { size: Size, viewport: { positions: { '/views/v/gene': [12, 15] } } });
    expect(r.issues).toEqual([]);
    expect(r.svg).toContain('data-section="/views/v/gene" data-direction="reversed" data-pieces="coding"');
    // positions 12 to 15 on the reversed coding axis are its last offsets but for the lower exon: [4, 7)
    const offsets = await renderSvg(s, { size: Size, viewport: { start: 4, end: 7 } });
    expect(recordDrawLists(r.lists, ['data'])).toBe(recordDrawLists(offsets.lists, ['data']));
  });
});

function chunks(v: readonly number[]): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i + 1 < v.length; i += 2) out.push([v[i], v[i + 1]]);
  return out;
}
