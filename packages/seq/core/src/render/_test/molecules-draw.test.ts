/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A complex drawn against its data (_spec/README.md "Render", Blocks; Layer types: letters, The pairs layer): every
// chip its label in its pill (a molecule `whole`, its residues one slot, never a position's share), every modified
// position in view its mark (at the fit view, zoomed, cut by the window's edge), every link with both ends in view its
// arc and its hit box (ends on per-row units, an end on a merged copy drawn on its block); canvas and SVG alike.

import { describe, expect, it } from 'vitest';
import {
  type BlockPlacement,
  DarkTokens,
  drawListsToSvg,
  type Frame,
  hitView,
  LightTokens,
  type LettersLay,
  type PairsLay,
  type RangesLay,
  recordDrawLists,
} from '../index.js';
import { settle } from './helpers.js';
import { moleculesSpec } from './molecules-fixtures.js';

const Size = { width: 900, height: 400 };

const texts = (f: Frame, pass = 'data') =>
  recordDrawLists(f.lists, [pass as 'data'])
    .split('\n')
    .filter((l) => l.startsWith('text '))
    .map((l) => /^text "([^"]*)"/.exec(l)![1]);

/** Each strip's letters lay in a section. */
function lettersLays(f: Frame, section = 0): { strip: number; lay: LettersLay }[] {
  const s = f.layout.sections[section];
  const nb = s.placement.strips.bands;
  const out: { strip: number; lay: LettersLay }[] = [];
  for (let i = 0; i < s.placement.strips.n; i++)
    for (let j = 0; j < nb; j++)
      for (const l of s.lays[i * nb + j] ?? [])
        if (l.type.name === 'letters') out.push({ strip: i, lay: l.lay as LettersLay });
  return out;
}

function pairsLay(f: Frame): PairsLay {
  for (const l of f.layout.sections[0].sectionLays) if (l.type.name === 'pairs') return l.lay as PairsLay;
  throw new Error('no pairs layer');
}

describe('a complex drawn against its data', () => {
  for (const tokens of [LightTokens, DarkTokens])
    it(`every chip its label in its pill, a molecule one unit (${tokens === LightTokens ? 'light' : 'dark'})`, async () => {
      const s = await settle(moleculesSpec());
      s.renderer.tokens = tokens;
      const f = await s.frame({ size: Size });
      const p = f.layout.sections[0].placement as BlockPlacement;
      const blocks = f.layout.sections[0].section.blocks!;
      const drawn = texts(f);
      const svg = drawListsToSvg(f.lists);
      for (let i = 0; i < p.strips.n; i++) {
        if (!p.strips.chip[i]) continue;
        const label = blocks[p.strips.block[i]].label;
        expect(drawn).toContain(label);
        expect(svg).toContain(`>${label}<`);
      }
      expect(blocks.filter((b) => b.form === 'whole').map((b) => b.label)).toEqual(['Glycan', 'ATP', 'NA']);
      // the glycan's five residues are one chip, one unit of the section (layouts.md): its label's width, whatever the
      // zoom; its strip shows its block's every residue (its pieces)
      const g = p.strips.axes.findIndex((a) => a.system === '/data/c/g');
      expect([p.strips.axes[g].u0, p.strips.axes[g].u1]).toEqual([0, 1]);
      expect(Array.from(p.strips.axes[g].pieces)).toEqual([0, 5]);
      expect(blocks.find((b) => b.label === 'Glycan')!.axis.ranges).toEqual([0, 5]);
      const zoomed = await s.frame({ size: Size, viewport: { start: 60, end: 89 } });
      const pz = zoomed.layout.sections[0].placement as BlockPlacement;
      const gz = pz.strips.axes.findIndex((a) => a.system === '/data/c/g');
      expect(pz.strips.axes[gz].width).toBeCloseTo(p.strips.axes[g].width, 6);
    });

  it('a chip is hit as its whole block: its every residue', async () => {
    const s = await settle(moleculesSpec());
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const g = p.strips.axes.findIndex((a) => a.system === '/data/c/g');
    const a = p.strips.axes[g];
    const y = p.strips.bandY[g * p.strips.bands] + 2;
    const h = hitView(f.layout, a.x0 + a.width - 2, y)!;
    expect(h.system).toBe('/data/c/g');
    expect(h.extent?.parts.map((q) => [q.system, Array.from(q.ranges)])).toEqual([['/data/c/g', [0, 5]]]);
  });

  for (const viewport of [undefined, { start: 25, end: 31.5 }])
    it(`every modified position in view has its mark (${viewport ? 'cut by the window' : 'fit'})`, async () => {
      const s = await settle(moleculesSpec());
      const f = await s.frame({ size: Size, ...(viewport ? { viewport } : {}) });
      const [fc] = lettersLays(f).filter((l) => l.lay.system === '/data/c/A');
      // MSE at 0, SEP at 20, CRO over 30–32: those in view
      const want = viewport ? [30] : [0, 20, 30];
      const rows = [];
      for (let i = 0; i + 2 < fc.lay.mods.length; i += 3) rows.push(fc.lay.mods[i + 2]);
      expect(rows.map((r) => [0, 20, 30][r])).toEqual(want);
      // drawn: a rect each over the cells (a dot, or the CRO bar)
      const rects = recordDrawLists(f.lists, ['data'])
        .split('\n')
        .filter((l) => l.startsWith('rect '));
      expect(rects.length).toBeGreaterThan(want.length);
    });

  it('the tRNA: a mark per modified position, alternatives one mark; `code` prints the component where it fits', async () => {
    const s = await settle(moleculesSpec({ modified: 'code' }));
    const f = await s.frame({ size: Size });
    const [t] = lettersLays(f, 1);
    expect(t.lay.mods.filter((_, i) => i % 3 === 2)).toEqual([0, 1]);
    // 25 letters over 900 px: every cell has room for three letters
    expect(texts(f)).toEqual(expect.arrayContaining(['2MG', 'H2U']));
  });

  it('every link with both ends in view has its arc and its hit box; the ends on per-row units', async () => {
    const s = await settle(moleculesSpec({ copies: 'each' }));
    const f = await s.frame({ size: Size });
    const lay = pairsLay(f);
    // A10–glycan, A4–B4, ATP–A50: all three drawn
    expect(lay.n).toBe(3);
    for (let k = 0; k < lay.n; k++) {
      const c = lay.curve.subarray(8 * k, 8 * k + 8);
      // the curve's middle (a cubic at t = 1/2)
      const x = (c[0] + 3 * c[2] + 3 * c[4] + c[6]) / 8;
      const y = (c[1] + 3 * c[3] + 3 * c[5] + c[7]) / 8;
      const h = hitView(f.layout, x, y);
      expect(h?.element?.owner).toBe('/data/c/links');
    }
  });

  it('merged copies: an end on a copy is drawn on the block it is merged into', async () => {
    const s = await settle(moleculesSpec());
    const f = await s.frame({ size: Size });
    expect(pairsLay(f).n).toBe(3);
  });

  it("in text, the molecules' chips flow like words on one line, each labeled (EC-99)", async () => {
    const s = await settle(moleculesSpec({ layout: 'text' }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement;
    const chips: number[] = [];
    for (let i = 0; i < p.strips.n; i++) if (p.strips.chip[i]) chips.push(i);
    expect(chips.length).toBe(3);
    // one line: the same top, left to right
    expect(new Set(chips.map((i) => p.strips.top[i])).size).toBe(1);
    for (let k = 1; k < chips.length; k++)
      expect(p.strips.axes[chips[k]].x0).toBeGreaterThan(p.strips.axes[chips[k - 1]].x0);
    expect(texts(f)).toEqual(expect.arrayContaining(['Glycan', 'ATP', 'NA']));
  });

  it('a ranges layer over the links: each end a mark at its site, labeled by the other end', async () => {
    const s = await settle(moleculesSpec({ copies: 'each', sites: true }));
    const f = await s.frame({ size: Size });
    expect(f.layout.sections[0].section.blocks!.length).toBe(5);
    const sec = f.layout.sections[0];
    const nb = sec.placement.strips.bands;
    const marks = new Map<string, number>();
    for (let i = 0; i < sec.placement.strips.n; i++)
      for (const l of sec.lays[i * nb + 1] ?? [])
        if (l.type.name === 'ranges') marks.set(sec.placement.strips.axes[i].system, (l.lay as RangesLay).items);
    // A: the glycan's site, the disulfide's end, the metal's site; B: the disulfide's other end; the glycan; ATP
    expect(Object.fromEntries(marks)).toEqual({ '/data/c/A': 3, '/data/c/B': 1, '/data/c/g': 1, '/data/c/l': 1 });
    expect(s.view.issues.filter((i) => i.code !== 'seq.view.layer.unplaced')).toEqual([]);
    // zoomed onto N11: its mark wide enough for its label, the other end
    const z = await s.frame({ size: Size, viewport: { start: 6, end: 16 } });
    expect(texts(z, 'labels')).toEqual(expect.arrayContaining(['→ Glycan 1']));
  });
});
