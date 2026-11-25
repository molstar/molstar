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
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
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
      // the glycan's five residues are one chip, one unit of the section (_spec/README.md "Render", Section layouts):
      // its label's width, whatever the zoom; its strip shows its block's every residue (its pieces)
      const g = p.strips.axes.findIndex((a) => a.system === '/data/c/g');
      expect([p.strips.axes[g].u0, p.strips.axes[g].u1]).toEqual([0, 1]);
      expect(Array.from(p.strips.axes[g].pieces)).toEqual([0, 5]);
      expect(blocks.find((b) => b.label === 'Glycan')!.axis.ranges).toEqual([0, 5]);
      const zoomed = await s.frame({ size: Size, viewport: { start: 60, end: 89 } });
      const pz = zoomed.layout.sections[0].placement as BlockPlacement;
      const gz = pz.strips.axes.findIndex((a) => a.system === '/data/c/g');
      expect(pz.strips.axes[gz].width).toBeCloseTo(p.strips.axes[g].width, 6);
    });

  it('a ligand chip colored on its own: a theme over the molecule’s fields (component, class), each chip its color', async () => {
    const color = theme.rules([
      theme.rule([theme.eq('component', 'ATP')], '#d62728'),
      theme.rule([theme.eq('class', 'ion')], '#9467bd'),
    ]);
    const s = await settle(moleculesSpec({ color }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const pill = (system: string) => {
      const i = p.strips.axes.findIndex((a) => a.system === system);
      return lettersLays(f).find((l) => l.strip === i)!.lay.fill[0];
    };
    expect([pill('/data/c/l'), pill('/data/c/n')]).toEqual([0xd62728ff, 0x9467bdff]);
    // a chip no rule names takes the rules' default, neither rule's color
    expect([0xd62728ff, 0x9467bdff]).not.toContain(pill('/data/c/g'));
    expect(recordDrawLists(f.lists, ['data'])).toContain('#d62728ff');
    // the polymers' residues keep the default scheme where letters show: rules no clause of which holds give its
    // `else`, the layer type's default
    const zoom = { size: Size, viewport: { start: 0, end: 20 } };
    const fills = async (spec: Parameters<typeof settle>[0]) => {
      const z = await (await settle(spec)).frame(zoom);
      const a = lettersLays(z).find((l) => l.lay.system === '/data/c/A')!.lay;
      expect(a.step).toBe('letters');
      return Array.from(a.fill.subarray(0, a.n));
    };
    expect(await fills(moleculesSpec({ color }))).toEqual(await fills(moleculesSpec()));
  });

  it('a unit’s atoms read their unit’s first component, not the molecule’s (a unit over its second residue)', async () => {
    const b = createSeqSpec();
    const m = b.data.molecule('m', {
      residues: ['NAG', 'BMA'],
      atoms: { residue: [0, 0, 1, 1], name: ['C1', 'O1', 'C1', 'O2'], element: ['C', 'O', 'C', 'O'] },
    });
    b.data.complex('c', {}, (c) => c.unit('u', { of: m, range: [1, 2] }));
    b.views.view('v', (v) => {
      const s = v.section('s', { of: '/data/c' as never });
      s.block('u', { form: 'atoms' });
      s.track('t').layer('l', {
        type: 'letters',
        color: theme.rules([theme.rule([theme.eq('component', 'BMA')], '#d62728')]),
      });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size });
    const [atoms] = lettersLays(f);
    expect(atoms.lay.n).toBeGreaterThan(0);
    expect(atoms.lay.fill[0]).toBe(0xd62728ff);
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

  it('the links’ ends plans are built where a ranges layer over them asks, on its systems only', async () => {
    const arcs = await settle(moleculesSpec({}));
    const ends = arcs.model.annotations.get('/data/c/links')!.ends!;
    expect([...ends.systems].sort()).toEqual(['/data/c/A', '/data/c/B', '/data/c/g', '/data/c/l']);
    // arcs only: none built
    expect([...ends.built()]).toEqual([]);
    // a ranges layer over them on a section of the chains only: the chains' two
    const sites = await settle(moleculesSpec({ copies: 'each', sites: true, classes: ['polymer'] }));
    const built = [...sites.model.annotations.get('/data/c/links')!.ends!.built()].map((p) => p.system).sort();
    expect(built).toEqual(['/data/c/A', '/data/c/B']);
  });

  it('in text, a block’s label inline on its first line, taking part in its breaking; its other lines none (EC-95)', async () => {
    const s = await settle(moleculesSpec({ layout: 'text', copies: 'each' }));
    const f = await s.frame({ size: { width: 520, height: 600 } });
    const p = f.layout.sections[0].placement;
    const a = (i: number) => p.strips.axes[i];
    // A's first line, then its continuation lines on the same system, then B's first line
    const lines = [];
    for (let i = 0; i < p.strips.n; i++) if (!p.strips.chip[i]) lines.push(i);
    const onA = lines.filter((i) => a(i).system === '/data/c/A');
    expect(onA.length).toBeGreaterThan(1);
    const first = a(onA[0]);
    const next = a(onA[1]);
    // the first line starts after its label: fewer units, its letters shifted right by the label
    expect(first.x0).toBeGreaterThan(next.x0);
    expect(first.u1 - first.u0).toBeLessThan(next.u1 - next.u0);
    const labels = texts(f, 'labels');
    expect(labels.filter((t) => t === 'A')).toHaveLength(1);
    expect(labels.filter((t) => t === 'B')).toHaveLength(1);
  });

  it('in text with numbers, a first line’s start number stays in the gutter, clear of its inline label (EC-95)', async () => {
    const s = await settle(moleculesSpec({ layout: 'text', copies: 'each', numbers: 'start' }));
    const f = await s.frame({ size: { width: 520, height: 600 } });
    // the labels pass's texts with their x, anchor and measured width
    const marks = recordDrawLists(f.lists, ['labels'])
      .split('\n')
      .filter((l) => l.startsWith('text '))
      .map((l) => {
        const m = /^text "([^"]*)" ([\d.-]+) ([\d.-]+) (\w+) .* width ([\d.]+)/.exec(l)!;
        const x = Number(m[2]);
        const w = Number(m[5]);
        const [x0, x1] = m[4] === 'end' ? [x - w, x] : m[4] === 'middle' ? [x - w / 2, x + w / 2] : [x, x + w];
        return { text: m[1], y: Number(m[3]), x0, x1 };
      });
    const label = marks.find((m) => m.text === 'A')!;
    const number = marks.find((m) => m.text === '1' && Math.abs(m.y - label.y) < 1)!;
    expect(label).toBeDefined();
    expect(number).toBeDefined();
    expect(number.x1).toBeLessThanOrEqual(label.x0);
    const svg = drawListsToSvg(f.lists);
    expect(svg).toContain('>A<');
    expect(svg).toContain('>1<');
  });
});
