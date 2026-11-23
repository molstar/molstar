/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Chips under motion and at their edges: two chips stacked draw each pill once, nothing painted over it; a pan through
// the row's inverse (continuous in the window, what is under the pointer moving with it); a zoom keeping its anchor's
// px as a chip leaves the view; a chip's label never lost (no pill: in the label row); a chip at its own first unit;
// chips in text flowing like words; a pill on a ring's letters; the chips' floor (half the width at most, their labels
// elided); the pill's hairline (dark), ink (4.5:1) and SVG rect; a label that does not fit its block elided or dropped.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { contrastRgba } from '@molstar/seqx-core/theme/index';
import { describe, expect, it } from 'vitest';
import {
  type BlockPlacement,
  BlockMetrics,
  bundledTextMetrics,
  DarkTokens,
  DefaultRenderConfig,
  drawListsToSvg,
  type Frame,
  LightTokens,
  PlacePoint,
  type RingPlacement,
  recordDrawLists,
  TextLayout,
  type TextPlacement,
  TracksLayout,
} from '../index.js';
import { seeded } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 800, height: 600 };
const text = bundledTextMetrics();
const chipWidth = (label: string) => text.measure('12px sans', label) + 2 * BlockMetrics.chipPad;

interface ChipsOptions {
  readonly arrange?: 'row' | 'stack';
  readonly layout?: 'rings' | 'text';
  /** The ligands, in order, after the chains (or `first`: before them). */
  readonly ligands?: readonly string[];
  readonly first?: boolean;
  /** The ligands between the chains. */
  readonly middle?: boolean;
  /** Draw no letters track (a chip with no pill). */
  readonly noLetters?: boolean;
  readonly cdk2?: number;
  readonly cdk2Label?: string;
}

/** CDK2 and cyclin A with ligands as chips, a features and a letters track. */
function chipsSpec(o: ChipsOptions = {}) {
  const b = createSeqSpec({ title: 'chips' });
  const r = seeded(7);
  const protein = (n: number) => {
    let s = '';
    for (let i = 0; i < n; i++) s += 'ACDEFGHIKLMNPQRSTVWY'[Math.floor(r() * 20)];
    return s;
  };
  const a = b.data.sequence('cdk2', { alphabet: 'protein', letters: protein(o.cdk2 ?? 120) });
  const c = b.data.sequence('cyclin', { alphabet: 'protein', letters: protein(80) });
  const names = o.ligands ?? ['ATP', 'MG'];
  // the first two by their names (ATP, MG), the rest by their index
  const ligands = names.map((l, k) =>
    b.data.sequence(k < 2 && /^[A-Z]+$/.test(l) ? l.toLowerCase() : `lig${k}`, {
      alphabet: 'protein',
      letters: 'X',
      name: l,
    }),
  );
  const fa = b.data.ranges('cdk2-features', {
    target: a,
    rows: { start: [40], end: [52], id: ['psta'], label: ['PSTAIRE'] },
  });
  b.views.view('v', (v) => {
    const s = v.section('complex', {
      ...(o.arrange ? { arrange: o.arrange } : {}),
      ...(o.layout ? { layout: o.layout } : {}),
    });
    const lig = () => ligands.forEach((l, k) => s.block(`l${k}`, { of: l, label: (o.ligands ?? ['ATP', 'MG'])[k] }));
    if (o.first) lig();
    s.block('a', { of: a, label: o.cdk2Label ?? 'CDK2' });
    if (o.middle) lig();
    s.block('b', { of: c, label: 'cyclin A' });
    if (!o.first && !o.middle) lig();
    s.track('features', {}, (t) => t.layer('cdk2', { type: 'ranges', data: fa }));
    if (!o.noLetters) s.track('letters').layer('letters', { type: 'letters' });
  });
  return b.build();
}

const rects = (f: Frame, pass: 'data' | 'labels' = 'data') =>
  recordDrawLists(f.lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('rect '));
const texts = (f: Frame, pass: 'data' | 'labels') =>
  recordDrawLists(f.lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('text '));
const placementOf = (f: Frame) => f.layout.sections[0].placement as BlockPlacement;
const chipsOf = (p: { strips: BlockPlacement['strips'] }) =>
  p.strips.axes.map((_, i) => i).filter((i) => p.strips.chip[i] === 1);

describe('two chips stacked', () => {
  it('each pill drawn once, no later rect over it; the track labels once a row', async () => {
    const s = await settle(chipsSpec({ arrange: 'stack' }));
    const f = await s.frame({ size: Size });
    const p = placementOf(f);
    const chips = chipsOf(p);
    expect(chips.length).toBe(2);
    expect(p.strips.top[chips[0]]).toBe(p.strips.top[chips[1]]);
    const lines = rects(f);
    // the chips' pills: a squircle a pixel inside each chip
    const pills = lines
      .map((l, k) => ({ k, r: l.split(' ').slice(1, 5).map(Number), pill: l.includes('squircle') }))
      .filter((x) => x.pill && chips.some((i) => Math.abs(x.r[0] - p.strips.axes[i].x0 - 1) < 0.5));
    expect(pills.length).toBe(2);
    for (const { k, r } of pills) {
      const [x, y, w, h] = r;
      // no rect after the pill covers its box
      for (let j = k + 1; j < lines.length; j++) {
        const [x2, y2, w2, h2] = lines[j].split(' ').slice(1, 5).map(Number);
        const covers = x2 <= x && y2 <= y && x2 + w2 >= x + w && y2 + h2 >= y + h;
        expect(covers).toBe(false);
      }
    }
    // each row's track labels once: 'letters' three times (two chains' rows and the chips')
    const labels = texts(f, 'labels').filter((l) => l.startsWith('text "letters"'));
    expect(labels.length).toBe(3);
  });
});

describe('a row under a pan', () => {
  it('a 0.01-unit pan moves no strip by more than a pixel, a chip crossing either edge', async () => {
    const s = await settle(chipsSpec());
    // ATP at 200, MG at 201 (total 202): windows of 30 units swept across the chips' edges
    let prev: Map<number, [number, number]> | undefined;
    for (let start = 168; start <= 172.2; start += 0.01) {
      const f = await s.frame({ size: Size, viewport: { start, end: start + 30 } });
      const p = placementOf(f);
      const now = new Map<number, [number, number]>();
      for (let i = 0; i < p.strips.n; i++)
        now.set(p.strips.block[i], [p.strips.axes[i].x0, p.strips.axes[i].x0 + p.strips.axes[i].width]);
      if (prev)
        for (const [b, [x0, x1]] of now) {
          const q = prev.get(b);
          if (!q) {
            expect(x1 - x0).toBeLessThanOrEqual(1);
            continue;
          }
          expect(Math.abs(x0 - q[0])).toBeLessThanOrEqual(1);
          expect(Math.abs(x1 - q[1])).toBeLessThanOrEqual(1);
        }
      prev = now;
    }
  });

  it('what is under the pointer moves with it, chips in view', async () => {
    // ATP at 120, MG at 121, between the chains: in view all along
    const s = await settle(chipsSpec({ middle: true }));
    const v0 = { start: 95, end: 150 };
    const f0 = await s.frame({ size: Size, viewport: v0 });
    const p0 = placementOf(f0);
    const axis = f0.layout.sections[0].axis!;
    const c = { total: 202, width: axis.width, height: 600, cyclic: false, maxScale: 40, placement: p0 };
    // a pointer on cyclin A, then pans of 30 and 8 px either way
    const y = p0.strips.bandY[1] + 2;
    const px = p0.strips.axes[3].x0 + 60;
    const u = p0.offsetAt(px, y);
    for (const dx of [-30, -8, 8, 30]) {
      const v = TracksLayout.nav.pan(c, v0, dx, 0);
      const f = await s.frame({ size: Size, viewport: v });
      const p = placementOf(f);
      const pt = new PlacePoint();
      expect(p.place(u, pt)).toBe(true);
      expect(pt.x).toBeCloseTo(px - dx, 1);
      // the chips keep their slots
      for (const i of chipsOf(p)) expect(p.strips.axes[i].width).toBeCloseTo(chipWidth(p.labels[i]), 6);
    }
  });
});

describe('a row panned about the pointer', () => {
  it('the content a pan brings under the pointer is the unit that was dx px right of it, chips and gaps between', async () => {
    // the chips (ATP, MG) at the end enter the view as the window pans right
    // the chips (ATP, MG) at the end enter the view as the window pans right
    const s = await settle(chipsSpec());
    const v0 = { start: 130, end: 185 };
    const f0 = await s.frame({ size: Size, viewport: v0 });
    const p0 = placementOf(f0);
    const axis = f0.layout.sections[0].axis!;
    const c = { total: 202, width: axis.width, height: 600, cyclic: false, maxScale: 40, placement: p0 };
    const y = p0.strips.bandY[1] + 2;
    // a pointer on CDK2, the unit dx px right of it (on a block, not in a gap) comes to it
    // (the placement is reused by the next frame: the units first)
    const moves = [
      [100, 150],
      [100, 200],
      [300, 100],
    ] as const;
    const units = moves.map(([x, dx]) => p0.offsetAt(axis.x0 + x + dx, y));
    for (const [i, [x, dx]] of moves.entries()) {
      const u = units[i];
      const v = TracksLayout.nav.pan(c, v0, dx, 0, x);
      const f = await s.frame({ size: Size, viewport: v });
      const p = placementOf(f);
      const pt = new PlacePoint();
      expect(p.place(u, pt)).toBe(true);
      expect(pt.x).toBeCloseTo(axis.x0 + x, 0);
    }
  });
});

describe('a zoom about a point in a row with a chip first', () => {
  it('the point keeps its px as the chip leaves the view', async () => {
    const s = await settle(chipsSpec({ first: true, ligands: ['ATP'] }));
    const f0 = await s.frame({ size: Size });
    const p0 = placementOf(f0);
    expect(chipsOf(p0)).toEqual([0]);
    const axis = f0.layout.sections[0].axis!;
    const y = p0.strips.bandY[1] + 2;
    let v = { start: 0, end: 201 };
    let p = p0;
    // a pointer on CDK2, a third of its way in
    const px = p0.strips.axes[1].x0 + p0.strips.axes[1].width / 3;
    const at = p0.offsetAt(px, y);
    for (let k = 0; k < 4; k++) {
      const c = { total: 201, width: axis.width, height: 600, cyclic: false, maxScale: 40, placement: p };
      v = TracksLayout.nav.zoom(c, v, 1.6, at) as typeof v;
      const f = await s.frame({ size: Size, viewport: v });
      p = placementOf(f);
      const pt = new PlacePoint();
      expect(p.place(at, pt)).toBe(true);
      // within the device grid's half pixel (the blocks' edges are snapped to it)
      expect(Math.abs(pt.x - px)).toBeLessThanOrEqual(0.5);
    }
    // the chip left the view
    expect(chipsOf(p)).toEqual([]);
  });
});

describe('a chip’s label never lost', () => {
  it('with no letters track, and before its letter arrives, the label is in the label row', async () => {
    const s = await settle(chipsSpec({ noLetters: true }));
    const f = await s.frame({ size: Size });
    const labels = texts(f, 'labels');
    expect(labels.some((l) => l.startsWith('text "ATP"'))).toBe(true);
    expect(labels.some((l) => l.startsWith('text "MG"'))).toBe(true);
    // a frame before any tile is held: no pill yet, the labels said above
    const fresh = await settle(chipsSpec());
    const first = fresh.renderer.frame({ size: Size });
    const said = recordDrawLists(first.lists, ['labels', 'data'])
      .split('\n')
      .filter((l) => l.startsWith('text "ATP"'));
    expect(said.length).toBe(1);
  });
});

describe('a chip at its own first unit', () => {
  it('place(o) is the chip, not the end of the block before it', async () => {
    const s = await settle(chipsSpec());
    const f = await s.frame({ size: Size });
    const p = placementOf(f);
    const pt = new PlacePoint();
    expect(p.place(200, pt)).toBe(true);
    expect(p.strips.chip[pt.strip]).toBe(1);
    expect(p.labels[pt.strip]).toBe('ATP');
    expect(p.place(201, pt) && p.labels[pt.strip]).toBe('MG');
    expect(p.place(120, pt) && p.labels[pt.strip]).toBe('cyclin A');
  });
});

describe('chips in text', () => {
  it('each chain starts a line on its own system, the chips flow like words on one line, their label’s width', async () => {
    const s = await settle(chipsSpec({ layout: 'text' }));
    const f = await s.frame({ size: Size });
    const sec = f.layout.sections[0];
    expect(sec.layout).toBe(TextLayout);
    const p = sec.placement as TextPlacement;
    const systems = p.strips.axes.map((a) => a.system);
    expect(systems).toContain('/data/cdk2');
    expect(systems).toContain('/data/cyclin');
    const chips = chipsOf(p);
    expect(chips.map((i) => p.strips.axes[i].system)).toEqual(['/data/atp', '/data/mg']);
    // on one line, side by side, each its label's width
    expect(p.strips.top[chips[0]]).toBe(p.strips.top[chips[1]]);
    expect(p.strips.axes[chips[1]].x0).toBeGreaterThan(p.strips.axes[chips[0]].x0 + p.strips.axes[chips[0]].width);
    expect(p.strips.axes[chips[0]].width).toBeCloseTo(chipWidth('ATP'), 6);
    expect(p.strips.axes[chips[1]].width).toBeCloseTo(chipWidth('MG'), 6);
    // a chain's line starts at its block's first unit: the cyclin's first line at section unit 120
    const cyc = systems.indexOf('/data/cyclin');
    expect(p.strips.first[cyc]).toBe(120);
    // placement is the inverse of hit-testing across the lines and the chips
    const pt = new PlacePoint();
    for (const u of [0, 50, 119, 120, 199, 200, 201]) {
      expect(p.place(u + 0.5, pt)).toBe(true);
      expect(p.at([pt.x, pt.y], pt) && Math.floor(pt.offset)).toBe(u);
    }
    // the pills drawn
    const labels = texts(f, 'data');
    expect(labels.some((l) => l.startsWith('text "ATP"'))).toBe(true);
    expect(labels.some((l) => l.startsWith('text "MG"'))).toBe(true);
  });
});

describe('a pill on a ring', () => {
  it('the chip’s letters ring draws its label in a pill, not a letter cell', async () => {
    const s = await settle(chipsSpec({ layout: 'rings', ligands: ['ATP'] }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as RingPlacement;
    expect(chipsOf(p).length).toBe(1);
    const data = recordDrawLists(f.lists, ['data']);
    expect(data).toMatch(/text "ATP"/);
    // once: the pill draws it, no outside label with a leader beside it
    const all = recordDrawLists(f.lists, ['data', 'labels']);
    expect(all.match(/text "ATP"/g)?.length).toBe(1);
  });
});

describe('the chips’ floor', () => {
  it('many chips take at most half the width, shrunk to no less than the least, their labels whole, cut keeping `minChars` letters, or none', async () => {
    const ligands = Array.from({ length: 12 }, (_, k) => `Ligand number ${k}`);
    const s = await settle(chipsSpec({ ligands }));
    const f = await s.frame({ size: Size });
    const p = placementOf(f);
    const axis = f.layout.sections[0].axis!;
    const chips = chipsOf(p);
    expect(chips.length).toBe(12);
    let sum = 0;
    for (const i of chips) {
      sum += p.strips.axes[i].width;
      expect(p.strips.axes[i].width).toBeGreaterThanOrEqual(BlockMetrics.chipMin - 1e-9);
    }
    expect(sum).toBeLessThanOrEqual(axis.width / 2 + 1e-6);
    // the polymers keep the rest
    expect(p.strips.axes[0].width + p.strips.axes[1].width).toBeGreaterThan(axis.width / 2 - 0.25 * axis.width);
    // each pill's label fits its slot, never past it; a cut keeps at least `labels.minChars` letters (_spec/README.md
    // "Render", Section layouts)
    for (const l of texts(f, 'data')) {
      const w = Number(/ width ([\d.]+)/.exec(l)![1]);
      const max = Number(/ max ([\d.]+)/.exec(l)![1]);
      expect(w).toBeLessThanOrEqual(max);
      const t = /^text "([^"]*)"/.exec(l)![1];
      if (t.endsWith('…')) expect(t.length - 1).toBeGreaterThanOrEqual(DefaultRenderConfig.labels.minChars);
    }
  });

  // a pill's label follows the rule of a label inside a mark (labelsInside's): whole, or cut keeping `minChars` letters,
  // or none (the pill alone, its label in the label row); never a cut of fewer, at any width, canvas and SVG alike
  for (const width of [560, 1400])
    it(`at ${width} px, seven chips: every pill's text whole or kept to at least \`minChars\` letters, the SVG the same`, async () => {
      const ligands = ['ATP', '2AN', 'PYZ', '4SP', 'Staurosporine', 'LZ1 analogue', 'B49'];
      const s = await settle(chipsSpec({ ligands }));
      const f = await s.frame({ size: { width, height: Size.height } });
      const pills = texts(f, 'data').map((l) => /^text "([^"]*)"/.exec(l)![1]);
      for (const t of pills) {
        if (t.endsWith('…')) expect(t.length - 1, t).toBeGreaterThanOrEqual(DefaultRenderConfig.labels.minChars);
        else expect(ligands, t).toContain(t);
      }
      // the short codes are whole or not drawn: never "AT…" or "2…"
      expect(pills.filter((t) => /^.{1,3}…$/.test(t))).toEqual([]);
      const svg = drawListsToSvg(f.lists);
      for (const t of pills) expect(svg).toContain(`>${t}<`);
      expect(svg).not.toMatch(/>[^<]{1,3}…</);
    });
});

describe('the pill', () => {
  for (const [name, tokens] of [
    ['light', LightTokens],
    ['dark', DarkTokens],
  ] as const)
    it(`${name}: its ink at 4.5:1 against its fill; the cells’ hairline in dark; its rect in the SVG`, async () => {
      const s = await settle(chipsSpec({ ligands: ['ATP'] }));
      s.renderer.tokens = tokens;
      const f = await s.frame({ size: Size });
      const p = placementOf(f);
      const x0 = p.strips.axes[chipsOf(p)[0]].x0;
      const pill = rects(f).find((l) => l.includes('squircle') && Math.abs(Number(l.split(' ')[1]) - x0 - 1) < 0.5)!;
      expect(pill).toBeDefined();
      // '#rrggbb' or '#rrggbbaa' as RGBA
      const rgba = (hex: string) => {
        const h = hex.replace('#', '');
        return Number.parseInt(h.length === 6 ? `${h}ff` : h, 16) >>> 0;
      };
      const fill = rgba(pill.split(' ')[5]);
      const label = texts(f, 'data').find((l) => l.startsWith('text "ATP"'))!;
      const ink = rgba(/ (#[0-9a-f]{6,8}) width/.exec(label)![1]);
      expect(contrastRgba(ink, fill)).toBeGreaterThanOrEqual(4.5);
      if (name === 'dark') expect(pill).toContain(' stroke ');
      else expect(pill).not.toContain(' stroke ');
      // the SVG draws the pill as a path or a rect with its fill, and its label
      const svg = drawListsToSvg(f.lists);
      expect(svg).toContain('>ATP<');
      expect(svg.toLowerCase()).toContain(pill.split(' ')[5].slice(0, 7).toLowerCase());
      expect(svg).toMatch(/<(rect|path)[^>]*>/);
    });
});

describe('a block label that does not fit', () => {
  it('a sliver block’s label is elided or dropped, never drawn over the next block’s', async () => {
    const s = await settle(chipsSpec({ cdk2Label: 'cyclin-dependent kinase 2 (A)' }));
    // CDK2 a sliver at the window's left
    const f = await s.frame({ size: Size, viewport: { start: 118, end: 200 } });
    const p = placementOf(f);
    expect(p.strips.axes[0].system).toBe('/data/cdk2');
    const cdk2 = p.strips.axes[0];
    const next = p.strips.axes[1].x0;
    for (const l of texts(f, 'labels')) {
      const x = Number(l.split('"')[2].trim().split(' ')[0]);
      const w = Number(/ width ([\d.]+)/.exec(l)![1]);
      if (Math.abs(x - cdk2.x0) > 0.5) continue;
      expect(x + w).toBeLessThanOrEqual(next + 0.5);
    }
    expect(p.shown[0].length).toBeLessThan(p.labels[0].length);
  });
});
