/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync, readdirSync } from 'node:fs';
import type { Node } from '@molstar/core/repr/tree';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore, SeqRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { openModel } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import { nodeFixtureAssets } from '../../node-assets.js';
import { ExampleViews } from '../../views.js';
import { readCoreSvg } from './core-svg.js';
import { type Drawn, drawSection, type Mark } from './draw.js';

// The second consumer of the resolved view (EC-183; rendering.md "A second consumer, in tests only"): a plain-SVG drawing
// of letters, ranges, values and variants from the resolved view and the store alone (`./draw.ts`, `./read.ts`: public
// entries of `@molstar/seqx-core` only). The field test: the same spec drawn by core and by the consumer agrees on what is
// where, read back from core's SVG as positions (its plot maps a window of positions to its width): which letter at which
// position in which colour, which range from where to where with which label, which bar at which height, which head at
// which site. Layouts differ (px, heights, lanes); positions, text and colours by value do not.
//
// What is independent: positions, text, heights (shares of the track), the dedupe of a range over a tile edge, the
// complement on a reversed axis. What is shared: the colours, through the theme compile (`compileChannel`,
// `builtInSchemes`, the token table) that core's channels use too, so a colour check shows the resolved view carries
// the forms and the consumer reads the data they run over, never that the theme compile is right (its own tests are).
// Data access (`layerSource`) and the token table are core's, public through `render`.

const Width = 720;

interface Both {
  readonly mine: Drawn;
  /** Core's drawing read back in positions. */
  readonly core: ReturnType<typeof readCoreSvg>;
  /** The continuous position a core x is at. */
  at(x: number): number;
  /** A core track by the start of its accessible name. */
  track(name: string): { x: number; y: number; width: number; height: number };
}

/** Draws a spec both ways over the same window of positions of its first section. */
async function both(spec: Node, window: readonly [number, number], reversed = false): Promise<Both> {
  const ctx = createSeqContext({
    registry: [{ kind: 'format-runtime', runtime: SeqRuntime }],
    ports: { assets: () => nodeFixtureAssets() },
  });
  const store = new DataStore();
  const opened = await openModel(spec, store, {
    runtimes: ctx.runtimes,
    config: ctx.resolveConfig,
    assets: ctx.assets,
  });
  const section = opened.view!.sections[0];
  const mine = await drawSection(opened.model, opened.view!, store, { width: Width, window });
  const r = await ctx.renderSvg(spec, {
    size: { width: Width, height: 320 },
    viewport: { positions: { [section.address]: [window[0], window[1]] } },
  });
  ctx.retire();
  const core = readCoreSvg(r.svg);
  const first = [...core.tracks.values()][0];
  const [ws, we] = window;
  const at = (x: number) =>
    reversed ? we - ((x - first.x) / first.width) * (we - ws) : ws + ((x - first.x) / first.width) * (we - ws);
  return {
    mine,
    core,
    at,
    track(name) {
      const t = [...core.tracks].find(([k]) => k.startsWith(name));
      if (!t) throw new Error(`no track ${name}`);
      return t[1];
    },
  };
}

const of = (b: Both, kind: Mark['kind']) => b.mine.marks.filter((m) => m.kind === kind);
const inTrack = <T extends { y: number }>(rows: readonly T[], t: { y: number; height: number }) =>
  rows.filter((r) => r.y >= t.y - 0.5 && r.y < t.y + t.height);

/** A mark's span in positions, clipped to the window: what core's shape snaps to pixels of. */
const span = (m: Mark, [ws, we]: readonly [number, number]) => [Math.max(m.position, ws), Math.min(m.end, we)] as const;

const p53 = (id: string) => ExampleViews.find((v) => v.dataset === 'p53' && v.id === id)!.build();

describe('EC-183 the second consumer draws what core draws', () => {
  it('EC-183 letters: each letter at its position, its cell in the colour core fills and its ink', async () => {
    const win = [60, 85] as const;
    const b = await both(p53('tracks'), win);
    const letters = of(b, 'letter');
    expect(letters.map((m) => m.position)).toEqual(Array.from({ length: 25 }, (_, i) => 60 + i));
    expect(b.mine.pending).toBe(0);
    // core's runs: a letter at its cell's centre, read back as the position it stands at
    expect(b.core.letters).toHaveLength(letters.length);
    const coreAt = new Map(b.core.letters.map((l) => [Math.round(b.at(l.x) - 0.5), l.c]));
    expect(coreAt.size).toBe(letters.length);
    for (const m of letters) expect(coreAt.get(m.position), `position ${m.position}`).toBe(m.text);
    // core's cells: the fill by position
    const cells = inTrack(b.core.rects, b.track('p53'));
    expect(cells).toHaveLength(letters.length);
    for (const c of cells) {
      const m = letters.find((l) => l.position === Math.floor(b.at((c.x0 + c.x1) / 2)))!;
      expect(c.fill, `cell of ${m.position}`).toBe(m.fill);
    }
    // the consumer's own drawing says the same letters, as text
    expect([...b.mine.svg.matchAll(/text-anchor="middle"[^>]*>(.)<\/text>/g)].map((x) => x[1]).join('')).toBe(
      letters.map((m) => m.text).join(''),
    );
  });

  it('EC-183 ranges: each domain from where to where, its label, its colour by its category', async () => {
    const win = [80, 130] as const;
    const b = await both(p53('tracks'), win);
    const ranges = of(b, 'range');
    // PRR (63–94) ends in the window and DBD (93–292) starts in it; TAD is out
    expect(ranges.map((m) => m.text)).toEqual(['Proline-rich', 'DNA-binding']);
    const shapes = inTrack(b.core.paths, b.track('Domains')).sort((p, q) => p.x0 - q.x0);
    expect(shapes).toHaveLength(ranges.length);
    ranges.forEach((m, i) => {
      const [s, e] = span(m, win);
      expect(b.at(shapes[i].x0), m.text).toBeCloseTo(s, 0);
      expect(b.at(shapes[i].x1), m.text).toBeCloseTo(e, 0);
      expect(shapes[i].fill, m.text).toBe(m.fill);
      // core names each item as `label, system from–to` (from is 1-based): the label and the range
      expect(
        b.core.items.find((x) => x.startsWith(`${m.text},`)),
        m.text,
      ).toBe(`${m.text}, p53 ${m.position + 1}–${m.end}`);
    });
    // the labels core writes inside its shapes are the consumer's
    for (const t of b.core.texts.filter((x) => x.size === 10)) expect(ranges.map((m) => m.text)).toContain(t.text);
  });

  it('EC-183 values: each bar at its position, in the value colour, as high as its value in the domain in view', async () => {
    const win = [60, 85] as const;
    const b = await both(p53('tracks'), win);
    const bars = of(b, 'value');
    expect(bars).toHaveLength(25);
    const t = b.track('Conservation');
    // core's bars: rounded shapes (the track's hairline at the baseline is a rect, no bar)
    const coreBars = inTrack(b.core.paths, t);
    expect(coreBars).toHaveLength(bars.length);
    for (const c of coreBars) {
      const m = bars.find((x) => x.position === Math.floor(b.at((c.x0 + c.x1) / 2)))!;
      expect(c.fill).toBe(m.fill);
      // the same share of the track less the headroom (a tenth): core snaps heights to the pixel grid
      expect(c.h / t.height / 0.9, `bar of ${m.position}`).toBeCloseTo(m.share!, 2);
      // core's bars are the series color at 60% (the SVG's `fill-opacity`)
      expect(c.alpha, `bar of ${m.position}`).toBeCloseTo(0.6, 2);
    }
    // the two missing values (100, 101) draw no bar
    const missing = await both(p53('tracks'), [95, 105]);
    const expected = [95, 96, 97, 98, 99, 102, 103, 104];
    expect(of(missing, 'value').map((m) => m.position)).toEqual(expected);
    // core leaves them a gap too: no mark in the missing colour, no stub at the floor
    const marks = inTrack(missing.core.paths, missing.track('Conservation')).map((c) =>
      Math.floor(missing.at((c.x0 + c.x1) / 2)),
    );
    expect(marks).toEqual(expected);
  });

  it('EC-183 variants: each head at its site, in the colour of its consequence', async () => {
    const win = [170, 200] as const;
    const b = await both(p53('variants'), win);
    const heads = of(b, 'variant');
    expect(heads.map((m) => m.text)).toEqual(['R175H', 'R196*']);
    // the heads are the shapes of the variants track (its only layer; core's ruler ticks are rects in another track)
    const coreHeads = inTrack(b.core.paths, b.track('Variants')).sort((p, q) => p.x0 - q.x0);
    expect(coreHeads).toHaveLength(heads.length);
    heads.forEach((m, i) => {
      expect(b.at((coreHeads[i].x0 + coreHeads[i].x1) / 2)).toBeCloseTo(m.position + 0.5, 0);
      expect(coreHeads[i].fill, m.text).toBe(m.fill);
    });
    // the second section follows the first through a map: a layer behind a path is left to a consumer that builds them
    expect(b.mine.skipped).toContain('/views/main/protein/ruler/ruler');
  });
});

// A system long enough to span two tiles of the store (1,024 positions), shown reversed and forward across the edge.
function long(direction?: 'reversed'): Node {
  const b = createSeqSpec({ title: 'long' });
  const n = 1500;
  const seq = b.data.sequence('s', {
    alphabet: 'dna',
    letters: Array.from({ length: n }, (_, i) => 'ACGT'[(i * 7 + (i >> 3)) & 3]).join(''),
  });
  const values = b.data.values('v', {
    target: seq,
    values: Float32Array.from({ length: n }, (_, i) => ((i * 37) % 100) / 100),
  });
  const ranges = b.data.ranges('r', {
    target: seq,
    rows: {
      start: [1005, 1020, 1030],
      end: [1015, 1040, 1100],
      id: ['a', 'b', 'c'],
      label: ['one', 'two', 'three'],
      category: ['x', 'y', 'x'],
    },
  });
  b.views.view('main', (v) => {
    const s = v.section('s', { of: seq, ...(direction ? { direction } : {}) });
    s.track('letters', { label: 'DNA' }).layer('letters', { type: 'letters' });
    s.track('ranges', { label: 'Ranges' }).layer('ranges', { type: 'ranges', data: ranges });
    s.track('values', { label: 'Values' }).layer('values', { type: 'values', data: values });
  });
  return b.build();
}

describe('EC-183 across a tile edge and on a reversed axis', () => {
  for (const reversed of [false, true]) {
    it(`EC-183 a window over the edge between two tiles, ${reversed ? 'reversed' : 'forward'}: the same letters, ranges and bars`, async () => {
      const win = [1010, 1040] as const;
      const b = await both(long(reversed ? 'reversed' : undefined), win, reversed);
      expect(b.mine.pending).toBe(0);
      const letters = of(b, 'letter');
      expect(letters).toHaveLength(30);
      // the order on the screen is the axis's: ascending forward, descending reversed
      const xs = letters.map((m) => m.x);
      expect(xs.every((x, i) => i === 0 || (reversed ? x < xs[i - 1] : x > xs[i - 1]))).toBe(true);
      const coreAt = new Map(b.core.letters.map((l) => [Math.round(b.at(l.x) - 0.5), l.c]));
      for (const m of letters) expect(coreAt.get(m.position), `position ${m.position}`).toBe(m.text);
      const bars = of(b, 'value');
      const t = b.track('Values');
      const coreBars = inTrack(b.core.paths, t);
      expect(coreBars).toHaveLength(bars.length);
      for (const c of coreBars) {
        const m = bars.find((x) => x.position === Math.floor(b.at((c.x0 + c.x1) / 2)))!;
        expect(c.fill).toBe(m.fill);
        expect(c.h / t.height / 0.9, `bar of ${m.position}`).toBeCloseTo(m.share!, 2);
        expect(c.alpha, `bar of ${m.position}`).toBeCloseTo(0.6, 2);
      }
      const ranges = of(b, 'range');
      expect(ranges.map((m) => m.text)).toEqual(['one', 'two', 'three']);
      for (const m of ranges)
        expect(
          b.core.items.some((x) => x.startsWith(`${m.text},`)),
          m.text,
        ).toBe(true);
      // and where they are: core's shapes, by the positions they span (a range over the tile edge once, clipped to the window)
      const shapes = inTrack(b.core.paths, b.track('Ranges'))
        .map((p) => [Math.min(b.at(p.x0), b.at(p.x1)), Math.max(b.at(p.x0), b.at(p.x1)), p.fill] as const)
        .sort((p, q) => p[0] - q[0] || p[1] - q[1]);
      expect(shapes).toHaveLength(ranges.length);
      const mine = ranges.map((m) => [...span(m, win), m.fill] as const).sort((p, q) => p[0] - q[0] || p[1] - q[1]);
      mine.forEach(([s, e, fill], i) => {
        expect(shapes[i][0]).toBeCloseTo(s, 0);
        expect(shapes[i][1]).toBeCloseTo(e, 0);
        expect(shapes[i][2]).toBe(fill);
      });
    });
  }
});

describe('EC-183 the consumer needs nothing private', () => {
  it('EC-183 imports only the package’s listed entries, Mol*’s tree and itself', () => {
    const pkg = JSON.parse(readFileSync(new URL('../../../../core/package.json', import.meta.url), 'utf8')) as {
      exports: Record<string, unknown>;
    };
    // the entries a package lists by name: not its `./*` wildcard (a path into `src/`), not a hidden folder (null)
    const entries = Object.entries(pkg.exports)
      .filter(([k, v]) => k !== './*' && v !== null)
      .map(([k]) => `@molstar/seqx-core${k.slice(1)}`);
    const dir = new URL('./', import.meta.url);
    const files = readdirSync(dir).filter((f) => f.endsWith('.ts') && f !== 'consumer.test.ts' && f !== 'core-svg.ts');
    expect(files.sort()).toEqual(['draw.ts', 'read.ts']);
    const specifiers: string[] = [];
    for (const f of files) {
      const text = readFileSync(new URL(f, dir), 'utf8');
      // static imports (type ones too), re-exports and dynamic `import('…')`
      for (const m of text.matchAll(/from '([^']+)'/g)) specifiers.push(m[1]);
      for (const m of text.matchAll(/\bimport\(\s*'([^']+)'/g)) specifiers.push(m[1]);
      // `render` is data access and the token table alone: no layout, draw list or backend
      for (const m of text.matchAll(/(?:import|export)[^;]*?'@molstar\/seqx-core\/render'/gs)) {
        const names = /\{([^}]*)\}/
          .exec(m[0])?.[1]
          .split(',')
          .map((n) => n.trim().replace(/^type\s+/, ''))
          .filter(Boolean);
        expect(names, m[0]).toBeDefined();
        for (const n of names!) expect(['LightTokens', 'layerSource', 'LayerSource'], n).toContain(n);
      }
      expect(text).not.toMatch(/import\(\s*[^'\s)]/);
    }
    for (const s of specifiers) {
      if (s.startsWith('./')) continue;
      if (s.startsWith('@molstar/seqx-core')) expect(entries, s).toContain(s);
      else expect(s, 'only Mol*’s tree besides').toBe('@molstar/core/repr/tree');
    }
    expect(specifiers.filter((s) => s.startsWith('@molstar/seqx-core')).length).toBeGreaterThan(4);
  });
});
