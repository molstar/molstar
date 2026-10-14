/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// N1: per type, the element description on a small resolved view (no DOM); every acceptance case materializes with
// exactly its expected issues (a golden list).

import { createSeqViewSpec, translate } from '@molstar/seq-core/builder';
import {
  clustalXRules,
  compile,
  createRegistries,
  type DocNode,
  KernelRegistry,
  MvsRegistry,
  projectionOf,
  type ResolvedView,
  register,
  resolveView,
  type SeqViewSpecDoc,
  type StreamRequest,
} from '@molstar/seq-core/format';
import { type Issue, OrderedSet, type SeqSet } from '@molstar/seq-core/model';
import { proceduralRow } from '@molstar/seq-core/stream/tiles/server';
import { describe, expect, it } from 'vitest';
import { CASES } from '../../core/_test/format/cases.js';
import { fixtureLoaders, registries } from '../../core/_test/format/load-fixtures.js';
import { InteractionsRegistry } from '../../viewer/src/interactions/index.js';
import {
  AxisRuns,
  axesOf,
  BlockRuns,
  extentOf as extentOfRuns,
  type Located,
  locate as locateAt,
} from '../src/axis.js';
import { rangesOf } from '../src/events.js';
import { type AxisMeta, type ElementDesc, extentLoci, materialize, NightingaleSupported } from '../src/index.js';
import { AXISLESS, LABELED, lanesOf } from '../src/materialize.js';
import { fetchStreamed, type StreamedValues, streamedLevel } from '../src/stream.js';
import { highlightOf as highlightOfProjection } from './highlight-of.js';

/** The block runs of an axis range text (`a:b,…`, 1-based), as the arrays events carry. */
const extentOf = (axis: AxisMeta, text: string) =>
  BlockRuns.toArray(extentOfRuns(axis, AxisRuns.parse(text, AxisRuns.create()), BlockRuns.create()));
/** A host's loci as the manager's highlight, through the view's shared projection. */
const highlightOf = (
  view: ResolvedView,
  set: SeqSet,
  axis: AxisMeta,
  loci: Parameters<typeof highlightOfProjection>[2],
) => highlightOfProjection(projectionOf(view, set), axis, loci);
/** The block and index at an axis index, or null. */
const locate = (axis: AxisMeta, at: number) => {
  const out: Located = { block: -1, index: -1 };
  return locateAt(axis, at, out) ? { ...out } : null;
};

const kernel = createRegistries();
register(kernel, [KernelRegistry, MvsRegistry]);
const none = { load: async () => Promise.reject(new Error('no files')) };

async function draw(doc: SeqViewSpecDoc, view: string, reg = kernel) {
  const { set, issues } = await compile(doc, reg, reg === kernel ? none : fixtureLoaders());
  expect(issues.filter((i) => i.severity === 'error')).toEqual([]);
  const r = resolveView(doc, set, view, reg, NightingaleSupported);
  return { set, view: r.view, resolved: r.issues, ...materialize(r.view, set) };
}
/** Every described element of a tag, depth first. */
const all = (d: ElementDesc, tag: string, out: ElementDesc[] = []): ElementDesc[] => {
  if (d.tag === tag) out.push(d);
  for (const c of d.children) all(c, tag, out);
  return out;
};
const one = (d: ElementDesc, tag: string) => all(d, tag)[0] as ElementDesc;
const data = (d: ElementDesc, tag: string) => one(d, tag).properties.data as Record<string, unknown>[];
const said = (issues: Issue[]) => issues.map((i) => `${i.code} ${i.path}`);

const SEQ = 'MKVLAAGIVGLLLAQ'; // 15
function small(build: (b: ReturnType<typeof createSeqViewSpec>, s: string) => void) {
  const b = createSeqViewSpec({ title: 'small' });
  b.data.sequence({ ref: 's', value: SEQ, label: 'S' });
  build(b, 's');
  return b.getDocument();
}

describe('per type', () => {
  it('ranges across a cyclic origin (CG1b): one feature of two fragments, each drawn, both lit on hover', async () => {
    const b = createSeqViewSpec({ title: 'cyclic' });
    const s = b.data.sequence({ ref: 'p', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna', cyclic: true });
    const a = b.data.annotation({
      ref: 'a',
      geometry: 'ranges',
      target: s,
      value: { id: ['ori', 'x'], start: [30, 10], stop: [6, 20] },
    });
    b.view({ ref: 'v', label: 'V' }).block({ sequence: s }).track({ label: 'A' }).ranges({ data: a });
    const { root, issues } = await draw(b.getDocument(), 'v');
    const [ori, x] = data(root, 'nightingale-track');
    // Its two runs, [30, 40) and [0, 6), as fragments of one feature (1-based, inclusive); the hull is Nightingale's.
    expect(ori.locations).toEqual([
      {
        fragments: [
          { start: 31, end: 40 },
          { start: 1, end: 6 },
        ],
      },
    ]);
    expect([ori.start, ori.end, x.start, x.end]).toEqual([1, 40, 11, 20]);
    // Hovered, both fragments light (in axis order), not the hull between them.
    const out = rangesOf({ feature: ori } as never, 0, AxisRuns.create());
    expect(Array.from({ length: out.count }, (_, k) => [out.start[k], out.stop[k]])).toEqual([
      [0, 6],
      [30, 40],
    ]);
    expect(said(issues)).toEqual([]);
  });

  it("ranges across a block's origin seam (NG-04): a fragment per run, not the whole axis", async () => {
    const b = createSeqViewSpec({ title: 'seam' });
    const s = b.data.sequence({ ref: 'p', value: 'ACGTACGTTA'.repeat(4), alphabet: 'dna', cyclic: true });
    const a = b.data.annotation({
      ref: 'a',
      geometry: 'ranges',
      target: s,
      value: { id: ['r'], start: [25], stop: [35] },
    });
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: s, origin: 30 } as never)
      .track({ label: 'A' })
      .ranges({ data: a });
    const [r] = data((await draw(b.getDocument(), 'v')).root, 'nightingale-track');
    // Axis 1..10 holds 30..39 and 31..40 holds 0..29: [25, 35) is axis 36..40 and 1..5.
    expect(r.locations).toEqual([
      {
        fragments: [
          { start: 1, end: 5 },
          { start: 36, end: 40 },
        ],
      },
    ]);
  });

  it('letters: the sequence, and a letter theme as a colored sequence (its losses known, not reported)', async () => {
    const doc = small((b) => b.view({ ref: 'v', label: 'V' }).track({ label: 'Res' }).letters());
    const { root, issues } = await draw(doc, 'v');
    const m = one(root, 'nightingale-manager');
    expect(m.attributes).toEqual({ 'reflected-attributes': 'display-start,display-end,highlight' });
    // `length` on each element, never the manager (it re-sets it on every element at each hover).
    // One gutter on every element; the sequence's own axis hidden even without a ruler (only a ruler draws positions).
    expect(one(root, 'nightingale-sequence').attributes).toMatchObject({
      length: '15',
      class: AXISLESS,
      'margin-left': '14',
      'margin-right': '14',
    });
    expect(one(root, 'nightingale-sequence').properties.data).toBe(SEQ);
    const c = one(root, 'nightingale-colored-sequence');
    expect(c.properties.data).toBe(SEQ);
    expect(c.attributes.scale).toMatch(/(^|,)M:\d+/);
    expect(c.attributes['color-range'].split(',').length).toBeGreaterThan(1);
    expect(said(issues)).toEqual([]);
  });

  it('ranges: an item per feature with fragments, display as a shape, labels inside on request, boundaries reported', async () => {
    const doc = small((b, s) => {
      const a = b.data.annotation({
        ref: 'a',
        geometry: 'ranges',
        target: s,
        value: { start: [2, 7, 11], stop: [5, 9, 12], id: [null, 'd', 'd'], label: ['one', null, null] },
      });
      const cut = b.data.annotation({ ref: 'cut', geometry: 'ranges', target: s, value: { boundary: [4] } });
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'A' }).ranges({
        data: a,
        display: 'brackets',
        labels: true,
        color: { theme: 'uniform', params: { color: 'red' } },
      });
      v.track({ label: 'Cut' }).ranges({ data: cut });
      v.track({ label: 'Plain' }).ranges({ data: a });
    });
    const { root, issues } = await draw(doc, 'v');
    const tracks = all(root, 'nightingale-track');
    const [a, cut, plain] = tracks.map((t) => t.properties.data as Record<string, unknown>[]);
    // Filled in the layer's color, outlined darker (features of one color that touch stay apart).
    expect(a.map((f) => [f.accession, f.start, f.end, f.shape, f.fill, f.color])).toEqual([
      ['v.ranges:0:0', 3, 5, 'bridge', '#ff0000', '#a50000'],
      ['d', 8, 12, 'bridge', '#ff0000', '#a50000'],
    ]);
    expect(a[1].locations).toEqual([
      {
        fragments: [
          { start: 8, end: 9 },
          { start: 12, end: 12 },
        ],
      },
    ]);
    expect(a[1].svs).toMatchObject({ item: 1 });
    expect(one(root, 'nightingale-track').attributes.layout).toBe('non-overlapping');
    expect(cut.map((f) => [f.start, f.end, f.shape])).toEqual([[4, 5, 'bridge']]);
    // A label is always the hover hint's; laid inside the feature (mount's layer) only with `labels`.
    expect(a.map((f) => f.svs)).toMatchObject([{ label: 'one', inside: true }, {}]);
    expect(plain.map((f) => f.svs)).toMatchObject([{ label: 'one' }, {}]);
    expect(plain[0].svs).not.toHaveProperty('inside');
    expect(tracks.map((t) => LABELED in t.attributes)).toEqual([true, false, false]);
    expect(said(issues)).toEqual(['svs.render.fallback v.ranges-2']);
    // Lanes as Nightingale places them (first fit, ends inclusive); a track grows with them, 14 px a lane.
    const f = (start: number, end: number) => ({ start, end });
    expect([lanesOf([]), lanesOf([f(1, 5), f(6, 9)]), lanesOf([f(1, 5), f(5, 9), f(2, 6), f(7, 8)])]).toEqual([
      0, 1, 3,
    ]);
    expect(tracks.map((t) => t.attributes.height)).toEqual(['24', '24', '24']);
  });

  it('points and emphasis: a feature per covered position, shapes from the glyph theme', async () => {
    const doc = small((b, s) => {
      const p = b.data.annotation({
        ref: 'p',
        geometry: 'ranges',
        target: s,
        value: { start: [1, 5], stop: [3, 6], category: ['x', 'y'] },
      });
      const glyph = { theme: 'field_map', params: { map: { x: 'diamond', y: 'star' } } };
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'P' }).points({ data: p, glyph });
      const t = v.track({ label: 'E' });
      t.letters();
      t.emphasis({ data: p });
    });
    const { root, issues, resolved } = await draw(doc, 'v');
    const [pts, emph] = all(root, 'nightingale-track').map((t) => t.properties.data as Record<string, unknown>[]);
    expect(pts.map((f) => [f.start, f.shape])).toEqual([
      [2, 'diamond'],
      [3, 'diamond'],
      [6, 'circle'],
    ]);
    expect(emph.map((f) => f.start)).toEqual([2, 3, 6]);
    expect(resolved.map((i) => i.code)).toEqual(['svs.render.fallback']); // emphasis drawn as points
    // Track e's two layers apart: said once for the section, naming it.
    expect(said(issues)).toEqual(['svs.render.glyph v.points', 'svs.render.fallback v.main']);
    expect(issues[1].message).toBe(
      'Nightingale: no layer compositing: each layer is its own element, in order (1 track(s): "v.e")',
    );
  });

  it('variants: single substitutions on the axis, the rest reported', async () => {
    const doc = small((b, s) => {
      const x = b.data.annotation({
        ref: 'x',
        geometry: 'ranges',
        target: s,
        value: { start: [0, 3, 6], stop: [1, 5, 7], alt: ['A', 'GG', '*'] },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'Var' }).variants({ data: x });
    });
    const { root, issues } = await draw(doc, 'v');
    const v = one(root, 'nightingale-variation-canvas').properties.data as {
      sequence: string;
      variants: Record<string, unknown>[];
    };
    expect(v.sequence).toBe(SEQ);
    expect(v.variants.map((x) => [x.start, x.variant])).toEqual([
      [1, 'A'],
      [7, '*'],
    ]);
    expect(v.variants[0].color).toMatch(/^#/);
    expect(said(issues)).toEqual(['svs.render.fallback v.variants']);
  });

  it('secondary: helix and strand shapes, a turn as a rectangle, DSSP letters reduced', async () => {
    const doc = small((b, s) => {
      const ss = b.data.annotation({
        ref: 'ss',
        geometry: 'ranges',
        target: s,
        value: { start: [0, 5, 10], stop: [4, 8, 12], category: ['G', 'strand', 'T'] },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'SS' }).secondary({ data: ss });
    });
    const { root, issues } = await draw(doc, 'v');
    expect(data(root, 'nightingale-track').map((f) => [f.start, f.end, f.shape])).toEqual([
      [1, 4, 'helix'],
      [6, 8, 'strand'],
      [11, 12, 'rectangle'],
    ]);
    expect(issues).toEqual([]);
  });

  it('values: one series over the axis, clamped to the domain, missing as null; bars a known loss', async () => {
    const doc = small((b, s) => {
      const d = b.data.annotation({
        ref: 'd',
        geometry: 'values',
        target: s,
        value: { index: [0, 1, 4], value: [0.5, 3, -1] },
      });
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'D' }).values({ data: d, domain: [0, 1] });
      v.track({ label: 'L' }).values({
        data: d,
        encoding: 'area',
        color: { theme: 'uniform', params: { color: '#00ff00' } },
      });
    });
    const { root, issues } = await draw(doc, 'v');
    const [bars, area] = all(root, 'nightingale-linegraph-track').map(
      (t) => (t.properties.data as Record<string, unknown>[])[0],
    );
    // Its own hover highlight off: the pointer's hover lights it, as every row.
    expect(all(root, 'nightingale-linegraph-track').every((t) => 'highlight-on-click' in t.attributes)).toBe(true);
    expect(bars.range).toEqual([0, 1]);
    expect((bars.values as { value: number | null }[]).slice(0, 5).map((x) => x.value)).toEqual([
      0.5,
      1,
      null,
      null,
      0,
    ]);
    expect([area.color, area.fill, area.range]).toEqual(['#00ff00', '#00ff00', [-1, 3]]);
    expect(said(issues)).toEqual([]);
  });

  it('color: a feature per run of one color', async () => {
    const doc = small((b, s) => {
      const d = b.data.annotation({
        ref: 'd',
        geometry: 'values',
        target: s,
        value: { index: [0, 1, 2, 5], value: [10, 10, 90, 90] },
      });
      const pal = {
        kind: 'discrete',
        mode: 'absolute',
        colors: [
          ['#ff0000', 0],
          ['#0000ff', 50],
        ],
      };
      b.view({ ref: 'v', label: 'V' })
        .track({ label: 'C' })
        .color({ data: d, color: { theme: 'field_palette', params: { field: 'value', palette: pal } } });
    });
    const { root, issues } = await draw(doc, 'v');
    expect(data(root, 'nightingale-track').map((f) => [f.start, f.end, f.color])).toEqual([
      [1, 2, '#ff0000'],
      [3, 3, '#0000ff'],
      [6, 6, '#0000ff'],
    ]);
    expect(issues).toEqual([]);
  });

  it('values are JSON numbers: decimals pass exactly, and key a palette by their text', async () => {
    const doc = small((b, s) => {
      const value = [0.15, 0.3, 0.1];
      const d = b.data.annotation({ ref: 'd', geometry: 'values', target: s, value: { value } });
      const r = b.data.annotation({ ref: 'r', geometry: 'values', target: s, value: { index: [0, 1, 2], value } });
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'D' }).values({ data: d, domain: [0, 1] });
      const pal = { kind: 'categorical', colors: { '0.15': '#ff0000', '0.3': '#0000ff' }, missing_color: '#808080' };
      v.track({ label: 'C' }).color({ data: r, color: { theme: 'field_palette', params: { palette: pal } } });
    });
    const { root, issues } = await draw(doc, 'v');
    const line = data(root, 'nightingale-linegraph-track')[0];
    expect((line.values as { value: number | null }[]).slice(0, 3).map((x) => x.value)).toEqual([0.15, 0.3, 0.1]);
    expect(data(root, 'nightingale-track').map((f) => [f.start, f.end, f.color])).toEqual([
      [1, 1, '#ff0000'],
      [2, 2, '#0000ff'],
      [3, 3, '#808080'],
    ]);
    expect(issues).toEqual([]);
  });

  it('ruler: navigation; labels other than the axis positions reported', async () => {
    const plain = small((b) => b.view({ ref: 'v', label: 'V' }).track({ label: '' }).ruler());
    const r = await draw(plain, 'v');
    expect(one(r.root, 'nightingale-navigation').attributes.height).toBe('40');
    expect(r.issues).toEqual([]);
    const sub = small((b, s) => {
      b.data.unit({ ref: 'u', sequence: s, range: [3, 12] });
      b.view({ ref: 'v', label: 'V' }).block({ unit: 'u' }).track({ label: '' }).ruler();
    });
    expect(said((await draw(sub, 'v')).issues)).toEqual(['svs.render.fallback v.ruler-2']); // the track is v.ruler
  });

  it("pairs: contacts on the axis (1-based, min-distance 0), the document's losses in one issue", async () => {
    const doc = small((b, s) => {
      const p = b.data.annotation({
        ref: 'p',
        geometry: 'pairs',
        target: s,
        value: { start: [0, 2], b_start: [10, 8], category: ['hb', 'hb'] },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'P' }).pairs({ data: p, side: 'below' });
    });
    const { root, issues } = await draw(doc, 'v');
    const links = one(root, 'nightingale-links');
    expect(links.attributes['min-distance']).toBe('0');
    expect(links.properties.contacts).toEqual([
      [1, 11, 1],
      [3, 9, 1],
    ]);
    expect(issues.map((i) => i.message)).toEqual(['Nightingale: nightingale-links does not draw categories, side']);
  });

  it("pairs over copies (NG-02): each copy's pair joined on its own block, side by side and stacked", async () => {
    const copies = (arrangement: 'side' | 'stacked') =>
      small((b, s) => {
        const p = b.data.annotation({ ref: 'p', geometry: 'pairs', target: s, value: { start: [0], b_start: [10] } });
        for (const u of ['u1', 'u2']) b.data.unit({ ref: u, sequence: s });
        const v = b.view({ ref: 'v', label: 'V', copies: 'each', arrangement } as never);
        for (const u of ['u1', 'u2']) v.block({ unit: u });
        v.track({ label: 'P' }).pairs({ data: p });
      });
    const side = await draw(copies('side'), 'v');
    expect(all(side.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([
      [
        [1, 11, 1],
        [SEQ.length + 1, SEQ.length + 11, 1],
      ],
    ]);
    const stacked = await draw(copies('stacked'), 'v');
    expect(all(stacked.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([
      [[1, 11, 1]],
      [[1, 11, 1]],
    ]);
    expect(stacked.issues.map((i) => i.message).filter((m) => m.includes('between blocks'))).toEqual([]);
    expect(all(stacked.root, 'nightingale-track')).toEqual([]);
  });

  it("pairs over copies (LA-05): a copy whose partner end is cut by its range takes no other copy's partner", async () => {
    const copies = (arrangement: 'side' | 'stacked') =>
      small((b, s) => {
        const p = b.data.annotation({ ref: 'p', geometry: 'pairs', target: s, value: { start: [0], b_start: [10] } });
        for (const u of ['u1', 'u2']) b.data.unit({ ref: u, sequence: s });
        const v = b.view({ ref: 'v', label: 'V', copies: 'each', arrangement } as never);
        v.block({ unit: 'u1', range: [0, 8] } as never);
        v.block({ unit: 'u2' });
        v.track({ label: 'P' }).pairs({ data: p });
      });
    const side = await draw(copies('side'), 'v');
    expect(all(side.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([[[9, 19, 1]]]);
    const stacked = await draw(copies('stacked'), 'v');
    expect(all(stacked.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([[], [[1, 11, 1]]]);
    expect(stacked.issues.map((i) => i.message).filter((m) => m.includes('between blocks'))).toEqual([]);
    expect(all(stacked.root, 'nightingale-track')).toEqual([]);
  });

  it('pairs across two blocks of one unit (LA-05) still join: a link side by side, points with the note stacked', async () => {
    const split = (arrangement: 'side' | 'stacked') =>
      small((b, s) => {
        const p = b.data.annotation({ ref: 'p', geometry: 'pairs', target: s, value: { start: [0], b_start: [10] } });
        b.data.unit({ ref: 'u1', sequence: s });
        const v = b.view({ ref: 'v', label: 'V', arrangement } as never);
        v.block({ unit: 'u1', range: [0, 8] } as never);
        v.block({ unit: 'u1', range: [8, SEQ.length] } as never);
        v.track({ label: 'P' }).pairs({ data: p });
      });
    const side = await draw(split('side'), 'v');
    expect(all(side.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([[[1, 11, 1]]]);
    const stacked = await draw(split('stacked'), 'v');
    expect(all(stacked.root, 'nightingale-links').map((x) => x.properties.contacts)).toEqual([[], []]);
    expect(stacked.issues.map((i) => i.message).filter((m) => m.includes('between blocks')).length).toBe(1);
    expect(all(stacked.root, 'nightingale-track').length).toBe(2);
  });

  it('matrix: setHeatmapData with unique state labels and null cells left out; the palette a known loss', async () => {
    const doc = small((b, s) => {
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: s,
        states: ['A', 'B', 'A'],
        value: {
          index: [0, 2],
          value: [
            [0.1, null, 0.3],
            [1, 2, 3],
          ],
        },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'M' }).matrix({ data: m });
    });
    const { root, issues } = await draw(doc, 'v');
    const h = one(root, 'nightingale-sequence-heatmap');
    expect(h.attributes['heatmap-id']).toBe('svs-heatmap-0-0-0');
    const [x, y, cells] = h.calls[0].args as [number[], string[], { xValue: number; yValue: string; score: number }[]];
    expect([h.calls[0].method, x.length, y]).toEqual(['setHeatmapData', 15, ['A', 'B', 'A (2)']]);
    expect(cells.map((c) => [c.xValue, c.yValue])).toEqual([
      [1, 'A'],
      [1, 'A (2)'],
      [3, 'A'],
      [3, 'B'],
      [3, 'A (2)'],
    ]);
    expect(said(issues)).toEqual([]);
  });

  it('matrix (NG-07): a domain clamps the scores and spans the scale (anchors under the first cell); a palette is said', async () => {
    const doc = small((b, s) => {
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: s,
        states: ['A', 'B'],
        value: {
          index: [0, 1],
          value: [
            [1, 2],
            [3, 400],
          ],
        },
      });
      const color = { theme: 'field_palette', params: { palette: { kind: 'continuous', colors: 'Viridis' } } };
      b.view({ ref: 'v', label: 'V' })
        .track({ label: 'M' })
        .matrix({ data: m, domain: [0, 100], color } as never);
    });
    const { root, issues } = await draw(doc, 'v');
    const cells = one(root, 'nightingale-sequence-heatmap').calls[0].args[2] as { xValue: number; score: number }[];
    // 400 clamped to 100; the domain's low end an anchor under the first cell (its high end a cell's).
    expect(cells.map((c) => [c.xValue, c.score])).toEqual([
      [1, 0],
      [1, 1],
      [1, 2],
      [2, 3],
      [2, 100],
    ]);
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      [
        'svs.render.theme',
        'Nightingale: nightingale-sequence-heatmap has its own colors: the color theme is not drawn',
      ],
    ]);
  });

  it('opacity (NG-08): on points and color features; said where Nightingale has none', async () => {
    const doc = small((b, s) => {
      const a = b.data.annotation({ ref: 'a', geometry: 'ranges', target: s, value: { start: [1], stop: [3] } });
      const d = b.data.annotation({
        ref: 'd',
        geometry: 'values',
        target: s,
        value: { start: [1], stop: [3], value: [2] },
      });
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'P' }).points({ data: a, opacity: 0.2 } as never);
      v.track({ label: 'C' }).color({ data: a, opacity: 0.4 } as never);
      v.track({ label: 'V' }).values({ data: d, opacity: 0.5 } as never);
    });
    const { root, issues } = await draw(doc, 'v');
    const [points, colors] = all(root, 'nightingale-track').map((t) => t.properties.data as { opacity?: number }[]);
    expect([points[0].opacity, colors[0].opacity]).toEqual([0.2, 0.4]);
    expect(issues.map((i) => i.message)).toEqual([
      'Nightingale: nightingale-linegraph-track has no opacity: drawn opaque',
    ]);
  });

  it('direction (NG-09) and a ruler of another system (NG-10) are said', async () => {
    const doc = small((b, s) => {
      const a = b.data.annotation({
        ref: 'a',
        geometry: 'ranges',
        target: s,
        value: { start: [1], stop: [3], strand: ['-'] },
      });
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'R' }).ranges({ data: a, direction_field: 'strand' } as never);
      v.track({ label: 'N' }).ruler({ of: s } as never);
    });
    const { issues } = await draw(doc, 'v');
    expect(issues.map((i) => [i.code, i.message])).toEqual([
      ['svs.render.direction', 'Nightingale: nightingale-track draws no direction: every feature reads as forward'],
      [
        'svs.render.fallback',
        'Nightingale: nightingale-navigation: the ruler shows axis positions, not those of the system it is of',
      ],
    ]);
  });

  it('an opacity on letters, a matrix or a ruler is said, as their elements draw none (LA-12)', async () => {
    const doc = small((b, s) => {
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: s,
        states: ['A', 'B'],
        value: {
          index: [0, 1],
          value: [
            [1, 2],
            [3, 4],
          ],
        },
      } as never);
      const v = b.view({ ref: 'v', label: 'V' });
      v.track({ label: 'L' }).letters({ opacity: 0.3 } as never);
      v.track({ label: 'M' }).matrix({ data: m, opacity: 0.3 } as never);
      v.track({ label: 'N' }).ruler({ opacity: 0.3 } as never);
    });
    const { issues } = await draw(doc, 'v');
    expect(issues.map((i) => i.message).filter((m) => m.includes('has no opacity'))).toEqual(
      ['nightingale-sequence', 'nightingale-sequence-heatmap', 'nightingale-navigation'].map(
        (t) => `Nightingale: ${t} has no opacity: drawn opaque`,
      ),
    );
  });

  it("matrix: states in their keyed order and kept by the filter (KS2): the heatmap's rows and cells", async () => {
    const doc = small((b, s) => {
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: s,
        states: ['A', 'B', 'C'],
        value: { index: [0], value: [[1, 2, 3]] },
      });
      b.data.annotation({ ref: 't', geometry: 'rows', value: { state: ['A', 'B', 'C'], k: [2, 1, 3] } } as never);
      b.view({ ref: 'v', label: 'V' })
        .track({ label: 'M' })
        .matrix({ data: m, order: { data: 't', field: 'k' }, filter: [{ data: 't', field: 'k', max: 2 }] } as never);
    });
    const { root, issues } = await draw(doc, 'v');
    const h = one(root, 'nightingale-sequence-heatmap');
    const [, y, cells] = h.calls[0].args as [number[], string[], { yValue: string; score: number }[]];
    expect(y).toEqual(['B', 'A']);
    // The layer's domain is every state's ([1, 3]), as the viewer's scale: an anchor at 3 under the first cell.
    expect(cells.map((c) => [c.yValue, c.score])).toEqual([
      ['A', 3],
      ['A', 1],
      ['B', 2],
    ]);
    expect(said(issues)).toEqual([]);
  });

  it("matrix: cells of one score get an anchor at 0 under the first (the heatmap's color scale throws on one value)", async () => {
    const doc = small((b, s) => {
      const m = b.data.annotation({
        ref: 'm',
        geometry: 'matrix',
        target: s,
        states: ['A', 'B'],
        value: {
          index: [0, 2],
          value: [
            [1, 1],
            [1, null],
          ],
        },
      });
      b.view({ ref: 'v', label: 'V' }).track({ label: 'M' }).matrix({ data: m });
    });
    const h = one((await draw(doc, 'v')).root, 'nightingale-sequence-heatmap');
    const cells = h.calls[0].args[2] as { xValue: number; yValue: string; score: number }[];
    // First the anchor at the first cell's place; the cells after it (the grid keeps the last item at a place).
    expect(cells.map((c) => [c.xValue, c.yValue, c.score])).toEqual([
      [1, 'A', 0],
      [1, 'A', 1],
      [1, 'B', 1],
      [3, 'A', 1],
    ]);
  });

  it("6: the fingerprint's matrix fallback with the interactions extension: a heatmap whose scores span two values", async () => {
    const reg = registries();
    register(reg, [InteractionsRegistry]);
    const doc = CASES.case6();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const root = materialize(resolveView(doc, set, 'poses-view', reg, NightingaleSupported).view, set).root;
    const cells = one(root, 'nightingale-sequence-heatmap').calls[0].args[2] as { score: number }[];
    const scores = cells.map((c) => c.score);
    expect(cells.length).toBeGreaterThan(1);
    expect(Math.max(...scores)).toBeGreaterThan(Math.min(...scores));
  });

  it('rows and logo: nightingale-msa rows in lane order, a scheme by dictionary; a profile as probabilities', async () => {
    const b = createSeqViewSpec({ title: 'aln' });
    const aln = b.data.alignment({
      ref: 'a',
      value: ['MK-LA', 'MKVLA', 'MR-LG'],
      labels: ['x', 'y', 'z'],
      parent_rows: [null, 0, null],
    });
    const prof = b.data.computed({ ref: 'prof', provider: 'profile', alignment: aln });
    const v = b.view({ ref: 'v', label: 'V' }).block({ alignment: aln });
    v.track({ label: 'Logo' }).logo({ data: prof, scale: 'probability' });
    v.track({ label: 'Rows' }).rows({ row_subset: [0, 1], color: { theme: 'letter', params: { colors: 'Taylor' } } });
    const { root, issues } = await draw(b.getDocument(), 'v');
    const msa = one(root, 'nightingale-msa');
    expect(msa.properties.data).toEqual([
      { name: 'x', sequence: 'MK-LA' },
      { name: 'y', sequence: 'MKVLA' },
    ]);
    expect(msa.attributes['color-scheme']).toBe('taylor');
    const logo = one(root, 'nightingale-conservation-track').properties.data as {
      index: number[];
      probabilities: Record<string, number[]>;
    };
    expect(logo.index).toEqual([1, 2, 3, 4, 5]);
    expect(logo.probabilities.M[0]).toBeCloseTo(1);
    expect(logo.probabilities.K[1]).toBeCloseTo(2 / 3);
    expect(said(issues)).toEqual(['svs.render.fallback v.rows-2']); // variant indentation
  });

  it("rows in a keyed set's order and kept mask (KS1): the permutation and the filter reach nightingale-msa", async () => {
    const rows = ['ACDEFGHIKL', 'ACDEFGHIKM', 'ACDEFGHIAM', 'ACDEFGHIKV', 'WYWYWYWYWY', 'WYWYWY----'];
    const msaOf = async (layer: Record<string, unknown>) => {
      const b = createSeqViewSpec({ title: 'aln' });
      const aln = b.data.alignment({ ref: 'a', value: rows, labels: rows.map((_, i) => `r${i}`) });
      b.data.computed({ ref: 'stats', provider: 'row_stats', alignment: aln } as never);
      b.data.computed({ ref: 'nr', provider: 'redundancy', alignment: aln, threshold: 0.8 } as never);
      b.view({ ref: 'v', label: 'V' })
        .block({ alignment: aln })
        .track({ label: 'Rows' })
        .rows(layer as never);
      const { root } = await draw(b.getDocument(), 'v');
      return (one(root, 'nightingale-msa').properties.data as { name: string }[]).map((r) => r.name);
    };
    expect(await msaOf({ order: { data: 'stats', field: 'identity', direction: 'desc' } })).toEqual([
      'r0',
      'r1',
      'r3',
      'r2',
      'r4',
      'r5',
    ]);
    expect(await msaOf({ filter: [{ data: 'nr', field: 'representative', values: [true] }] })).toEqual(['r0', 'r4']);
    expect(
      await msaOf({
        order: { data: 'stats', field: 'identity' },
        filter: [{ data: 'stats', field: 'identity', min: 0.85 }],
      }),
    ).toEqual(['r1', 'r3', 'r0']);
  });

  it("rows: Clustal X known by the table it was expanded from, not by its rules' shape", async () => {
    const scheme = async (color: unknown) => {
      const b = createSeqViewSpec({ title: 'aln' });
      const aln = b.data.alignment({ ref: 'a', value: ['MKVLA', 'MRVLG'], labels: ['x', 'y'] });
      b.view({ ref: 'v', label: 'V' })
        .block({ alignment: aln })
        .track({ label: 'Rows' })
        .rows({ color: color as never });
      const { root, issues } = await draw(b.getDocument(), 'v');
      return [one(root, 'nightingale-msa').attributes['color-scheme'], issues.map((i) => i.message).join('; ')];
    };
    const [clustal, said] = await scheme({ theme: 'clustal_x' });
    expect(clustal).toBe('clustal2');
    expect(said).toContain('Clustal X');
    // A hand-written `rules` theme with as many `share` rules as the table is not the table.
    const share = {
      when: [{ measure: { share: 'letter' }, min: 0.5 }],
      theme: { theme: 'uniform', params: { color: 'red' } },
    };
    const look = Array.from({ length: clustalXRules().length }, () => share);
    expect((await scheme({ theme: 'rules', params: { rules: look } }))[0]).not.toBe('clustal2');
    // A table over another names itself in the issue: `identity_to_row` over Clustal X draws its base, Clustal X.
    const [base, over] = await scheme({ theme: 'identity_to_row', params: { base: { theme: 'clustal_x' } } });
    expect(base).toBe('clustal2');
    expect(over).toContain('the "identity_to_row" theme (its base is drawn)');
    // Color schemes (seqspec:format#color-schemes): nightingale-msa's own where it has one, else the theme the scheme
    // stands for.
    const named = async (params: Record<string, unknown>) => (await scheme({ theme: 'scheme', params }))[0];
    expect(await named({ name: 'clustal_x' })).toBe('clustal2');
    expect(await named({ name: 'taylor' })).toBe('taylor');
    expect(await named({ name: 'hydrophobicity' })).toBe('hydro');
    const [own, lost] = await scheme({ theme: 'scheme', params: { name: 'taylor', overrides: { K: '#000' } } });
    expect(own).toBe('zappo');
    expect(lost).toContain('does not draw these letter colors');
    // A theme nightingale-msa has no scheme for is drawn in a default one, and said (NG-03).
    for (const theme of [
      { theme: 'percent_identity' },
      { theme: 'uniform', params: { color: 'red' } },
      { theme: 'scheme', params: { name: 'percent_identity' } },
    ]) {
      const [drawn, note] = await scheme(theme);
      expect([drawn, note]).toEqual(['zappo', expect.stringMatching(/has no "[a-z_]+" scheme: drawn in zappo/)]);
    }
  });

  it('letters (NG-06): a letter the dictionary leaves uncolored takes the ground, not its first color', async () => {
    const b = createSeqViewSpec({ title: 'x' });
    const s = b.data.sequence({ ref: 's', value: 'MKXLAZ*B' });
    b.view({ ref: 'v', label: 'V' })
      .block({ sequence: s })
      .track({ label: 'L' })
      .letters({ color: { theme: 'letter', params: { colors: 'Zappo' } } } as never);
    const c = one((await draw(b.getDocument(), 'v')).root, 'nightingale-colored-sequence');
    const range = c.attributes['color-range'].split(',');
    const at = Object.fromEntries(c.attributes.scale.split(',').map((x) => [x.slice(0, x.lastIndexOf(':')), x]));
    const ground = range.findIndex((x) => x.startsWith('#ffffff:'));
    expect(ground).toBeGreaterThan(0);
    for (const k of ['X', 'Z', '*', 'B']) expect(at[k]).toBe(`${k}:${ground}`);
    expect(at.M).not.toBe(`M:${ground}`);
  });
});

describe('view features', () => {
  const twoBlocks = (params: Record<string, unknown>) =>
    small((b, s) => {
      b.data.sequence({ ref: 't', value: 'ACDEFG' });
      const r = b.data.annotation({ ref: 'r', geometry: 'ranges', target: 't', value: { start: [1], stop: [3] } });
      const v = b
        .view({ ref: 'v', label: 'V', ...params })
        .block({ sequence: s })
        .block({ sequence: 't' });
      v.focus({ target: 't', start: 2, stop: 4 });
      v.track({ label: 'Res' }).letters();
      v.group({ label: 'G' }).track({ label: 'R' }).ranges({ data: r });
      v.track({ label: 'Hidden', hidden: true }).letters();
    });

  it('side by side: one manager over the concatenated blocks, a block label track, focus as the display range', async () => {
    const { root, issues, view, set } = await draw(twoBlocks({}), 'v');
    const managers = all(root, 'nightingale-manager');
    expect(managers.length).toBe(1);
    // The focus opens every element, never the manager, which would re-apply it after each zoom (NG-01).
    expect(managers[0].attributes['display-start']).toBeUndefined();
    const opened = all(root, 'nightingale-track').map((t) => [
      t.attributes['display-start'],
      t.attributes['display-end'],
    ]);
    expect(new Set(opened.map((x) => x.join(':')))).toEqual(new Set(['18:19']));
    // The block labels, group G's merged row, its track R.
    expect(all(root, 'nightingale-track').map((t) => t.attributes.length)).toEqual(['21', '21', '21']);
    expect(one(root, 'nightingale-sequence').properties.data).toBe(`${SEQ}ACDEFG`);
    // Block names are laid over the bars (mount's label layer).
    expect(data(root, 'nightingale-track').map((f) => [f.start, f.end, f.tooltipContent, f.svs])).toMatchObject([
      [1, 15, 'S', { layer: -1, label: 'S', inside: true }],
      [16, 21, 't', { layer: -1, label: 't', inside: true }],
    ]);
    expect(LABELED in one(root, 'nightingale-track').attributes).toBe(true);
    expect(all(root, 'nightingale-track')[1].properties.data).toMatchObject([{ start: 17, end: 18 }]);
    expect(all(root, 'span').map((d) => d.properties.textContent)).toEqual(['', 'Res', '▾', 'G', 'R']);
    expect(said(issues)).toEqual(['svs.render.fallback v.main']);
    // The axis back to blocks and indices; a host's loci to `highlight`.
    const [axis] = axesOf(view, set, 0);
    expect([locate(axis, 0), locate(axis, 16)]).toEqual([
      { block: 0, index: 0 },
      { block: 1, index: 1 },
    ]);
    const unit = view.blocks.unit[1];
    const loci = {
      kind: 'seq-loci' as const,
      set,
      elements: [{ unit: unit as never, indices: OrderedSet.ofRange(1, 4) as never }],
    };
    expect(highlightOf(view, set, axis, loci)).toBe('17:19');
    // A hovered item's extent (its fragments, across the block edge) as block runs, and back as one highlight.
    const extent = extentOf(axis, '14:17,20:20');
    expect(extent).toEqual([
      { block: 0, start: 13, stop: 15 },
      { block: 1, start: 0, stop: 2 },
      { block: 1, start: 4, stop: 5 },
    ]);
    expect(highlightOf(view, set, axis, extentLoci(view, set, extent))).toBe('14:17,20:20');
    expect([extentOf(axis, ''), extentLoci(view, set, [])]).toEqual([[], null]);
  });

  it("stacked: one manager per block, the longest's length; positions cut the axis", async () => {
    const stacked = await draw(twoBlocks({ arrangement: 'stacked' }), 'v');
    expect(all(stacked.root, 'nightingale-manager').length).toBe(2);
    expect(all(stacked.root, 'nightingale-sequence').map((m) => m.attributes.length)).toEqual(['15', '15']);
    expect(said(stacked.issues)).toEqual(['svs.render.fallback v.main']);
    const cut = await draw(twoBlocks({ positions: 'r' }), 'v');
    expect(one(cut.root, 'nightingale-sequence').attributes.length).toBe('2');
    expect(one(cut.root, 'nightingale-sequence').properties.data).toBe('CD');
    // A cut axis back to block indices: the shown subset's, split where they jump.
    expect(extentOf(axesOf(cut.view, cut.set, 0)[0], '1:2')).toEqual([{ block: 0, start: 1, stop: 3 }]);
    const shown = OrderedSet.ofSortedArray(Int32Array.of(2, 3, 7));
    const b = {
      block: 0,
      offset: 0,
      count: 3,
      first: 0,
      ligand: false,
      shown,
      reversed: false,
      layout: null,
      unit: 0,
      members: null,
      alignment: -1,
    };
    expect(extentOf({ section: 0, scrollGroup: null, length: 3, blocks: [b] }, '1:3')).toEqual([
      { block: 0, start: 2, stop: 4 },
      { block: 0, start: 7, stop: 8 },
    ]);
  });

  it("pieces (AX1): a range and a reversed piece, filtered and not, drawn in the piece's order", async () => {
    const doc = (o: { direction?: 1 | -1; positions?: string }) =>
      small((b, s) => {
        b.data.annotation({ ref: 'f', geometry: 'ranges', target: s, value: { start: [4, 9], stop: [6, 10] } });
        const v = b.view({ ref: 'v', label: 'V', ...(o.positions ? { positions: o.positions } : {}) });
        v.block({ sequence: s, range: [3, 11], direction: o.direction ?? 1 });
        v.track({ label: 'L' }).letters();
        v.track({ label: 'F' }).ranges({ data: 'f' });
      });
    const fwd = await draw(doc({}), 'v');
    expect(one(fwd.root, 'nightingale-sequence').properties.data).toBe(SEQ.slice(3, 11));
    const rev = await draw(doc({ direction: -1 }), 'v');
    expect(one(rev.root, 'nightingale-sequence').properties.data).toBe([...SEQ.slice(3, 11)].reverse().join(''));
    // [4, 6) and [9, 10) of [3, 11) reversed: axis positions 6..7 and 2 (1-based).
    const frags = (r: typeof rev) => data(r.root, 'nightingale-track').map((f) => f.locations);
    expect(frags(rev)).toEqual([[{ fragments: [{ start: 6, end: 7 }] }], [{ fragments: [{ start: 2, end: 2 }] }]]);
    const axis = axesOf(rev.view, rev.set, 0)[0];
    expect([axis.length, locate(axis, 0)]).toEqual([8, { block: 0, index: 10 }]);
    expect(extentOf(axis, '1:3')).toEqual([{ block: 0, start: 8, stop: 11 }]);
    // Filtered: the filter's positions in the piece, from the last.
    const cut = await draw(doc({ direction: -1, positions: 'f' }), 'v');
    expect(one(cut.root, 'nightingale-sequence').properties.data).toBe(SEQ[9] + SEQ[5] + SEQ[4]);
    const cutAxis = axesOf(cut.view, cut.set, 0)[0];
    expect(extentOf(cutAxis, '1:3')).toEqual([
      { block: 0, start: 9, stop: 10 },
      { block: 0, start: 4, stop: 6 },
    ]);
    expect(frags(cut)).toEqual([[{ fragments: [{ start: 2, end: 3 }] }], [{ fragments: [{ start: 1, end: 1 }] }]]);
  });

  it("letter maps (LM1): letters, an alignment's rows and a polymer's letters `of` are the letters shown, reversed too", async () => {
    const doc = (map: unknown, direction: 1 | -1 = 1) =>
      small((b, s) => {
        const aln = b.data.alignment({ ref: 'a', value: ['MKVLA', 'MRVLG'], labels: ['x', 'y'] });
        const v = b.view({ ref: 'v', label: 'V' });
        v.section({ ref: 'v.s', direction })
          .block({ sequence: s, range: [3, 11] })
          .track({ label: 'L' })
          .letters({
            letter_map: map,
          } as never);
        v.section({ ref: 'v.r' })
          .block({ alignment: aln })
          .track({ label: 'R' })
          .rows({ letter_map: map } as never);
      });
    const lower = SEQ.slice(3, 11).toLowerCase();
    const fwd = await draw(doc('lower'), 'v');
    expect(all(fwd.root, 'nightingale-sequence').map((e) => e.properties.data)).toContain(lower);
    expect(data(fwd.root, 'nightingale-msa' as never).map((r) => r.sequence)).toEqual(['mkvla', 'mrvlg']);
    const rev = await draw(doc({ L: 'l', V: 'v', A: 'a' }, -1), 'v');
    expect(one(rev.root, 'nightingale-sequence').properties.data).toBe(
      [...SEQ.slice(3, 11)].reverse().join('').replace(/L/g, 'l').replace(/V/g, 'v').replace(/A/g, 'a'),
    );
    const of = small((b, s) => {
      b.data.sequence({ ref: 'other', value: 'GGGKKK' });
      b.data.mapping({ ref: 'm', from: 'other', to: s, method: 'offset', offset: 2 } as never);
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: s })
        .track({ label: 'L' })
        .letters({
          of: 'other',
          letter_map: { K: 'k' },
        } as never);
    });
    expect(one((await draw(of, 'v')).root, 'nightingale-sequence').properties.data).toBe('--GGGkkk-------');
    // Through a codon step, each residue a feature over its codon: its color, tooltip and letter by the letter shown.
    const codons = (map: unknown, dna = 'ATGAAGGAC') => {
      const b = createSeqViewSpec();
      const g = b.data.sequence({ ref: 'g', alphabet: 'dna', value: dna });
      const p = translate(b.data, { ref: 'p', from: g, cds: [[0, 9]] });
      b.view({ ref: 'v', label: 'V' })
        .block({ sequence: g })
        .track({ label: 'P' })
        .letters({ of: p, letter_map: map } as never);
      return b.getDocument();
    };
    const residues = async (d: SeqViewSpecDoc) =>
      (data((await draw(d, 'v')).root, 'nightingale-track') as { tooltipContent: string; fill: string }[]).map((f) => [
        f.tooltipContent,
        f.fill,
      ]);
    const mapped = await residues(codons({ K: 'D' }));
    expect(mapped.map((f) => f[0])).toEqual(['Met 1', 'D 2', 'Asp 3']);
    expect(mapped[1][1]).toBe(mapped[2][1]);
    expect(mapped[1][1]).not.toBe((await residues(codons(null)))[1][1]);
    // Each residue's one-letter code laid inside its codon by mount (a labeled track), the hint naming the residue.
    const letters = async (map: unknown) => {
      const track = one((await draw(codons(map), 'v')).root, 'nightingale-track');
      const svs = (track.properties.data as { svs: { label?: string; inside?: boolean; hint?: string } }[]).map(
        (f) => f.svs,
      );
      expect(track.attributes[LABELED]).toBe('');
      expect(svs.every((x) => x.inside)).toBe(true);
      return svs.map((x) => [x.label, x.hint]);
    };
    expect(await letters(null)).toEqual([
      ['M', 'Met 1'],
      ['K', 'Lys 2'],
      ['D', 'Asp 3'],
    ]);
    expect(await letters({ K: 'D' })).toEqual([
      ['M', 'Met 1'],
      ['D', 'D 2'],
      ['D', 'Asp 3'],
    ]);
    // No map: as stored.
    expect(one((await draw(doc(null), 'v')).root, 'nightingale-sequence').properties.data).toBe(SEQ.slice(3, 11));
  });

  it('text: one unwrapped axis, reported once for the view', async () => {
    const doc = small((b) => b.view({ ref: 'v', label: 'V', layout: 'text' }).track({ label: 'S' }).letters());
    expect(said((await draw(doc, 'v')).issues)).toEqual(['svs.render.fallback v']);
  });
});

/**
 * Every case's views as `code ref` per issue materialize reports: what Nightingale cannot draw of what the document
 * chose. Losses every layer of a type has (seqspec:libraries#known-losses) are not reported, and the
 * resolution's own issues (types and glyphs drawn as their fallbacks) are resolveView's and not listed here.
 */
const GOLDEN: Record<string, Record<string, string[]>> = {
  case1: {
    text: ['fallback text', 'fallback text.main'],
    tracks: [],
  },
  case2: { fv: ['fallback fv.main', 'fallback fv.ruler-2'] },
  case3: {
    family: ['theme family.rows', 'fallback family.main'],
  },
  case4: { rna: ['fallback rna.pairs', 'fallback rna.pairs-2'] },
  case5: {
    map: ['fallback map.ranges-2'],
    seq: ['fallback seq'],
  },
  case6: {
    chains: ['fallback chains.main', 'fallback chains.ruler-2', 'fallback chains.pairs', 'direction chains.pairs'],
    'poses-view': ['fallback poses-view.ruler-2'], // the pocket's numbers are not its axis positions
  },
  case7: Object.fromEntries(
    ['complex-one', 'complex-each'].map((v) => [v, [`fallback ${v}.main`, `fallback ${v}.ruler-2`]]),
  ),
  // The pLDDT line in its bands: a linegraph has one color.
  case8: { chains: ['fallback chains.main', 'fallback chains.ruler-2', 'theme chains.values'] },
  // Conservation: streamed, fetched by mount; the rows track's two layers (rows, the bound rows' helices) apart, said once for the section.
  case9: {
    msa: ['svs.stream.fallback msa.values', 'fallback msa.rows-2', 'fallback msa.main', 'fallback msa.chains'],
  },
  case10: {
    chains: ['fallback chains.main', 'fallback chains.pairs', 'direction chains.pairs'],
    alignment: ['theme alignment.rows', 'fallback alignment.main'],
    text: ['fallback text', 'fallback text.main'],
    strip: [],
    both: ['fallback both.rows-2'],
  },
  case12: { 'design-view': ['fallback design-view.main'] },
  case13: { features: [] },
  case14: { 'gene-view': ['fallback gene-view.letters-2'] }, // no three-letter residues: residues over their codons are features
  // Indels: each design's row is its layers apart (no compositing), the insertions' and junctions' boundaries bridges
  // over k..k+1, the edits' variants without their deletions and insertions or the split by design, a design's `wt`
  // numbering ruler as axis positions, and `ruler` with `of` (the alignment's WT numbering) as column positions (NG-10).
  case15: {
    wt_frame: [
      'fallback wt_frame.ranges-7',
      'fallback wt_frame.ranges-9',
      'fallback wt_frame.variants',
      'fallback wt_frame.edits-as-variants',
      'fallback wt_frame.main',
    ],
    aligned: ['fallback aligned.ruler-3'],
    own_axes: [
      'fallback own_axes.wt',
      'fallback own_axes.ruler-5',
      'fallback own_axes.d1',
      'fallback own_axes.ruler-8',
      'fallback own_axes.ranges-6',
      'fallback own_axes.d2',
      'fallback own_axes.ruler-11',
      'fallback own_axes.ranges-8',
      'fallback own_axes.d3',
      'fallback own_axes.ruler-14',
      'fallback own_axes.ranges-10',
      'fallback own_axes.d4',
    ],
  },
  // A plasmid map: the ring drawn as a line (`circular` is a style, a hint), the cut sites as bridges over k..k+1 and
  // in the component's colors, the repeats' chords as links without their category, the features apart from the sites.
  case16: {
    pmap: [
      'direction pmap.ranges',
      'theme pmap.ranges-2',
      'fallback pmap.ranges-2',
      'fallback pmap.pairs-2',
      'fallback pmap.ring',
      'direction pmap.ranges-3',
      'fallback pmap.ranges-4',
      'fallback pmap.pairs-4',
    ],
  },
  // A genotype matrix: a heatmap, its rows the states in their keyed order and filter; `density` (a hint) not used; its
  // palette not drawn (the heatmap's own colors, NG-07).
  case17: { 'genotypes-view': ['theme genotypes.matrix'] },
  // Finding things in p53: drawn whole; its location sets light as its elements' `highlight` (the Search part is a host's).
  case18: { tracks: [] },
  // A chromosome streamed by region: 83 Mb, past Nightingale's limit (GM-05): its section left out, said once.
  case19: { chr17: ['omitted chr17.main'] },
  // P1 on its genome: one view, its chr2 section a whole streamed chromosome (242 Mb): left out (GM-05); the protein's
  // and the transcript's render.
  case20: {
    coding: [
      'svs.stream.fallback coding.ranges-2',
      'svs.stream.fallback coding.variants-2',
      'svs.stream.fallback coding.values',
      'omitted coding.genome',
    ],
  },
  // A gene on its reverse strand: rulers as axis positions (a descending piece), the protein's codons as features, the
  // junctions' arcs as links.
  case21: {
    gene2: [
      'fallback gene2.ruler-2',
      'direction gene2.ranges',
      'fallback gene2.letters-2',
      'fallback gene2.pairs-2',
      'fallback gene2.ruler-4',
      'fallback gene2.letters-4',
      'direction gene2.ranges-3',
      'fallback gene2.pairs-4',
    ],
  },
  // A family triaged: the rows with their data columns and tree as Nightingale's MSA (no columns, no tree).
  case22: { family: ['fallback family.members'] },
};

describe('the cases, drawn', () => {
  const reg = registries();
  const drawn = async (name: string, view: string, edit?: (doc: SeqViewSpecDoc) => SeqViewSpecDoc) => {
    const doc = (edit ?? ((d) => d))(CASES[name]());
    const { set } = await compile(doc, reg, fixtureLoaders());
    return materialize(resolveView(doc, set, view, reg, NightingaleSupported).view, set).root;
  };

  it("7: merged copies' features once each, helices as outlined coils, one axis under the ruler", async () => {
    const root = await drawn('case7', 'complex-one');
    const ss = all(root, 'nightingale-track').find((t) => t.attributes['data-svs-layer'] === '2');
    const features = (ss?.properties.data ?? []) as { start: number; end: number; shape: string; fill?: string }[];
    const places = features.map((f) => `${f.start}-${f.end}`);
    expect(new Set(places).size).toBe(places.length); // chains A, B, C: one feature per place, not three lanes
    expect(features.filter((f) => f.shape === 'helix').every((f) => f.fill === 'transparent')).toBe(true);
    expect(features.some((f) => f.shape === 'strand')).toBe(true);
    // The interface in the document's one color (no categories to tell its items apart by).
    const iface = all(root, 'nightingale-track').find((t) => t.attributes['data-svs-layer'] === '3');
    const colors = ((iface?.properties.data ?? []) as { fill: string }[]).map((f) => f.fill);
    expect(new Set(colors)).toEqual(new Set(['#d1495b']));
    // The contacts' one color, for mount's style over nightingale-links' own ramp.
    expect(one(root, 'nightingale-links').attributes.style).toBe('--svs-pair-color: #4a6fa5');
    expect(one(root, 'nightingale-sequence').attributes.class).toBe(AXISLESS);
    // The heatmap on the shown chain (the effects fixture), within a bounded height.
    const h = one(root, 'nightingale-sequence-heatmap');
    const cells = h.calls[0].args[2] as { xValue: number }[];
    expect(cells.length).toBeGreaterThan(0);
    expect(Number(h.attributes.height)).toBeLessThanOrEqual(160);
    // Hydropathy's provider is not registered: the layer is placed nowhere and draws nothing.
    expect(all(root, 'nightingale-linegraph-track')).toEqual([]);
  });

  it("7: the UniProt track split by feature type: a group's rows, collapsed at open, the layer whole beside its header", async () => {
    const root = await drawn('case7', 'complex-one');
    const g = all(root, 'div').find((d) => d.attributes['data-svs-group'] === 'complex-one.uniprot') as ElementDesc;
    const [header, body] = g.children;
    expect('data-collapsed' in g.attributes).toBe(true);
    const text = (d: ElementDesc) => all(d, 'span').map((x) => String(x.properties.textContent));
    expect(text(header).slice(0, 2)).toEqual(['▸', 'UniProt']);
    type F = { svs: { item: number }; tooltipContent: string };
    const whole = all(header, 'nightingale-track')[0].properties.data as F[];
    // A row per type, labeled by it, indented; each its type's features alone, together the layer's.
    const rows = body.children.map((r) => [text(r)[0], all(r, 'nightingale-track')[0].properties.data as F[]] as const);
    expect(rows.map(([label]) => label)).toContain('Chain');
    expect(body.children.every((r) => r.attributes.style === '--svs-depth: 1')).toBe(true);
    for (const [label, fs] of rows) expect(fs.every((f) => f.tooltipContent.includes(label))).toBe(true);
    expect(rows.reduce((n, [, fs]) => n + fs.length, 0)).toBe(whole.length);
  });

  it('split `empty: "hide"`: a split track per value present, the variants grid condensed', async () => {
    const hide = (order: string[] | null) => (doc: SeqViewSpecDoc) => {
      const walk = (n: DocNode): DocNode => {
        const sp = n.params?.split as { by: string } | null | undefined;
        const split = sp ? { ...sp, order: order ?? (sp as { order?: string[] }).order, empty: 'hide' } : sp;
        return { ...n, params: sp ? { ...n.params, split } : n.params, children: n.children?.map(walk) };
      };
      return { ...doc, root: walk(doc.root) };
    };
    const g = (root: ElementDesc) =>
      all(root, 'div').find((d) => d.attributes['data-svs-group'] === 'complex-one.uniprot') as ElementDesc;
    const labels = (root: ElementDesc) =>
      g(root).children[1].children.map((r) => all(r, 'span')[0].properties.textContent);
    const shown = labels(await drawn('case7', 'complex-one', hide(['No such type', 'Chain'])));
    expect(shown[0]).toBe('Chain');
    expect(shown).not.toContain('No such type');
    const grid = all(await drawn('case13', 'features', hide(null)), 'nightingale-variation-canvas')[0];
    expect(grid.attributes['condensed-view']).toBe('');
    expect(all(await drawn('case13', 'features'), 'nightingale-variation-canvas')[0].attributes['condensed-view']).toBe(
      undefined,
    );
  });

  it('5: each kind of feature its own color by default; a table without categories colored per item', async () => {
    const root = await drawn('case5', 'map');
    const [features, sites] = all(root, 'nightingale-track').map((t) => t.properties.data as Record<string, unknown>[]);
    const color = Object.fromEntries(features.map((f) => [f.tooltipContent, f.color]));
    expect(color['start codon · start_codon']).toBeDefined();
    expect(color['EcoRI · restriction_site']).toBe(color['XhoI · restriction_site']);
    expect(new Set(Object.values(color)).size).toBe(3);
    expect(new Set(sites.map((f) => f.color)).size).toBe(sites.length);
    expect(sites.every((f) => typeof f.color === 'string')).toBe(true);
    // `underline` is the text layout's: on the plain axis the features are rectangles, no loss of their own.
    const text = all(await drawn('case5', 'seq'), 'nightingale-track')[0].properties.data as { shape: string }[];
    expect(new Set(text.map((f) => f.shape))).toEqual(new Set(['rectangle']));
  });

  it("14: residues as features over their codons, protein data as whole codons on the gene's axis", async () => {
    const root = await drawn('case14', 'gene-view');
    const tracks = all(root, 'nightingale-track');
    const layer = (k: number) => tracks.find((t) => t.attributes['data-svs-layer'] === String(k))?.properties.data;
    type F = { start: number; end: number; locations: { fragments: { start: number; end: number }[] }[] };
    const residues = layer(2) as (F & { tooltipContent: string; fill: string })[];
    expect(residues.length).toBe(20);
    expect([residues[0].start, residues[0].end, residues[0].tooltipContent]).toEqual([13, 15, 'Met 1']);
    expect([residues[11].start, residues[11].end, residues[11].locations[0].fragments]).toEqual([
      46,
      78,
      [
        { start: 46, end: 46 },
        { start: 77, end: 78 },
      ],
    ]);
    expect(residues.every((f) => f.fill?.startsWith('#'))).toBe(true);
    // Each residue's letter, laid inside its codon by mount.
    const svs = (residues as unknown as { svs: { label: string; inside: boolean } }[]).map((f) => f.svs);
    expect([svs[0].label, svs.every((x) => x.inside && x.label.length === 1)]).toEqual(['M', true]);
    const [domain] = layer(4) as F[];
    expect([domain.start, domain.end, domain.locations[0].fragments]).toEqual([
      25,
      84,
      [
        { start: 25, end: 46 },
        { start: 77, end: 84 },
      ],
    ]);
    const [site] = layer(5) as F[];
    expect([site.start, site.end]).toEqual([91, 93]);
    expect(all(root, 'nightingale-sequence').length).toBe(1); // the gene's own letters; none over codons
  });

  it("13: groups as header rows, collapsed ones with their tracks' features merged beside them, tracks indented", async () => {
    const root = await drawn('case13', 'features');
    const groups = all(root, 'div').filter((d) => d.attributes.class === 'svs-group');
    const text = (d: ElementDesc) => all(d, 'span').map((x) => x.properties.textContent);
    expect(groups.map((g) => [text(g.children[0]).slice(0, 2), 'data-collapsed' in g.attributes])).toEqual([
      [['▾', 'Molecule processing'], false],
      [['▸', 'PTM'], true],
      [['▸', 'Sequence information'], true],
      [['▾', 'Domains & sites'], false],
      [['▾', 'Variants'], false],
    ]);
    // The variants grid (split by alt): nightingale-variation's own, a row per residue, not a track per value.
    expect(all(groups[4], 'nightingale-variation-canvas').length).toBe(1);
    expect(all(groups[4], 'div').filter((d) => d.attributes['data-svs-group']).length).toBe(1);
    const [header, body] = groups[1].children;
    expect(header.children[0].children[0].attributes).toMatchObject({ class: 'svs-toggle', 'aria-expanded': 'false' });
    expect(body.children.map((r) => [r.attributes.style, text(r)[0]])).toEqual([
      ['--svs-depth: 1', 'Modified residue'],
      ['--svs-depth: 1', 'Cross-link'],
    ]);
    // The merged row: every member's features (copies), each still naming its layer and item.
    const members = all(body, 'nightingale-track').flatMap(
      (x) => x.properties.data as { start: number; svs: object }[],
    );
    const [merged] = header.children[1].children;
    expect(merged.attributes).toMatchObject({ length: '393', layout: 'non-overlapping' });
    expect(merged.attributes['data-svs-layer']).toBeUndefined();
    const data = merged.properties.data as { start: number; svs: object }[];
    expect(data).toEqual(members);
    expect(data[0]).not.toBe(members[0]);
    // Variants merged as points; a track's hint is its data's provenance.
    const variants = groups[4].children[0].children[1].children[0].properties.data as {
      start: number;
      shape: string;
    }[];
    expect(variants.map((f) => [f.start, f.shape]).slice(0, 2)).toEqual([
      [47, 'circle'],
      [72, 'circle'],
    ]);
    const hints = all(root, 'span').filter((x) => x.attributes.class === 'svs-hint');
    expect(hints.at(-1)?.attributes['data-svs-hint']).toBe(
      'Proteins API · imported · polymorphisms and somatic hotspots',
    );
    // Domains, regions and motifs ask for their labels inside (`labels`); the other tracks keep them for the hint.
    const rowOf = (label: string) =>
      all(root, 'div').find((d) => d.attributes.class === 'svs-track' && text(d)[0] === label);
    const labeled = (label: string) => all(rowOf(label) as ElementDesc, 'nightingale-track')[0].attributes;
    expect(['Domain', 'Region', 'Motif', 'Chain'].map((l) => LABELED in labeled(l))).toEqual([true, true, true, false]);
    // A variation track: one wider left gutter on every element, for its residue axis, so all keep one axis.
    const margins = new Set(
      ['nightingale-navigation', 'nightingale-sequence', 'nightingale-track', 'nightingale-variation-canvas']
        .flatMap((t) => all(root, t))
        .map((x) => `${x.attributes['margin-left']}/${x.attributes['margin-right']}`),
    );
    expect([...margins]).toEqual(['30/14']);
  });

  it("9: a lazy alignment's rows drawn as their fallback, the reference row's letters (fetched whole)", async () => {
    const root = await drawn('case9', 'msa');
    const rows = all(root, 'nightingale-sequence').find((e) => e.attributes['data-svs-layer'] === '2');
    expect(rows?.properties.data).toBe(proceduralRow(0, 320));
    // Labeled by what it draws, the details in its hint; no axis of its own (the view has no ruler).
    const row = all(root, 'div').find((d) => d.attributes.class === 'svs-track' && all(d, 'nightingale-sequence')[0]);
    const spans = all(row as ElementDesc, 'span');
    expect(spans.map((x) => x.properties.textContent)).toEqual([
      'SEQ0/1-265 (reference row; 17,999 rows not drawn)',
      'i',
    ]);
    expect(spans[1].attributes['data-svs-hint']).toMatch(/^Track "Rows": the alignment's 18,000 rows load as/);
    expect(rows?.attributes.class).toBe(AXISLESS);
  });

  it("2: hovering a CDR lights all of it: its fragments as the extent, the extent as every block's highlight", async () => {
    const doc = CASES.case2();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const view = resolveView(doc, set, 'fv', reg, NightingaleSupported).view;
    const root = materialize(view, set).root;
    const cdrs = all(root, 'nightingale-track').find((t) => t.attributes['data-svs-layer'] === '2');
    const features = (cdrs?.properties.data ?? []) as { start: number; end: number; tooltipContent: string }[];
    const h3 = features.find((f) => f.tooltipContent.startsWith('CDR-H3'));
    const axis = one(root, 'nightingale-manager').properties.svs as Parameters<typeof extentOf>[0];
    const extent = extentOf(axis, `${h3?.start}:${h3?.end}`);
    expect(extent.length).toBe(1);
    expect(extent[0].stop - extent[0].start).toBe((h3?.end ?? 0) - (h3?.start ?? 0) + 1);
    expect(highlightOf(view, set, axis, extentLoci(view, set, extent))).toBe(`${h3?.start}:${h3?.end}`);
    // Every CDR, of either block (VL's past VH's on the one axis), and a residue of each: what the story lights.
    const blocks = new Set<number>();
    for (const f of features) {
      const runs = extentOf(axis, `${f.start}:${f.end}`);
      for (const r of runs) blocks.add(r.block);
      expect(highlightOf(view, set, axis, extentLoci(view, set, runs))).toBe(`${f.start}:${f.end}`);
      expect(highlightOf(view, set, axis, extentLoci(view, set, extentOf(axis, `${f.end}:${f.end}`)))).toBe(
        `${f.end}:${f.end}`,
      );
    }
    expect(blocks.size).toBe(2);
  });

  it('9: conservation, streamed, described for mount; its bins fetched at the finest level a pixel holds', async () => {
    const root = await drawn('case9', 'msa');
    const line = one(root, 'nightingale-linegraph-track');
    const p = line.properties.svs as StreamedValues;
    expect([
      p.ref,
      p.domain,
      p.parts.map((x) => [x.stream.source.data_kind, x.stream.source.declared.measure]),
    ]).toEqual(['msa.values', [0, 1], [['values', 'conservation']]]);
    // Nothing until mount fetches: every value missing.
    const [series] = line.properties.data as { values: { value: number | null }[] }[];
    expect(series.values.every((v) => v.value === null)).toBe(true);
    // The finest level whose bins are a pixel wide, else the coarsest.
    const widths = (n: number) => Array.from({ length: n }, (_, k) => 2 ** k);
    expect([
      streamedLevel(widths(2), 320, 1000),
      streamedLevel(widths(2), 320, 100),
      streamedLevel(widths(8), 18000, 1000),
    ]).toEqual([0, 1, 5]);
    const fine = await fetchStreamed(p, 320, 1000);
    expect([fine.level, fine.issues, fine.values.slice(0, 12)]).toEqual([
      0,
      [],
      [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0, 0.1].map((x) => expect.closeTo(x, 5)),
    ]);
    const coarse = await fetchStreamed(p, 320, 100); // a bin of two columns: the middle of its [min, max]
    expect([coarse.level, coarse.values[0], coarse.values[1], coarse.values[2]]).toEqual([
      1,
      expect.closeTo(0.05, 5),
      expect.closeTo(0.05, 5),
      expect.closeTo(0.25, 5),
    ]);
    // A region that fails is an issue; the others are drawn.
    const stream = p.parts[0].stream;
    const broken = {
      ...stream,
      provider: {
        fetch: (r: StreamRequest, signal: AbortSignal) =>
          r.start === 0 ? Promise.reject(new Error('503')) : stream.provider.fetch(r, signal),
      },
    };
    const half = await fetchStreamed({ ...p, parts: [{ ...p.parts[0], stream: broken }] }, 320, 1000);
    expect(half.issues.map((i) => [i.code, i.path, i.message])).toEqual([
      ['svs.stream.load', 'msa.values', `Nightingale: 1 tile(s) failed: ${stream.source.uri} pf00042 [0, 256): 503`],
    ]);
    expect([half.values[0], half.values[256]]).toEqual([null, expect.closeTo(0.6, 5)]);
    // A reversed piece (NG-05): its first index is its highest; the same tiles, drawn right to left.
    const b0 = p.parts[0].block;
    const reversed = { ...b0, shown: OrderedSet.ofRange(b0.first, b0.first + b0.count), reversed: true };
    const back = await fetchStreamed({ ...p, parts: [{ ...p.parts[0], block: reversed }] }, 320, 1000);
    expect([back.issues, back.values.slice(0, 3)]).toEqual([[], [319, 318, 317].map((k) => fine.values[k])]);
  });
});

describe('the cases (golden issues)', () => {
  const reg = registries();
  it.each(Object.entries(CASES))('%s', async (name, make) => {
    const doc = make();
    const { set } = await compile(doc, reg, fixtureLoaders());
    const views = (doc.root.children ?? []).filter((c) => c.kind === 'view').map((c) => c.ref as string);
    const got = Object.fromEntries(views.map((v) => [v, issuesOf(doc, set, v)]));
    expect(got).toEqual(GOLDEN[name]);
  });
  const issuesOf = (doc: SeqViewSpecDoc, set: SeqSet, v: string) => {
    const view = resolveView(doc, set, v, reg, NightingaleSupported).view;
    const r = materialize(view, set);
    expect(r.issues.every((i) => i.severity === 'warning' && i.message.startsWith('Nightingale: '))).toBe(true);
    return r.issues.map((i) => `${i.code.replace('svs.render.', '')} ${i.path}`);
  };
});

describe("a section's highlights (LS1)", () => {
  it("are its elements' `highlight` (not the manager's), every entry whatever its style; a section without them has none", async () => {
    const doc = small((b, s) => {
      b.data.annotation({
        ref: 'sites',
        geometry: 'ranges',
        target: s,
        value: { start: [2, 9], stop: [4, 12], id: ['a', 'b'] },
      });
      const v = b.view({ ref: 'v', label: 'V' });
      const highlights = [
        { data: 'sites', items: ['b'] },
        { data: 'sites', items: ['a'], style: 'accent' },
      ];
      v.section({ ref: 'one', highlights } as never)
        .block({ sequence: s })
        .track({ label: 'R' })
        .letters();
      v.section({ ref: 'two' }).block({ sequence: s }).track({ label: 'R' }).letters();
    });
    const { root } = await draw(doc, 'v');
    const [one, two] = all(root, 'nightingale-manager');
    const lit = (m: ElementDesc) => [
      m.attributes.highlight,
      ...m.children.flatMap((c) => all(c, 'nightingale-sequence')).map((x) => x.attributes.highlight),
    ];
    expect([lit(one), lit(two)]).toEqual([
      [undefined, '3:4,10:12'],
      [undefined, undefined],
    ]);
    const kept = (one.properties.svs as AxisMeta).highlights as AxisRuns;
    expect(AxisRuns.text(kept)).toBe('3:4,10:12');
    expect((two.properties.svs as AxisMeta).highlights).toBeUndefined();
  });
});
