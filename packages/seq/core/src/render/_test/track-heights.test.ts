/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The defaults decided from the audit (_spec/README.md "Render", Track heights and Labels): a track fits its content in
// steps of a lane, tracks of one kind share the tallest, the section's `heights` names the sizing and the host's
// `trackSizing` port has the last word; labels inside their mark when they fit whole, else in a label lane above the
// track, else dropped; the label column measured; sections 20 px apart; a colourless layer the next categorical color;
// an item across a cyclic origin said as its two runs. The canvas's lists and the SVG alike.

import type { Node } from '@molstar/core/repr/tree';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  DefaultRenderConfig,
  hitView,
  itemText,
  lanesHeight,
  type RangesLay,
  recordDrawLists,
  type TrackSizing,
} from '../index.js';
import { describeItems } from '../items.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';

type Rows = { start: number[]; end: number[]; id: string[]; label?: string[] };

/** One section of ranges tracks, each its rows; the section's `heights`, each track's own height when given. */
function spec(tracks: { rows: Rows; height?: number; color?: string }[], heights?: string, letters = 200): Node {
  const b = createSeqSpec({ title: 'heights' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(letters / 4) });
  b.views.view('main', (v) => {
    const sec = v.section('s', { of: s, ...(heights ? { heights } : {}) });
    tracks.forEach((t, k) => {
      const r = b.data.ranges(`r${k}`, { target: s, rows: t.rows });
      sec
        .track(`t${k}`, { label: `T${k}`, ...(t.height !== undefined ? { height: t.height } : {}) })
        .layer('l', { type: 'ranges', data: r, ...(t.color ? { color: t.color } : {}) });
    });
  });
  return b.build();
}

/** Rows that need `n` lanes (each over the one before), short ids (their labels fit inside). */
const stacked = (n: number): Rows => ({
  start: Array.from({ length: n }, (_, i) => i),
  end: Array.from({ length: n }, () => 150),
  id: Array.from({ length: n }, (_, i) => `r${i}`),
});

const Size = { width: 800, height: 600 };

/** A frame of the spec settled (its layout and lists), with the renderer's sizings and port as given. */
async function frameOf(node: Node, options: Parameters<typeof settle>[3] = {}, config = {}, runtimes = [] as never[]) {
  const s = await settle(node, runtimes, config, options);
  const f = await s.frame({ size: Size, settled: true });
  return { f, layout: f.layout, issues: f.layout.issues };
}

async function heightsOf(node: Node, options: Parameters<typeof settle>[3] = {}) {
  const r = await frameOf(node, options);
  return { r, heights: r.layout.sections[0].tracks.map((t) => t.box.height) };
}

describe('a track’s height fits its content, in steps', () => {
  it('one lane 22 px, two 36, five 78, then the lanes compact; an empty track one lane', async () => {
    expect([1, 2, 3, 4, 5, 6].map((n) => lanesHeight(n, DefaultRenderConfig))).toEqual([22, 36, 50, 64, 78, 78]);
    const { heights } = await heightsOf(spec([{ rows: stacked(1) }], 'uniform'));
    expect(heights).toEqual([22]);
    for (const [n, h] of [
      [2, 36],
      [5, 78],
      [7, 78],
    ]) {
      const one = await heightsOf(spec([{ rows: stacked(n) }]));
      expect(one.heights, `${n} lanes`).toEqual([h]);
    }
    const empty = await heightsOf(spec([{ rows: { start: [], end: [], id: [] } }]));
    expect(empty.heights).toEqual([22]);
  });

  it('tracks of one kind share the tallest; `uniform` and `fixed`; a track’s own height overrides', async () => {
    const tracks = [{ rows: stacked(1) }, { rows: stacked(3) }, { rows: stacked(2), height: 30 }];
    expect((await heightsOf(spec(tracks))).heights).toEqual([50, 50, 30]);
    expect((await heightsOf(spec(tracks, 'uniform'))).heights).toEqual([50, 50, 30]);
    expect((await heightsOf(spec(tracks, 'fixed'))).heights).toEqual([36, 36, 30]);
  });

  it('marks centred in the band by the whole data’s lanes; a pan never changes a track’s height', async () => {
    const node = spec([{ rows: { start: [0, 10, 160], end: [100, 150, 190], id: ['a', 'b', 'c'] } }]);
    const s = await settle(node);
    const at = async (viewport: { start: number; end: number }) => {
      const f = await s.frame({ size: Size, viewport });
      const t = f.layout.sections[0].tracks[0];
      return { h: t.box.height, lay: t.layers[0].lay as RangesLay, y: t.box.y };
    };
    const fit = await at({ start: 0, end: 200 });
    // two lanes over the whole data: 36 px, the lanes 4 px below the band's top
    expect(fit.h).toBe(36);
    expect(fit.lay.laneY).toBe(fit.y + 4);
    // a window where one lane is in view: the same height, the same lanes
    for (const viewport of [
      { start: 150, end: 200 },
      { start: 120, end: 170 },
    ]) {
      const w = await at(viewport);
      expect(w.h).toBe(36);
      expect(w.lay.laneY).toBe(fit.lay.laneY);
    }
  });

  it('an `ext:` sizing from the registry; one not registered is said, and the tracks fit their content', async () => {
    const half: TrackSizing = { name: 'ext:half', heights: (tracks) => tracks.map((t) => (t.content ?? t.fixed) / 2) };
    const tracks = [{ rows: stacked(1) }, { rows: stacked(3) }];
    const reg = await heightsOf(spec(tracks, 'ext:half'), { sizings: { 'ext:half': half } });
    expect(reg.heights).toEqual([11, 25]);
    const not = await heightsOf(spec(tracks, 'ext:half'));
    expect(not.heights).toEqual([50, 50]);
    expect(not.r.issues.map((i) => i.code)).toContain('seq.render.heights');
  });

  it('the host’s `trackSizing` port over all: given each track’s content and the named answer', async () => {
    const seen: number[][] = [];
    const { heights } = await heightsOf(spec([{ rows: stacked(1) }, { rows: stacked(3) }]), {
      trackSizing: (_section: unknown, tracks: readonly { content?: number }[], named: readonly number[]) => {
        seen.push(
          tracks.map((t) => t.content ?? -1),
          [...named],
        );
        return named.map(() => 40);
      },
    });
    expect(seen[0]).toEqual([22, 50]);
    expect(seen[1]).toEqual([50, 50]);
    expect(heights).toEqual([40, 40]);
    // an answer of another length is ignored
    const bad = await heightsOf(spec([{ rows: stacked(1) }]), { trackSizing: () => [1, 2] });
    expect(bad.heights).toEqual([22]);
  });

  it('the language checks the param: the three names and `ext:` ones', async () => {
    const { checkSeqSpec } = await import('@molstar/seqx-core/lang');
    expect(checkSeqSpec(spec([{ rows: stacked(1) }], 'ext:mine'))).toEqual([]);
    expect(checkSeqSpec(spec([{ rows: stacked(1) }], 'tall')).length).toBeGreaterThan(0);
  });
});

/** Two labeled items: a long one whose label fits inside, a short one whose label goes beside. */
const Labeled: Rows = {
  start: [0, 120, 124],
  end: [100, 122, 126],
  id: ['long', 'short', 'next'],
  label: ['Long domain', 'Short motif', 'Next motif'],
};

describe('labels inside if they fit, else beside', () => {
  it('a long item’s label inside its mark, a short one’s in the label lane above, the next moved with a leader', async () => {
    const node = spec([{ rows: Labeled }]);
    const r = await frameOf(node);
    const svg = await renderSvg(node, { size: Size });
    const t = r.layout.sections[0].tracks[0];
    const lay = t.layers[0].lay as RangesLay;
    // one lane and a label lane: 36 px
    expect(t.box.height).toBe(36);
    expect(lay.labelY).toBe(t.box.y + 4);
    expect(lay.laneY).toBe(lay.labelY + 14);
    // inside: the long one, whole
    const inside = Array.from({ length: lay.m }, (_, j) => j).filter((j) => lay.mlabel[j] > 0);
    expect(inside.map((j) => lay.keys[lay.mitem[j]])).toEqual(['long']);
    // beside: the two short ones, the longer first at its start, the other moved past it
    const beside = Array.from({ length: lay.ln }, (_, k) => ({
      key: lay.keys[lay.litem[k]],
      x: lay.lx[k],
      want: lay.lwant[k],
      w: lay.lw[k],
    }));
    expect(beside.map((b) => b.key).sort()).toEqual(['next', 'short']);
    const [a, b] = [...beside].sort((p, q) => p.x - q.x);
    expect(a.x).toBe(a.want);
    expect(b.x).toBeGreaterThanOrEqual(a.x + a.w);
    // the canvas's lists and the SVG alike: each label once, whole; a leader for the moved one
    for (const text of ['Long domain', 'Short motif', 'Next motif'])
      expect(svg.svg.split(`>${text}<`).length - 1).toBe(1);
    expect(recordDrawLists(svg.lists, ['labels'])).toBe(recordDrawLists(r.f.lists, ['labels']));
    const labels = recordDrawLists(r.f.lists, ['labels']);
    expect(labels.split('\n').filter((l) => l.startsWith('path ')).length).toBe(1);
    // a placed label's hit is its item's
    for (const p of beside) {
      const hit = hitView(r.layout, p.x + p.w / 2, lay.labelY + 7);
      expect(hit?.element?.key).toBe(p.key);
    }
  });

  it('labels that all fit inside need no label lane; one the lane cannot hold is dropped, never cut', async () => {
    const fits = await frameOf(spec([{ rows: { start: [0], end: [100], id: ['a'], label: ['A'] } }]));
    expect(fits.layout.sections[0].tracks[0].box.height).toBe(22);
    // many short items in a row: their labels cannot all stand beside
    const n = 40;
    const rows: Rows = {
      start: Array.from({ length: n }, (_, i) => 5 * i),
      end: Array.from({ length: n }, (_, i) => 5 * i + 2),
      id: Array.from({ length: n }, (_, i) => `i${i}`),
      label: Array.from({ length: n }, (_, i) => `A long label ${i}`),
    };
    const r = await frameOf(spec([{ rows }]));
    const svg = await renderSvg(spec([{ rows }]), { size: Size });
    const lay = r.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    expect(lay.ln).toBeGreaterThan(0);
    expect(lay.ln).toBeLessThan(n);
    // placed ones never overlap, and each drawn whole
    const xs = Array.from({ length: lay.ln }, (_, k) => [lay.lx[k], lay.lx[k] + lay.lw[k]]).sort((p, q) => p[0] - q[0]);
    for (let k = 1; k < xs.length; k++) expect(xs[k][0]).toBeGreaterThanOrEqual(xs[k - 1][1]);
    expect(svg.svg).not.toContain('…');
  });
});

describe('the label column, the gaps', () => {
  it('measured to the widest label, between 80 and 160 px, an ellipsis past it; fixed when a width is given', async () => {
    const one = (label: string) => {
      const b = createSeqSpec({ title: 'labels' });
      const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
      b.views.view('main', (v) => v.section('s', { of: s }).track('t', { label }).layer('l', { type: 'letters' }));
      return b.build();
    };
    const x = async (label: string, config = {}) => {
      const r = await frameOf(one(label), {}, config);
      return { x: r.layout.sections[0].tracks[0].box.x, text: r.layout.sections[0].tracks[0].labelText };
    };
    const { inset, labelGap } = DefaultRenderConfig;
    expect((await x('p53')).x).toBe(inset + 80 + labelGap);
    const mid = await x('Secondary structure');
    expect(mid.x).toBeGreaterThan(inset + 80 + labelGap);
    expect(mid.x).toBeLessThan(inset + 160 + labelGap);
    expect(mid.text).toBe('Secondary structure');
    const long = await x('A track label far longer than the column could ever hold');
    expect(long.x).toBe(inset + 160 + labelGap);
    expect(long.text.endsWith('…')).toBe(true);
    expect((await x('p53', { labelWidth: 112 })).x).toBe(inset + 112 + labelGap);
  });

  it('sections 20 px apart; text lines 10 apart, their bands 3', async () => {
    expect(DefaultRenderConfig.sectionGap).toBe(20);
    expect(DefaultRenderConfig.text.lineGap).toBe(10);
    expect(DefaultRenderConfig.text.bandGap).toBe(3);
    const b = createSeqSpec({ title: 'two' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT' });
    b.views.view('main', (v) => {
      v.section('a', { of: s }).track('t').layer('l', { type: 'letters' });
      v.section('b', { of: s }).track('t').layer('l', { type: 'letters' });
    });
    const r = await frameOf(b.build());
    const [a, c] = r.layout.sections;
    expect(c.tracks[0].well.y - (a.tracks[0].well.y + a.tracks[0].well.height)).toBe(20);
  });
});

describe('a layer with no color takes the next categorical color', () => {
  it('by layer order: two colourless layers differ; a spec’s color wins; one data shown twice one color', async () => {
    const r = await frameOf(spec([{ rows: stacked(1) }, { rows: stacked(1) }, { rows: stacked(1), color: '#ff0000' }]));
    const fill = (k: number) => (r.layout.sections[0].tracks[k].layers[0].lay as RangesLay).fills[0];
    const cat = (n: number) => r.layout.tokens.colors[`seq-cat-${n}` as 'seq-cat-0'];
    expect([fill(0), fill(1)]).toEqual([cat(0), cat(1)]);
    expect(fill(2)).toBe(0xff0000ff);
    // the same data in two sections: one color
    const b = createSeqSpec({ title: 'twice' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(10) });
    const d = b.data.ranges('d', { target: s, rows: stacked(1) });
    b.views.view('main', (v) => {
      v.section('a', { of: s }).track('t').layer('l', { type: 'ranges', data: d });
      v.section('b', { of: s }).track('t').layer('l', { type: 'ranges', data: d });
    });
    const two = await frameOf(b.build());
    const fills = two.layout.sections.map((x) => (x.tracks[0].layers[0].lay as RangesLay).fills[0]);
    expect(fills).toEqual([cat(0), cat(0)]);
  });
});

describe('an item across a cyclic origin', () => {
  it('says its two runs from its start, not the circle’s hull', async () => {
    const gff = ['##gff-version 3', 'ring\tt\tCDS\t91\t110\t.\t+\t.\tID=across'].join('\n');
    const b = createSeqSpec();
    const ring = b.data.sequence('ring', { alphabet: 'dna', letters: 'ACGT'.repeat(25), cyclic: true, name: 'ring' });
    const file = b.data.source('file', { format: 'gff3', value: gff });
    const f = b.data.ranges('features', { target: ring, source: file });
    b.views.view('main', (v) => v.section('s', { of: ring }).track('t').layer('l', { type: 'ranges', data: f }));
    const s = await settle(b.build(), [RecordsRuntime]);
    const frame = await s.frame({ size: Size, settled: true });
    const lay = frame.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    const [item] = describeItems(lay);
    expect(item.runs).toEqual([90, 100, 0, 10]);
    expect(itemText('ring', item)).toBe('across, ring 91–100, 1–10');
  });
});
