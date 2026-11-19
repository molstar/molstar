/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Labels beside their marks (_spec/README.md "Render", Labels inside if they fit, else beside; layouts.md "Labels
// inside if they fit, else beside"): a leader never crosses another label or a mark (a label that cannot be placed so
// is dropped), neighbouring items of one layer with the same label are labeled once at the run's start (the hover
// still names each). A geometric check over the draw lists, the canvas's and the SVG's.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type DrawLists, hitView, type RangesLay, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { rng } from './small-helpers.js';

type Rows = { start: number[]; end: number[]; id: string[]; label?: string[] };
const Size = { width: 800, height: 400 };

function spec(layers: Rows[], length = 400): Node {
  const b = createSeqSpec({ title: 'label runs' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(length / 4) });
  b.views.view('main', (v) => {
    const t = v.section('s', { of: s }).track('t', { label: 'T' });
    layers.forEach((rows, k) =>
      t.layer(`l${k}`, { type: 'ranges', data: b.data.ranges(`r${k}`, { target: s, rows }) }),
    );
  });
  return b.build();
}

interface Box {
  x0: number;
  y0: number;
  x1: number;
  y1: number;
}

/** The leaders (a path of the labels pass), the labels' boxes and the data pass's rects. */
function geometry(lists: DrawLists) {
  const leaders: number[][] = [];
  const labels: Box[] = [];
  const marks: Box[] = [];
  for (const pass of ['data', 'labels'] as const)
    for (const line of recordDrawLists(lists, [pass]).split('\n')) {
      if (pass === 'data' && line.startsWith('rect ')) {
        const [x, y, w, h] = line.split(' ').slice(1, 5).map(Number);
        marks.push({ x0: x, y0: y, x1: x + w, y1: y + h });
      } else if (pass === 'labels' && line.startsWith('path ')) {
        const p = line.replace(/[ML]/g, ' ').split(/\s+/).filter(Boolean).slice(1).map(Number);
        leaders.push(p.slice(0, 6));
      } else if (pass === 'labels' && line.startsWith('text ')) {
        const m = /^text ".*" (\S+) (\S+) \S+ ".*" #\w+ width (\S+)/.exec(line)!;
        const [x, y, w] = [Number(m[1]), Number(m[2]), Number(m[3])];
        labels.push({ x0: x, y0: y - 8, x1: x + w, y1: y + 3 });
      }
    }
  return { leaders, labels, marks };
}

/** Whether the axis-aligned segment crosses the interior of a box. */
function hits(seg: number[], b: Box) {
  const [ax, ay, bx, by] = seg;
  const e = 0.01;
  return (
    Math.min(ax, bx) < b.x1 - e &&
    Math.max(ax, bx) > b.x0 + e &&
    Math.min(ay, by) < b.y1 - e &&
    Math.max(ay, by) > b.y0 + e
  );
}

/** The leaders' crossings of a label or a mark (a leader's own label starts 2 px past its end). */
function crossings(lists: DrawLists) {
  const { leaders, labels, marks } = geometry(lists);
  let n = 0;
  for (const p of leaders)
    for (const seg of [p.slice(0, 4), p.slice(2, 6)]) for (const b of [...labels, ...marks]) if (hits(seg, b)) n++;
  return { n, leaders: leaders.length, labels: labels.length };
}

/** A dense mixed layer: short items over each other, several lanes, a few labels repeated. */
function dense(): Rows {
  const r = rng(7);
  const n = 60;
  const start: number[] = [];
  const end: number[] = [];
  const label: string[] = [];
  for (let i = 0; i < n; i++) {
    const s = Math.floor(r() * 380);
    start.push(s);
    end.push(s + 1 + Math.floor(r() * 6));
    label.push(['site', 'motif', 'a long label', 'x'][Math.floor(r() * 4)] + (r() < 0.5 ? '' : ` ${i % 5}`));
  }
  return { start, end, label, id: start.map((_, i) => `i${i}`) };
}

describe('a leader never crosses a label or a mark', () => {
  it('a dense mixed layer: no crossing, canvas and SVG alike', async () => {
    const node = spec([dense()]);
    const s = await settle(node);
    const f = await s.frame({ size: Size, settled: true });
    const svg = await renderSvg(node, { size: Size });
    const c = crossings(f.lists);
    expect(c.labels).toBeGreaterThan(5);
    expect(c.n).toBe(0);
    expect(crossings(svg.lists).n).toBe(0);
    expect(recordDrawLists(svg.lists, ['labels'])).toBe(recordDrawLists(f.lists, ['labels']));
  });

  it('two layers of one track: the later layer’s leader never crosses the earlier layer’s marks', async () => {
    const a: Rows = { start: [100], end: [160], id: ['wide'], label: ['wide'] };
    const s = await settle(spec([a, dense()]));
    const f = await s.frame({ size: Size, settled: true });
    expect(crossings(f.lists).n).toBe(0);
  });

  it('a label that would need a leader through a neighbouring label is dropped, not drawn through it', async () => {
    const rows: Rows = {
      start: [0, 120, 124],
      end: [100, 122, 126],
      id: ['long', 'short', 'next'],
      label: ['Long domain', 'Short motif', 'Next motif'],
    };
    const s = await settle(spec([rows], 200));
    const f = await s.frame({ size: Size, settled: true });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    const placed = Array.from({ length: lay.ln }, (_, k) => lay.keys[lay.litem[k]]);
    expect(placed).toEqual(['short']);
    expect(crossings(f.lists).n).toBe(0);
  });
});

describe('a run of the same label is labeled once', () => {
  /** Eight neighbouring items with one label, a lone one far off, another run of a second label. */
  const rows: Rows = {
    start: [10, 12, 14, 16, 18, 20, 22, 24, 150, 300, 303, 306],
    end: [11, 13, 15, 17, 19, 21, 23, 25, 151, 301, 304, 307],
    id: Array.from({ length: 12 }, (_, i) => `i${i}`),
    label: [...Array.from({ length: 9 }, () => 'interface'), 'other', 'other', 'other'],
  };

  it('one label at the run’s start, no leader to each; far apart items keep their own', async () => {
    const node = spec([rows]);
    const s = await settle(node);
    const f = await s.frame({ size: Size, settled: true });
    const svg = await renderSvg(node, { size: Size });
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
    const placed = Array.from({ length: lay.ln }, (_, k) => lay.keys[lay.litem[k]]).sort();
    // the run's start, the lone item, the other run's start
    expect(placed).toEqual(['i0', 'i8', 'i9']);
    const labels = recordDrawLists(f.lists, ['labels']);
    expect(labels.split('text "interface"').length - 1).toBe(2);
    expect(labels.split('text "other"').length - 1).toBe(1);
    expect(svg.svg.split('>interface<').length - 1).toBe(2);
    expect(labels.split('\n').filter((l) => l.startsWith('path ')).length).toBe(0);
    expect(recordDrawLists(svg.lists, ['labels'])).toBe(labels);
    // the run's label is hit as its first item
    const k = Array.from({ length: lay.ln }, (_, q) => q).find((q) => lay.keys[lay.litem[q]] === 'i0')!;
    expect(hitView(f.layout, lay.lx[k] + 2, lay.labelY + 7)?.element?.key).toBe('i0');
  });
});
