/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Values look story, headless on a fake clock, composed as the story composes it (`compose.ts`): each cell draws, at
// the range its title says (zoomed in: bars; out: a profile), a missing value is a gap, and hover outlines a value.

import { readFileSync } from 'node:fs';
import { type ValuesLay, recordDrawLists } from '@molstar/seqx-core/render';
import { PackageFolder } from '@molstar/seqx-examples/node-assets';
import { DefaultSeqViewerSpec } from '@molstar/seqx-viewer/default-spec';
import { createSeqViewer, type SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { withViewRegistry } from '../../../examples/view-registry.js';
import { AllCells, Cells, ColorCells, zoomed } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const bytes = (uri: string) => new Uint8Array(readFileSync(`${PackageFolder}${uri}`));
const reader = {
  readText: async (uri: string) => new TextDecoder().decode(bytes(uri)),
  readBytes: async (uri: string) => bytes(uri),
  openRange: async (uri: string) => {
    const file = bytes(uri);
    return {
      size: file.length,
      read: async (ranges: readonly { offset: number; length: number }[]) =>
        ranges.map((r) => file.subarray(r.offset, r.offset + r.length)),
      dispose() {},
    };
  },
};

async function open(id: string, firstWidth = 480) {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
  const cell = AllCells.find((c) => c.id === id)!;
  const base = DefaultSeqViewerSpec();
  const v = createSeqViewer({
    ...base,
    ui: undefined,
    registry: withViewRegistry(base.registry ?? [], cell.view.registry?.({ delayMs: 100 }) ?? []),
    ports: { assets: () => ({ reader: reader as never }) },
  });
  viewer = v;
  v.resize(firstWidth);
  void v.load(cell.view.build()).then(() => v.settle());
  await vi.advanceTimersByTimeAsync(3000);
  if (firstWidth !== 480) v.resize(480);
  zoomed(cell, v);
  await vi.advanceTimersByTimeAsync(3000);
  v.draw();
  return v;
}

const values = (v: SeqViewer) => {
  for (const s of v.frame!.layout.sections)
    for (const t of s.tracks) for (const l of t.layers) if (l.type.name === 'values') return l.lay as ValuesLay;
  throw new Error('no values layer');
};

/** The data pass as text: its lines per primitive. */
const dataLines = (v: SeqViewer) => recordDrawLists(v.frame!.lists, ['data']).split('\n');

/** The values layer's bars: rounded rects of the data pass (the series color at 60%, tops rounded where they have height). */
const bars = (v: SeqViewer) => dataLines(v).filter((l) => l.startsWith('rect ') && / #[0-9a-f]{6}99( |$)/.test(l));

describe('the values look story, headless', () => {
  it.each(Cells.map((c) => c.id))('%s: draws, at the range its title says, and hover outlines a value', async (id) => {
    const cell = Cells.find((c) => c.id === id)!;
    const v = await open(id);
    const lay = values(v);
    expect(lay.n, id).toBeGreaterThan(10);
    if (cell.zoom) {
      // zoomed in: the view holds the range of the title, wide enough for bars, and the layer draws them
      const r = v.view.get(cell.zoom.section)!;
      expect([r.start, r.end], id).toEqual([cell.zoom.start, cell.zoom.end]);
      expect(lay.axis.scale, id).toBeGreaterThanOrEqual(4);
      expect(bars(v).length, id).toBeGreaterThanOrEqual(cell.zoom.end - cell.zoom.start - 3);
    } else {
      // zoomed out: a profile (paths), no bar
      expect(lay.axis.scale, id).toBeLessThan(4);
      expect(dataLines(v).filter((l) => l.startsWith('path ')).length, id).toBeGreaterThan(0);
      expect(bars(v), id).toEqual([]);
    }
    const hit = v.hit(lay.box.x + lay.box.width / 2, lay.box.y + lay.box.height / 2);
    expect(hit?.layer, id).toBe(lay.key);
    v.hover.set(hit, 'test');
    const f = v.draw() ?? v.frame!;
    expect(
      recordDrawLists(f.lists, ['hover'])
        .split('\n')
        .filter((l) => l.startsWith('rect ')).length,
      id,
    ).toBeGreaterThan(0);
  });

  it('a missing value is a gap: no bar and no mark in the missing color at residues 100 and 101', async () => {
    const v = await open('conservation-out');
    v.view.setRange('/views/main/protein', 95, 105);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    // the bars' left edges: one pitch apart but across the two missing positions, three pitches
    const xs = bars(v).map((l) => Number(l.split(' ')[1]));
    expect(xs.length).toBeGreaterThan(6);
    const pitch = xs[1] - xs[0];
    const steps = xs.slice(1).map((x, i) => Math.round((x - xs[i]) / pitch));
    expect(steps.filter((k) => k !== 1)).toEqual([3]);
    const missing = (v.frame!.layout.tokens.colors['seq-missing'] >>> 0).toString(16).padStart(8, '0');
    expect(dataLines(v).filter((l) => l.includes(`#${missing}`))).toEqual([]);
  });
});

describe('a viewer mounted before it has a width', () => {
  it.each(['conservation-out', 'bigwig-out', 'bigwig-in'])('%s: draws its data once it is sized', async (id) => {
    const v = await open(id, 0);
    expect(values(v).n, id).toBeGreaterThan(10);
  });
});

/** The values layers of the viewer's view by key. */
const layers = (v: SeqViewer) => {
  const out = new Map<string, ValuesLay>();
  for (const s of v.frame!.layout.sections)
    for (const t of s.tracks)
      for (const l of t.layers) if (l.type.name === 'values') out.set(l.lay.key.split('/').pop()!, l.lay as ValuesLay);
  return out;
};

describe('the values look story: uniform against by value', () => {
  it.each(ColorCells.map((c) => c.id))(
    '%s: a layer a color and encoding, drawn at the density its title says',
    async (id) => {
      const cell = ColorCells.find((c) => c.id === id)!;
      const v = await open(id);
      const lays = layers(v);
      expect([...lays.keys()].sort()).toEqual(
        ['bars', 'line', 'points'].flatMap((e) => [`uniform-${e}`, `value-${e}`]).sort(),
      );
      for (const e of ['bars', 'line', 'points']) {
        expect(lays.get(`uniform-${e}`)!.perBar, `${id} uniform ${e}`).toBe(false);
        expect(lays.get(`value-${e}`)!.perBar, `${id} by value ${e}`).toBe(true);
      }
      const lay = lays.get('value-line')!;
      const ground = (v.frame!.layout.tokens.colors['seq-track'] >>> 0).toString(16).padStart(8, '0');
      // dots: a ring of the well's ground a value, only where a value has 6 px
      const rings = dataLines(v).filter((l) => l.startsWith('path ') && l.includes(`#${ground}`));
      if (cell.zoom) {
        expect(lay.axis.scale, id).toBeGreaterThanOrEqual(6);
        expect(rings.length, id).toBeGreaterThanOrEqual(2 * (cell.zoom.end - cell.zoom.start - 2));
      } else {
        expect(lay.axis.scale, id).toBeLessThan(6);
        expect(rings, id).toEqual([]);
      }
    },
  );
});
