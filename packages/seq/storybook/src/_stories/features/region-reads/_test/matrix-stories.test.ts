/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { matrixCellHover, recordDrawLists } from '@molstar/seqx-core/render';
import { methylationItems } from '@molstar/seqx-examples/regions/methylation';
import { GenotypeDeletion, GenotypeSamples, genotypeItems } from '@molstar/seqx-examples/regions/genotypes';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openDeletion } from '../openings.js';
import { advance, drawn, exported, openRegion } from './region-twin.js';

// Features/Region reads/Matrix, headless on a fake clock (the story's composition): samples × positions (the matrix
// layer, a long deletion solid, EC-220, EC-222), a matrix on a ring (EC-220), and methylation with values missing
// (EC-221): the coarse level at the fit, the items zoomed in, each in the canvas and the SVG, light and dark.

let viewer: SeqViewer | undefined;
afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

interface MatrixLay {
  readonly images: { readonly length: number };
  readonly pend: { readonly length: number };
  readonly level: number;
  readonly items: boolean;
  readonly states: number;
  readonly n: number;
  readonly tables: readonly unknown[];
}

const lay = (v: SeqViewer, label: string) =>
  v.frame!.layout.sections[0].tracks.find((t) => t.track.label === label)?.layers[0].lay as unknown as MatrixLay;

/** The data pass's rects as `[x, y, width, height, fill]`. */
function rects(v: SeqViewer): number[][] {
  return recordDrawLists(v.frame!.lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => l.split(' ').slice(1, 5).map(Number));
}

describe('Features/Region reads/Matrix, headless', () => {
  it('samples (EC-220): at the fit the coarse level, one image of bins × 40 samples, waiting first', async () => {
    const v = (viewer = await openRegion('matrix-samples'));
    await advance(v, 16);
    expect(lay(v, 'Genotypes').pend.length).toBeGreaterThan(0);
    await advance(v, 8000);
    const l = lay(v, 'Genotypes');
    expect(l.items).toBe(false);
    expect(l.level).toBeGreaterThan(1);
    expect(l.images.length).toBe(1);
    expect(l.pend.length).toBe(0);
    expect(l.states).toBe(40);
    expect(drawn(v, 'image')).toBe(1);
    // its card: a bin, the sample, a majority genotype with its share
    const card = matrixCellHover(l as never, 100_000, 3, (s) => s)!;
    expect(card.rows.map((r) => r[0])).toEqual(['Bin', 'Sample', 'Bins']);
    expect(card.rows[1][1]).toMatch(new RegExp(`^${GenotypeSamples[3]} · \\S+ \\(.+\\), \\d+% of the bin$`));
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
  });

  it('samples zoomed in (EC-222): a row a sample, the long deletion one solid bar a row across the SNVs it covers', async () => {
    const v = (viewer = await openRegion('matrix-samples'));
    await advance(v, 8000);
    openDeletion(v);
    await advance(v, 8000);
    const l = lay(v, 'Genotypes');
    expect(l.items).toBe(true);
    expect(l.level).toBe(1);
    // the items drawn: the SNVs in the window and the deletion, a rect a state each
    expect(l.n).toBeGreaterThan(40 * 20);
    // 60 kb of the 140 kb window across the band: the deletion a bar at least a third of it wide in every row, never
    // broken by the SNVs inside it
    const bars = rects(v).filter(([, , w]) => w > 150);
    expect(bars.length).toBe(40);
    expect(new Set(bars.map(([x]) => x)).size).toBe(1);
    expect(new Set(bars.map(([, y]) => y)).size).toBe(40);
    // hovering inside it where no SNV is: the deletion, the sample and its genotype (a carrier's `0/1`)
    const g = genotypeItems();
    const snvs = new Set(Array.from(g.start, (s, k) => (k === g.deletion ? -1 : s)));
    let at = GenotypeDeletion.start + 1000;
    while (snvs.has(at)) at++;
    const carrier = GenotypeDeletion.carriers[0];
    const card = matrixCellHover(l as never, at, carrier, (s) => s)!;
    expect(card.rows[0][1]).toMatch(/^del1 \(.+:200,001–260,000\)$/);
    expect(card.rows.map((r) => r[0])).toEqual(['Item', 'Sample']);
    expect(card.rows[1][1]).toBe(`${GenotypeSamples[carrier]} · 0/1 (heterozygous)`);
    expect(card.text).toContain(`Sample ${GenotypeSamples[carrier]} · 0/1 (heterozygous)`);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
  });

  it('on a ring (EC-220): 40 thin rings of runs, a primitive count bounded by the rings and the device columns', async () => {
    const v = (viewer = await openRegion('matrix-ring'));
    await advance(v, 8000);
    const l = v.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Genotypes')!.layers[0]
      .lay as unknown as MatrixLay;
    expect(l.states).toBe(40);
    expect(l.images.length).toBe(1);
    // an image on a ring is runs, never a picture: paths a run, at most a run a device column a ring
    expect(drawn(v, 'image')).toBe(0);
    const runs = drawn(v, 'path');
    expect(runs).toBeGreaterThan(40);
    expect(runs).toBeLessThanOrEqual(40 * 1500);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
  });

  it('methylation (EC-221): values by the field’s colours, a missing one in the missing colour, bins of their mean', async () => {
    const v = (viewer = await openRegion('matrix-methylation'));
    await advance(v, 8000);
    const l = lay(v, 'Methylation');
    expect(l.states).toBe(8);
    expect(l.images.length).toBe(1);
    expect(drawn(v, 'image')).toBe(1);
    const card = matrixCellHover(l as never, 250_000, 2, (s) => s)!;
    expect(card.rows[1][1]).toMatch(/^brain 1 · mean (0|1|\d\.\d+)/);
    // zoomed in: items, a missing value in the missing colour
    v.runCommand('address.focus', { address: '/data/chr/(range :from 245000 :to 255000)' });
    await advance(v, 8000);
    const items = lay(v, 'Methylation');
    expect(items.items).toBe(true);
    expect(items.n).toBeGreaterThan(0);
    // every state's cell of every CpG in view, the missing ones in the missing colour (never the ramp's low end, never a zero)
    const m = methylationItems();
    let cells = 0;
    let missing = 0;
    for (let k = 0; k < m.start.length; k++)
      if (m.start[k] >= 245_000 && m.start[k] < 255_000)
        for (let s = 0; s < m.states; s++) Number.isNaN(m.value![k * m.states + s]) ? missing++ : cells++;
    expect(missing).toBeGreaterThan(0);
    expect(items.n).toBe(cells + missing);
    const gone = items.tokens!.colors['seq-missing'];
    expect(Array.from(items.rc.subarray(0, items.n)).filter((c) => c === gone)).toHaveLength(missing);
    for (const dark of [false, true]) {
      const svg = await exported(v, dark);
      expect(svg.settled).toBe(true);
      expect(svg.issues).toEqual([]);
    }
    expect((await exported(v, true)).svg).not.toBe((await exported(v, false)).svg);
  });
});
