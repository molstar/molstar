/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The canvas and the SVG of one spec lay their lanes out alike (core _spec/README.md "Render", Lanes and Hits): every
// view of the examples (the stories' specs), opened as its story composes it and settled on the fake clock, its every
// item in the lane its SVG export puts it in. The canvas lays out as the data arrives (its lanes sticky), the export
// once it has all of it: at rest an item packed first (a canonical transcript) arrived late takes the packing order's
// lane, as the export's. The quick lane holds TP53's GFF3 transcripts, read a page at a time (its canonical transcript
// first in its gene's stack in both). Each story answers after its own delays (`'staggered'`).

import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import {
  AcrossViews,
  type ExampleView,
  ExampleViews,
  GenomeFormatsViews,
  RegionReadsViews,
} from '@molstar/seqx-examples';
import type { RangesLay } from '@molstar/seqx-core/render';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { GenomeViews } from '../../features/genome/views.js';
import { openExample } from './draw-checks.js';
import { lanesTwin } from './lanes-twin.js';

/** The views the sweep leaves out, each with why. */
const Skipped: ReadonlyMap<string, string> = new Map([
  ['vcf-zoom-in', 'its variants come after a zoom in (the story’s point): until then its band says so, no lane'],
]);

const Views: readonly ExampleView[] = [...ExampleViews, ...GenomeViews, ...RegionReadsViews, ...AcrossViews].filter(
  (v) => !Skipped.has(v.id),
);

let viewer: SeqViewer | undefined;

beforeAll(async () => {
  await loadGenome();
  await loadBridges();
  await import('@molstar/seqx-viewer/bridge-view');
});

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

describe('the canvas and the SVG of a settled view: each item in one lane', () => {
  it('GFF3 TP53: its transcripts, read a page at a time, in the SVG’s lanes; TP53-201 first in its gene’s stack', async () => {
    const view = GenomeFormatsViews.find((v) => v.id === 'gff-tp53')!;
    viewer = await openExample(view, { fixtures: true, delayMs: 150, wait: 5000, awaitOpens: 'opens' });
    const r = await lanesTwin(viewer, view.id);
    expect(r.wrong).toEqual([]);
    expect(r.compared).toBeGreaterThan(20);
    const tx = viewer.frame!.layout.sections[0].tracks.find((t) => t.track.label === 'Transcripts')!.layers[0]
      .lay as RangesLay;
    const at = Array.from(tx.keys).indexOf('transcript:ENST00000269305');
    // the first lane: the canonical transcript is packed first (the gene line's `marked`)
    expect(tx.lanes[at]).toBe(0);
  });

  for (const view of Views)
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 'staggered', wait: 5000, awaitOpens: 'opens' });
      const r = await lanesTwin(viewer, `${view.dataset}/${view.id}`);
      expect(r.wrong).toEqual([]);
    });
});
