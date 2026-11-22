/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// No track is drawn as an empty well (_spec/README.md "Render", Tracks layout, Empty tracks): every example view, laid
// out settled as its story composes it, draws no track of ranges whose layers hold no mark at more than a one-line
// row's height. A track collapsed to its "no items" row counts too: a view that means one says so in `ExpectEmpty`, so a fixture gap
// (a declared track with no rows) fails here.
// Medium lane: it opens every view, headless on a fake clock.

// @module-tag medium

import { type RangesLay } from '@molstar/seqx-core/render';
import { ExampleViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

/** The tracks (`<view id>:<section id>/<track id>`) a story means to leave empty, with why. */
const ExpectEmpty: Readonly<Record<string, string>> = {};

/** The views the sweep leaves out, each with why (the sources are read settled: their opens and region reads awaited). */
const Skipped: ReadonlyMap<string, string> = new Map([
  [
    'vcf-zoom-in',
    'its variants come after a zoom in (the story’s point): the band says so until then, its regions refused',
  ],
]);

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** The tracks of a settled frame that are a blank well: ranges layers all holding no mark, taller than a one-line row. */
function blankWells(v: SeqViewer, id: string): string[] {
  const out: string[] = [];
  const row = v.frame!.layout.config.track.emptyPx;
  for (const s of v.frame!.layout.sections) {
    if (s.layout.name !== 'tracks') continue;
    for (const [j, t] of s.tracks.entries()) {
      // a track collapsed to its note is no well, but a story that draws one means it (`ExpectEmpty`): a fixture gap
      if (s.empties[j]) {
        out.push(`${id}:${t.track.address}`);
        continue;
      }
      if (t.box.height <= row) continue;
      const lays = t.layers.filter((l) => l.type.name === 'ranges');
      if (lays.length === 0 || lays.length !== t.layers.length) continue;
      // marks, or a density strip's bins (a source drawn at its summary level) holding counts
      const marks = lays.reduce((n, l) => {
        const y = l.lay as RangesLay;
        let bins = 0;
        for (let k = 0; k < y.strip.counts.length; k++) bins += y.strip.counts[k];
        return n + (y.m ?? 0) + (y.items ?? 0) + bins;
      }, 0);
      if (marks === 0)
        out.push(
          `${id}:${t.track.address} h=${t.box.height} ${lays
            .map((l) => {
              const y = l.lay as any;
              return `mode=${y.mode} m=${y.m} items=${y.items} pend=${y.pending?.length} tn=${y.tn} cnt=${y.strip?.counts?.length}`;
            })
            .join(';')}`,
        );
    }
  }
  return out.filter((a) => !(a in ExpectEmpty));
}

describe('no track of an example view is a blank well', () => {
  for (const view of ExampleViews.filter((x) => !Skipped.has(x.id))) {
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      expect(blankWells(viewer, view.id)).toEqual([]);
    });
  }
});
