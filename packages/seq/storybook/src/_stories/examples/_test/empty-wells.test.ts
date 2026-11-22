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

/** The views the sweep leaves out: their layout is not the tracks layout's page of tracks, or too large to open here. */
const Skipped: ReadonlySet<string> = new Set([]);

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
    for (const t of s.tracks) {
      // a track collapsed to its note is no well, but a story that draws one means it (`ExpectEmpty`): a fixture gap
      if (t.empty) {
        out.push(`${id}:${t.track.address}`);
        continue;
      }
      if (t.box.height <= row) continue;
      const lays = t.layers.filter((l) => l.type.name === 'ranges');
      if (lays.length === 0 || lays.length !== t.layers.length) continue;
      const marks = lays.reduce((n, l) => n + ((l.lay as RangesLay).m ?? 0) + ((l.lay as RangesLay).items ?? 0), 0);
      // a density strip (a source drawn at its summary level) draws its bins, not marks
      if (lays.some((l) => (l.lay as RangesLay).mode === 'density')) continue;
      const pending = lays.some((l) => ((l.lay as RangesLay).pending?.length ?? 0) > 0);
      if (marks === 0 && !pending) out.push(`${id}:${t.track.address}`);
    }
  }
  return out.filter((a) => !(a in ExpectEmpty));
}

describe('no track of an example view is a blank well', () => {
  for (const view of ExampleViews.filter((x) => !Skipped.has(x.id))) {
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      viewer = await openExample(view, { fixtures: true, gene: '' });
      expect(blankWells(viewer, view.id)).toEqual([]);
    });
  }
});
