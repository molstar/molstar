/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover contract (_spec/README.md "Render", Emphasis, The hover contract), swept: every view of the examples (the
// stories' specs), opened as its story composes it, settled on the fake clock; for every layer instance, points on the
// marks its data pass drew (`probesOf`: read from the draw lists, never guessed); the pointer moved there answers with
// a hit, a hover mark in the hovered section and in every section the hover is carried to with positions in view, drawn
// by both backends (the canvas's and the SVG's), and a card with words. A rail's grab zone over a ruler is the rail's
// (its hover is its own, Navigators): counted, not probed. A new layer type, layout or view is swept by adding nothing.
// The quick lane holds one view a layout (`hover-layouts.test.ts`).

// @module-tag medium

import { appendFileSync } from 'node:fs';
import { hashString } from '@molstar/core/data/util/hash-functions';
import { BuiltInLayerTypes, BuiltInLayouts } from '@molstar/seqx-core/render/built-ins';
import { loadBridges } from '@molstar/seqx-core/render/bridge-load';
import { loadGenome } from '@molstar/seqx-core/resolve';
import { AcrossViews, type ExampleView, ExampleViews, RegionReadsViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, beforeAll, describe, expect, it, vi } from 'vitest';
import { GenomeViews } from '../../features/genome/views.js';
import { CardDigests } from './card-digests.js';
import { openExample } from './draw-checks.js';
import { type Swept, sweepHover, sweepLegends } from './hover-sweep.js';

/** The views the sweep leaves out, each with why. */
const Skipped: ReadonlyMap<string, string> = new Map([
  ['vcf-zoom-in', 'its variants come after a zoom in (the story’s point): until then its band says so, no mark'],
]);

/** The built-in layer types no story draws, each with where its hover is held instead. */
const NoStory: ReadonlyMap<string, string> = new Map([
  ['frame', 'no story outlines a region: its outline and title hit as its item (core `text-layout.test.ts`)'],
]);

const Views: readonly ExampleView[] = [...ExampleViews, ...GenomeViews, ...RegionReadsViews, ...AcrossViews].filter(
  (v) => !Skipped.has(v.id),
);

let viewer: SeqViewer | undefined;

/** A view's cards as `card-digests.ts` holds them: their `hashString`, unsigned, in hex. */
const cardsDigest = (cards: string) => (hashString(cards) >>> 0).toString(16).padStart(8, '0');

// the genome and bridge chunks landed first: their imports are real I/O, never on the fake clock
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

describe('the hover contract: every drawn mark answers a hover with a hit, a mark here and linked, and a card', () => {
  const all: Swept[] = [];
  for (const view of Views)
    it(`${view.dataset}/${view.id}`, { tags: ['medium'] }, async () => {
      viewer = await openExample(view, { fixtures: true, gene: '', delayMs: 0, wait: 5000, awaitOpens: true });
      const id = `${view.dataset}/${view.id}`;
      const r = sweepHover(viewer, id);
      all.push(r);
      // the cards as composed by the providers: today's, byte for byte (`card-digests.ts`; a dump refreshes them)
      const dump = process.env.SEQ_CARDS_DUMP;
      if (dump) appendFileSync(dump, `${JSON.stringify([id, r.cards])}\n`);
      else if (CardDigests[id] !== undefined) expect(cardsDigest(r.cards), id).toBe(CardDigests[id]);
      expect(r.wrong).toEqual([]);
      // a view whose hover reaches another section checks it there: a projection dropped fails here
      if (r.links) expect(r.linked).toBeGreaterThan(0);
      expect(sweepLegends(viewer, `${view.dataset}/${view.id}`, r.entries)).toEqual([]);
    });

  it('reaches every built-in layer type and layout, and the bridges', { tags: ['medium'] }, () => {
    const types = new Set([...all.flatMap((r) => r.types), ...NoStory.keys()]);
    const layouts = new Set(all.flatMap((r) => r.layouts));
    expect(Object.keys(BuiltInLayerTypes).filter((t) => !types.has(t))).toEqual([]);
    expect(Object.keys(BuiltInLayouts).filter((l) => !layouts.has(l))).toEqual([]);
    for (const l of ['blocks', 'rows', 'bridge']) expect(layouts.has(l), l).toBe(true);
    // the sweep is a sample: a few marks a layer, still thousands of probes
    expect(all.reduce((n, r) => n + r.probes, 0)).toBeGreaterThan(2000);
  });
});
