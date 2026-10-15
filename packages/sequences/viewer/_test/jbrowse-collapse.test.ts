/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST4b follow-up: a gene track whose lanes do not fit draws one lane a gene (JBrowse 2's "collapse"), its transcripts'
// exons folded; the track menu's "Expand rows" shows every transcript, with the lane cap's hint. Over the demo's
// recorded slices (jbrowse-sources.ts, the TP53 view of the navigation test).

import { afterAll, expect, it } from 'vitest';
import { Dirty } from '../src/draw/composer.js';
import type { Streams } from '../src/streams.js';
import { composed, recording } from './composer-fixtures.js';
import { demoSources } from './jbrowse-sources.js';
import { JBROWSE_SECTION, jbrowseDoc } from './jbrowse-story-fixtures.js';

const src = demoSources('nav');
afterAll(src.save);
const TP53: [number, number] = [7_661_779, 7_687_538];
const CAPPED = "more lanes than fit: raise the track's height";

it("folds a gene's transcripts into a lane while the lanes do not fit, and Expand rows shows them all", async () => {
  const t = await composed(jbrowseDoc('chr17', ...TP53), { width: 800, height: 400, dpr: 1, loaders: src.loaders });
  const m = t.ctx.streams as Streams;
  t.layout(JBROWSE_SECTION);
  let done = false;
  void m.settled().then(() => {
    done = true;
  });
  for (let k = 0; k < 2000 && !done; k++) {
    t.scheduler.advance(16);
    await new Promise((ok) => setTimeout(ok, 0));
  }
  expect(done).toBe(true);
  const texts = () => {
    const be = t.frame(recording(), Dirty.tracks);
    return be.textItems.text.slice(0, be.textItems.n);
  };
  // The genes layer and its track.
  const view = t.ctx.graph.assembled.value?.view;
  if (!view) throw new Error('no view');
  let [layer, track] = ['', ''];
  for (let k = 0; k < view.resolved.layers.count; k++)
    if (view.set.strings.get(view.resolved.layers.type[k]) === 'ranges') {
      layer = view.set.strings.get(view.resolved.layers.ref[k]);
      track = view.set.strings.get(view.resolved.tracks.ref[view.resolved.layers.track[k]]);
    }
  /** The lanes the layer's regions give where the view starts (the lanes of the items the draw reads now). */
  const lanes = () => {
    const out = new Int32Array(6);
    out[0] = 0;
    return t.ctx.streams.get(layer)?.items?.(0, TP53[0], TP53[0] + 1, 0, out) ? out[5] : -1;
  };
  // Folded by default (TP53 has dozens of transcripts): a lane a gene, the gene's name once, no hint.
  const folded = texts();
  expect(folded).not.toContain(CAPPED);
  expect(folded.filter((x) => x === 'TP53').length).toBe(1);
  expect(lanes()).toBeLessThanOrEqual(3);
  // Expand rows: every transcript, as many as the band holds, and the hint.
  expect(t.ctx.actions.can('group.toggle', track)).toBe(true);
  t.ctx.actions.run('group.toggle', track);
  const expanded = texts();
  expect(expanded).toContain(CAPPED);
  expect(expanded.filter((x) => x === 'TP53').length).toBeGreaterThan(5);
  expect(lanes()).toBeGreaterThan(20);
  // And back.
  t.ctx.actions.run('group.toggle', track);
  expect(texts()).toEqual(folded);
  t.ctx.dispose();
}, 600_000);
