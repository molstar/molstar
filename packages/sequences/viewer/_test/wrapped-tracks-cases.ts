/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tracks drawn in full in the text view as a step, for the heap guard and its quick twins (wrapped-tracks.heap.ts,
// wrapped-tracks.quick.ts; seqspec:viewer#the-text-view, seqspec:viewer#performance): on the text preset with
// `wrappedTracks: 'full'`, the axes and the labels on, a section whose tracks are drawn whole each step as the rows
// scroll down and back. `kind` selects the tracks: `all` (the guard's: a tall values line and an area, sites as glyphs,
// domains with their labels and secondary structure over every line, prepare-fixtures.ts `publicationDoc`), or the
// letters and one track alone: a values `line`, `area` or `bars` strip, `points` (sites as glyphs) or labeled `ranges`.
// Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export type TracksKind = 'all' | 'line' | 'area' | 'bars' | 'points' | 'ranges';
export const TRACKS_KINDS: readonly TracksKind[] = ['all', 'line', 'area', 'bars', 'points', 'ranges'];

/**
 * The view for `kind`. A single track has a `height` of 2 rows (a values strip) or 1 (the glyphs, the labeled ranges)
 * beside the letters, over the same 120 positions as `publicationDoc`.
 */
async function tracksDoc(kind: TracksKind) {
  const { publicationDoc } = await import('./prepare-fixtures.js');
  if (kind === 'all') return publicationDoc();
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', alphabet: 'protein', value: 'MSEEKLKELAEKAIRLAKEGNPEEVLRALQ'.repeat(4) });
  const value = Array.from({ length: 120 }, (_, i) => Math.round(50 + 40 * Math.sin(i / 6)));
  const d = b.data.annotation({ ref: 'd', geometry: 'values', target: a, value: { value } });
  const sites = b.data.annotation({ ref: 's', geometry: 'ranges', target: a, value: { start: [5, 33, 77, 101] } });
  const dom = b.data.annotation({
    ref: 'dom',
    geometry: 'ranges',
    target: a,
    value: { start: [2, 60], stop: [40, 110], label: ['Kinase', 'Linker region'] },
  } as never);
  const v = b.view({ ref: 'v', label: 'V', layout: 'text' });
  v.track({ ref: 'v.letters', label: 'Letters' }).letters();
  const t = v.track({ ref: `v.${kind}`, label: kind, height: kind === 'points' || kind === 'ranges' ? 1 : 2 });
  if (kind === 'points') t.points({ data: sites });
  else if (kind === 'ranges') t.ranges({ data: dom, labels: true } as never);
  else t.values({ data: d, encoding: kind, domain: [0, 100] });
  return b.getDocument();
}

export async function wrappedTracks(kind: TracksKind = 'all') {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { DefaultPresets, KernelRenderers } = await import('../src/index.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const S = 'v.main';
  const text = DefaultPresets.find((p) => p.name === 'wrapped');
  // A single track's lines are shorter: a shorter section, so that the rows still scroll.
  const t = await driven(await tracksDoc(kind), {
    size: { [S]: [500, kind === 'all' ? 400 : 100] },
    registry: [
      {
        presets: [
          {
            ...(text as NonNullable<typeof text>),
            name: 'lab:text',
            render: { ...text?.render, wrappedTracks: 'full', wrappedTrackAxes: true, wrappedTrackLabels: true },
          },
        ],
      },
    ],
    renderers: KernelRenderers as never,
    labels: false,
  });
  t.ctx.ui.setStyle(S, 'lab:text');
  await t.ctx.graph.update();
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const l = t.layout(S);
  if (l.lines.count < 3 || l.scrolledHeight <= l.height) throw new Error('the text does not scroll');
  const composer = new SectionComposer(t.ctx, S);
  const be = recording();
  composer.frame(be, Dirty.all);
  const vp = t.vp(S);
  let paths = 0;
  let rects = 0;
  const step = (i: number) => {
    // Scrolled down and back a step each (a frame with nothing moved draws nothing): the lines redrawn whole.
    vp.panBy(0, i & 1 ? -11 : 11, 'pointer');
    be.reset();
    composer.frame(be, Dirty.tracks);
    paths += be.pathItems.n;
    rects += be.rectItems.n;
  };
  return {
    step,
    /** What the frames drew: the paths and rects over the steps run, the last frame's glyphs and texts. */
    state: () => ({ paths, rects, glyphs: be.glyphItems.n, texts: be.textItems.n }),
    /** Counts anew (a guard counts over its measured steps alone). */
    resetCounts: () => {
      paths = rects = 0;
    },
  };
}
