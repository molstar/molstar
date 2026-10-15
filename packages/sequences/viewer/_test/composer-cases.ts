/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's frame cases (CS1a), for its heap guard and its quick twin (composer.heap.ts, composer.quick.ts), on a
// recording backend that reports `Copy` and `Back`: pans and scrolls through the copy at a zoom (every type), every
// pass through a pan with a navigator's window lit, a wrapped section's scroll frames under a group label mode, pans
// with range labels sticking at the view's edges, and whole frames and pans over 400,000 positions zoomed out. Each
// returns its frame and, where the copy is the point, the copies made so far. Each call loads the viewer afresh
// (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/** The zooms (CSS px a cell) of the copy's pans, the navigators, the group label modes and the range labels' zooms. */
export const COPY_CELLS = [13, 3, 0.4] as const;
export const NAVIGATORS = ['band', 'rulerRail'] as const;
export const LABEL_MODES = ['aboveEnd', 'inlineStart', 'none'] as const;
export const LABEL_CELLS = [13, 4] as const;
/** The zooms over 400,000 positions (`test:heap-deep` adds 0.03). */
export const ZOOMED_OUT_CELLS = [0.0019, 0.4] as const;
/**
 * The quick twin's counts: a frame is 0.1 to 0.7 ms (M4), so windows of 250 frames.
 */
export const FRAME: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };
/**
 * The navigators' pans warm longer: under the quick flags the band's still compiled and allocated now and then (1 to
 * 3.5 bytes a frame, 15 of code) to about 7,000 frames, two quiet windows among them (2026-10-14, Node 24, M4).
 */
export const NAVIGATOR: AllocOptions = { warmSteps: 2_000, steps: 250, windows: 3 };
/**
 * The quick counts over 400,000 positions, by zoom: a frame (a whole redraw every third) is about 4.5 ms at 0.0019
 * (every position's class binned), so windows of 100 frames there; 0.7 ms at 0.4 (M4).
 */
export const ZOOMED_OUT_COUNTS: Record<(typeof ZOOMED_OUT_CELLS)[number], AllocOptions> = {
  0.0019: { warmSteps: 500, steps: 100, windows: 3, maxWarmSteps: 3_000 },
  0.4: FRAME,
};

async function load() {
  vi.resetModules();
  const fx = await import('./composer-fixtures.js');
  const docs = await import('./prepare-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  return { ...fx, ...docs, createSeqViewSpec };
}

/** Pans and scrolls through the copy at `cell`, every type: right, left, down, up. */
export async function copyPans(cell: number) {
  const { composed, allTypes, recording, COPY_BACK } = await load();
  const t = await composed(allTypes({ length: 2000 }), { width: 1000, height: 300, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(cell, 200, 'api');
  t.frame(be);
  const frame = (i: number) => {
    const k = i & 3;
    t.vp?.panBy(k === 0 ? 7 : k === 1 ? -7 : 0, k === 2 ? 11 : k === 3 ? -11 : 0, 'pointer');
    t.frame(be);
  };
  return { frame, copies: () => t.composer.stats.copies };
}

/** Every pass through a pan, `nav`'s window lit. */
export async function navigatorPans(nav: (typeof NAVIGATORS)[number]) {
  const { composed, allTypes, recording, COPY_BACK } = await load();
  const preset = { name: `lab:${nav}`, geometry: 'linear', navigator: nav } as const;
  const t = await composed(allTypes({ length: 2000 }), {
    width: 1000,
    height: 300,
    dpr: 2,
    registry: [{ presets: [preset] }],
  });
  t.ctx.ui.setStyle(t.section, preset.name);
  await t.ctx.graph.update();
  const be = recording(COPY_BACK);
  t.vp?.setView(3, 200, 'api');
  t.frame(be, 7);
  const frame = (i: number) => {
    t.vp?.panBy((i & 1) === 0 ? 7 : -7, 0, 'pointer');
    t.frame(be, 7);
  };
  return { frame };
}

/** A wrapped section's scroll frames (letters, emphasis, the gutter's numbers) with `mode` group labels. */
export async function wrappedScroll(mode: (typeof LABEL_MODES)[number]) {
  const { composed, recording, COPY_BACK, createSeqViewSpec } = await load();
  const b = createSeqViewSpec();
  const a = b.data
    .sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(150) })
    .numbering({ name: 'auth', method: 'offset', first: 1001 });
  const starts = Array.from({ length: 60 }, (_, k) => k * 50);
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: { start: starts, stop: starts.map((s) => s + 30) },
  });
  const v = b.view({ ref: 'v', label: 'V', layout: 'text', numbering: 'auth', height: 'box' });
  v.block({ sequence: a });
  const track = v.track({ label: 'L' });
  track.letters();
  track.emphasis({ data: f });
  const t = await composed(b.getDocument(), {
    width: 900,
    height: 400,
    dpr: 2,
    render: { wrappedGroupLabels: mode },
  });
  const be = recording(COPY_BACK);
  t.vp?.scrollTo(0, 2000, 'api');
  t.frame(be);
  const frame = (i: number) => {
    t.vp?.panBy(0, (i & 1) === 0 ? 37 : -37, 'pointer');
    t.frame(be);
  };
  return { frame };
}

/** Pans out and back over 160 px with range labels centered, sticking at the view's edges and cut, at `cell`. */
export async function rangeLabels(cell: number) {
  const { composed, recording, COPY_BACK, createSeqViewSpec } = await load();
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'MKVLATGHEKWRSPQIDNCY'.repeat(30) });
  const r = b.data.annotation({
    ref: 'r',
    geometry: 'ranges',
    target: a,
    value: {
      start: [20, 250, 300, 330, 420],
      stop: [260, 290, 380, 336, 590],
      label: ['Kinase domain', 'Linker', 'SH2 domain', 'Site', 'C-terminal tail'],
    },
  });
  const v = b.view({ ref: 'v', label: 'V' });
  v.block({ sequence: a });
  v.track({ label: 'Features' }).ranges({ data: r, labels: true });
  const t = await composed(b.getDocument(), { width: 1000, height: 200, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(cell, cell * 200, 'api');
  t.frame(be);
  // Out and back over 160 px: labels move with their ranges, stick, and are cut as their ranges leave the view.
  const frame = (i: number) => {
    t.vp?.panBy((i & 32) === 0 ? 5 : -5, 0, 'pointer');
    t.frame(be);
  };
  return { frame, copies: () => t.composer.stats.copies };
}

/** Whole frames (a redraw every third) and pans over 400,000 positions at `cell`. */
export async function zoomedOut(cell: number) {
  const { composed, scaleDoc, recording, COPY_BACK } = await load();
  const t = await composed(scaleDoc(), { width: 900, height: 420, dpr: 2 });
  const be = recording(COPY_BACK);
  t.vp?.setView(cell, 0, 'api');
  t.frame(be);
  const frame = (i: number) => {
    if (i % 3 === 0) t.composer.redraw();
    else t.vp?.panBy(i % 2 === 0 ? 37 : -37, 0, 'pointer');
    t.frame(be);
  };
  return { frame };
}
