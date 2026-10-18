/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The line styles' cases, for their heap guard and its quick twin (line-styles.heap.ts, line-styles.quick.ts): a ranges
// layer whose items cycle through every line style (`field_line_style`), drawn whole each step by the section's
// composer. `tracks`: a linear section zoomed in, panned left and right a step each (dashed and dotted strokes, waves,
// doubles, boxes, highlights over the rows); `text`: the text view scrolled down and back (the styled underline lanes,
// and the highlights and boxes behind each line's letters). Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/** The quick counts: a frame redraws the section (well under a millisecond), so windows of 250 steps. */
export const STYLED: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

const STYLES = ['solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight'];

export async function lineStyles(mode: 'tracks' | 'text') {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { createSeqViewSpec } = await import('@molstar/seq-core/builder');
  const { KernelRenderers } = await import('../src/index.js');
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const S = 'v.main';
  const n = 2000;
  const b = createSeqViewSpec();
  const a = b.data.sequence({ ref: 'a', value: 'ACDEFGHIKL'.repeat(n / 10) });
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: a,
    value: Array.from({ length: 200 }, (_, i) => ({
      start: i * 10,
      stop: i * 10 + 7 + (i % 3) * 5,
      label: `f${i}`,
      line_style: STYLES[i % STYLES.length],
    })),
  } as never);
  const v = b.view({ ref: 'v', label: 'V', ...(mode === 'text' ? { layout: 'text', height: 'box' } : {}) } as never);
  v.block({ sequence: a });
  v.track({ ref: 'v.letters' }).letters();
  v.track({ ref: 'v.features' }).ranges({
    data: f,
    labels: true,
    line_style: { theme: 'field_line_style' },
  } as never);
  const t = await driven(b.getDocument(), {
    size: { [S]: [600, 800] },
    renderers: KernelRenderers as never,
    labels: false,
  });
  await t.ctx.graph.update();
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const vp = t.vp(S);
  // Linear: 3 px a position, about 20 items in view.
  if (mode === 'tracks') vp.setView(3, 400, 'api');
  await t.ctx.graph.update();
  const composer = new SectionComposer(t.ctx, S);
  const be = recording();
  composer.frame(be, Dirty.all);
  const step = (i: number) => {
    be.reset();
    // Moved a step each way (a frame with nothing moved draws nothing): the section redrawn whole.
    if (mode === 'tracks') vp.panBy(i & 1 ? -11 : 11, 0, 'pointer');
    else vp.panBy(0, i & 1 ? -11 : 11, 'pointer');
    composer.frame(be, Dirty.tracks);
  };
  const state = () => ({ rects: be.rectItems.n, paths: be.pathItems.n });
  return { step, state };
}
