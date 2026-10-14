/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ring window's heap guards' run (CG3; ringWindow.heap.ts, ring-window-drag.heap.ts): a ring and its linked linear
// section through the real dispatcher and behaviors, the ring's composer drawing. `frames`: the linear section panned
// and zoomed each step and the pointer round the rings and over the labels; `drag`: a press on the window, then moves
// over 32 points along its body. Each run loads the viewer afresh (`vi.resetModules`). `measure` is the guards' warm-up,
// settle and measurement, which their positive controls run on `leaking(step)`.

import { vi } from 'vitest';
import { FRAME_N, FRAME_SETTLE_MS, FRAME_WARM, FRAME_WARM_MS, heapScaling, settleFrame } from './heap.js';

/** A guard's run: its steps after the shared warm-up and settle, measured (`measure`), and what they drew and hit. */
export async function run(mode: 'frames' | 'drag', n: number, steps = FRAME_N) {
  const s = await setup(mode, n);
  const r = measure(s.step, steps);
  return { r, owner: s.owner(), sectors: s.sectors(), labelHits: s.labelHits() };
}

/**
 * The guards' measurement: the shared warm-up (by count and by time), then until a step stops allocating (past the
 * timed warm-up the frames' climb ran about 12,000 steps more, at 4 to 47 bytes a step: V8 13.6, the pinned Node,
 * 2026-10-13, 120,000 positions after 600 failed at 28 bytes a step), then `steps` steps and 4× measured.
 */
export function measure(step: (i: number) => void, steps = FRAME_N) {
  const t0 = performance.now();
  let i = 0;
  for (; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) step(i);
  settleFrame(() => step(i++), FRAME_SETTLE_MS);
  return heapScaling(steps, () => step(i++));
}

/** The ring and its linked linear section driven, and a step of `mode`. */
export async function setup(mode: 'frames' | 'drag', n: number) {
  vi.resetModules();
  const { driven } = await import('./input-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { ringWindowDoc } = await import('./ring-window-fixtures.js');
  const { DefaultPresets } = await import('../src/presets.js');
  const { HoverKind, KernelRenderers, Inputs, PointerType } = await import('../src/index.js');
  const ITEM = HoverKind.item;
  const { Dirty, SectionComposer } = await import('../src/draw/index.js');
  const { placedRingLabels } = await import('../src/draw/ring-labels.js');
  const registry = [{ presets: DefaultPresets.filter((p) => p.name === 'circular') }];
  const t = await driven(ringWindowDoc(n - 50, n), {
    registry,
    renderers: KernelRenderers as never,
    size: { ring: [600, 600], lin: [1000, 200] },
  });
  // The frames alone: the hover label's strings are one allocation per hovered element by design (IN2's guard).
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const l = t.ctx.layout.get('ring')?.current;
  const lin = t.ctx.viewports.get('lin') as import('../src/viewport.js').Viewport | null;
  if (!l || !lin) throw new Error('no layout');
  const composer = new SectionComposer(t.ctx, 'ring');
  const be = recording();
  composer.frame(be, Dirty.all);
  const area = lin.live.current.width - (lin.axis?.insetX ?? 0);
  // The labels' texts, a point in each (every third frame step is on one: hit as its item).
  const labels = placedRingLabels.get(l);
  const anchors = new Float64Array(2 * (labels?.n ?? 0));
  for (let k = 0; labels && k < labels.n; k++) {
    anchors[2 * k] = labels.pts[10 * k + 8] + labels.side[k] * 8;
    anchors[2 * k + 1] = labels.pts[10 * k + 9];
  }
  let labelHits = 0;
  const p = Inputs.pointer();
  p.section = t.ctx.layout.sectionIndex('ring');
  const view = new Float64Array(2);
  const set = lin as unknown as { setViewFrom(a: Float64Array, o: string): void };
  const at = (a: number, r: number) => {
    p.x = l.ringCx + r * Math.cos(a);
    p.y = l.ringCy + r * Math.sin(a);
  };
  if (mode === 'drag') {
    // The linear section at a fifth of the ring, centered at the top; a press on the window's body.
    view[0] = (5 * area) / n;
    view[1] = n * 0.4 * view[0];
    set.setViewFrom(view, 'api');
    // First the pointer over all the rings, unpressed: every branch of the window's hit (its grips, its body, off
    // it) warmed, as a pointer reaches a window in use. A branch first taken mid-measure deoptimized it (213 bytes a
    // step until it was optimized again).
    p.type = PointerType.move;
    for (let i = 0; i < FRAME_WARM; i++) {
      at((i * 0.0731) % (2 * Math.PI), l.ringR * (0.6 + 0.45 * ((i * 0.17) % 1)));
      t.ctx.input.pointer(p);
    }
    at(Math.PI / 2, l.ringR * 0.8);
    p.type = PointerType.down;
    p.button = 0;
    p.buttons = 1;
    t.ctx.input.pointer(p);
    p.type = PointerType.move;
  } else p.type = PointerType.move;
  const step = (i: number) => {
    if (mode === 'frames') {
      // The linear section panned and zoomed (its window round the ring, across the origin), the pointer round the
      // rings and out over the labels.
      view[0] = (area / n) * (2 + (i % 7));
      view[1] = ((i * 37) % n) * view[0];
      set.setViewFrom(view, 'api');
      at((i * 0.37) % (2 * Math.PI), l.ringR * (0.5 + 0.7 * ((i * 0.13) % 1)));
      if (i % 3 === 0 && anchors.length > 0) {
        const k = (i / 3) % (anchors.length >> 1);
        p.x = anchors[2 * k];
        p.y = anchors[2 * k + 1];
      }
      t.ctx.input.pointer(p);
      if (i % 3 === 0 && t.ctx.hover.live.current.kind === ITEM) labelHits++;
      composer.frame(be, Dirty.overlay | Dirty.tracks);
    } else {
      // Back and forth along the window's body, never the same point twice in a row.
      at(Math.PI / 2 + 0.3 * Math.sin(i * 0.05), l.ringR * 0.8);
      t.ctx.input.pointer(p);
      composer.frame(be, Dirty.overlay);
    }
  };
  return {
    step,
    owner: () => t.ctx.ui.drag.value?.owner ?? '',
    sectors: () => be.sectorItems.n,
    labelHits: () => labelHits,
  };
}
