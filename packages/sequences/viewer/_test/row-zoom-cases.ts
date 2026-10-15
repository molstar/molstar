/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The row zoom's cases (KS2), for its heap guard and its quick twin (row-zoom.heap.ts, row-zoom.quick.ts), each on a
// viewer loaded afresh (`vi.resetModules`): the Alt wheel over case 3's rows, its layout's hit tests counted (`altWheel`);
// the Alt wheel over a dense genotype matrix, whose track follows the row scale (`matrixWheel`); the vertical
// overview's bottom grip dragged (`bottomGrip`). Each returns its step and the layouts' build and rescale counts. The
// wheel input is written inline by each step: a double passed to a helper is boxed, in the test too.

import { vi } from 'vitest';

const S = 'family.main';

async function open() {
  vi.resetModules();
  const { case3 } = await import('../../core/_test/format/cases.js');
  const { driven } = await import('./input-fixtures.js');
  const api = await import('../src/index.js');
  const { WindowKind } = await import('../src/windows.js');
  const t = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
  return { t, api, WindowKind };
}

/**
 * The Alt wheel over case 3's rows, below their full height both ways (so the vertical overview stays: a change there
 * rebuilds, by design), the pointer resting over them; `hits` counts the layout's `hitTestAt` calls.
 */
export async function altWheel() {
  const { t, api } = await open();
  t.vp(S).setRowView(0.5, 0, 'api');
  const y = t.layout(S).scrolledTop + 9.5;
  const counter = { hits: 0, sink: 0 };
  const layout = t.ctx.layout;
  const hitTestAt = layout.hitTestAt.bind(layout);
  layout.hitTestAt = (section, out) => {
    counter.hits++;
    return hitTestAt(section, out);
  };
  t.pointer('move', S, 400.5, y);
  const w = api.Inputs.wheel();
  w.section = t.ctx.layout.sectionIndex(S);
  w.x = 400.5;
  w.y = y;
  w.deltaX = 0;
  w.modifiers = api.Modifier.alt;
  const step = (i: number) => {
    w.deltaY = i & 1 ? -4.5 : 4.5;
    counter.sink = (counter.sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
  };
  return { t, step, counter, stats: t.layouts.stats };
}

/** A dense genotype matrix above letters tracks, in a box it overflows, the Alt wheel over it; `heights`: the track's. */
export async function matrixWheel() {
  vi.resetModules();
  const { genotypeDoc } = await import('./genotype-fixtures.js');
  const { driven } = await import('./input-fixtures.js');
  const api = await import('../src/index.js');
  const doc = genotypeDoc({ samples: 3000, positions: 200, below: 6 });
  // The kernel's renderers: the matrix's own height, which follows the row scale.
  const renderers = [...api.KernelRenderers] as never;
  const t = await driven(doc, { size: { 'v.main': [1000, 200] }, labels: false, renderers });
  const M = 'v.main';
  t.vp(M).setRowView(1.5, 0, 'api');
  const y = t.layout(M).scrolledTop + 30.5;
  t.pointer('move', M, 400.5, y);
  const w = api.Inputs.wheel();
  w.section = t.ctx.layout.sectionIndex(M);
  w.x = 400.5;
  w.y = y;
  w.deltaX = 0;
  w.modifiers = api.Modifier.alt;
  let sink = 0;
  const step = (i: number) => {
    w.deltaY = i & 1 ? -4.5 : 4.5;
    sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
  };
  // The matrix's track (the first scrolled entry) follows the zoom: its heights over four calls.
  const heights = new Set<number>();
  for (let i = 0; i < 4; i++) {
    step(i);
    heights.add(t.layout(M).entryHeight[t.layout(M).pinnedCount]);
  }
  return { step, heights: heights.size, stats: t.layouts.stats, sink: () => sink };
}

/** The vertical overview's bottom grip, pressed and moved by a step (one `PointerInput` reused, as an engine's). */
export async function bottomGrip() {
  const { t, api, WindowKind } = await open();
  t.vp(S).zoomTo(0, 0, 22);
  const o = t.frame(S, WindowKind.overview);
  const ox = (o[8] + o[9]) / 2;
  const idle = t.ctx.ui.drag.value === null;
  t.pointer('down', S, ox, o[11] - 1);
  const e = api.Inputs.pointer();
  e.section = t.ctx.layout.sectionIndex(S);
  e.buttons = 1;
  e.button = -1;
  e.type = api.PointerType.move;
  e.x = ox;
  const step = (i: number) => {
    e.y = o[11] - 1 + ((i & 31) - 8) * 0.4;
    t.ctx.input.pointer(e);
  };
  step(0);
  step(1);
  return {
    step,
    idle,
    owner: t.ctx.ui.drag.value?.owner,
    stats: t.layouts.stats,
    release: () => t.pointer('up', S, ox, o[11]),
  };
}
