/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The row zoom's heap guards (seqspec:viewer#performance; KS2): a row scale alone rescales the rows (`Layouts.rescale`:
// the rows' height and the entries below them), never rebuilding the layout, so the Alt wheel's zoom and the vertical
// overview's bottom grip are flat per event, with no build counted, and so is the Alt wheel over a dense matrix, whose
// track follows the row scale. Before KS2 each event rebuilt the layout (about
// 2,120 bytes). The positive control: the Alt wheel's calls allocating an object each (`leaking`), through the same
// warm-up (`prime`) and measurement, are caught. Each case loads the viewer afresh (`vi.resetModules`): after
// input.heap's other cases had run, the grip allocated 31 bytes a move.

import { describe, expect, it, vi } from 'vitest';
import { controlCaught, describeScaling, HEAP_DEEP, heapScaling, isFlat, leaking } from './heap.js';

const N = HEAP_DEEP ? 50_000 : 25_000;
const S = 'family.main';

/** As input.heap's: at least 300,000 steps, then until three chunks in a row allocate under 2 KB. */
function prime(step: (i: number) => void) {
  const gc = (globalThis as { gc?: () => void }).gc;
  let k = 0;
  for (let chunk = 0, quiet = 0; chunk < 100 && (quiet < 3 || chunk < 15); chunk++) {
    gc?.();
    const before = process.memoryUsage().heapUsed;
    for (let j = 0; j < 20_000; j++) step(k++);
    quiet = process.memoryUsage().heapUsed - before < 2048 ? quiet + 1 : 0;
  }
}

async function open() {
  vi.resetModules();
  const { case3 } = await import('../../core/_test/format/cases.js');
  const { driven } = await import('./input-fixtures.js');
  const api = await import('../src/index.js');
  const { WindowKind } = await import('../src/windows.js');
  const t = await driven(case3(), { size: { [S]: [1000, 150] }, labels: false });
  return { t, api, WindowKind };
}
/** A dense genotype matrix above letters tracks, in a box it overflows: its track follows the row scale. */
async function openMatrix() {
  vi.resetModules();
  const { genotypeDoc } = await import('./genotype-fixtures.js');
  const { driven } = await import('./input-fixtures.js');
  const api = await import('../src/index.js');
  const doc = genotypeDoc({ samples: 3000, positions: 200, below: 6 });
  // The kernel's renderers: the matrix's own height, which follows the row scale.
  const renderers = [...api.KernelRenderers] as never;
  const t = await driven(doc, { size: { 'v.main': [1000, 200] }, labels: false, renderers });
  return { t, api };
}

describe('the row zoom rescales the rows: no rebuild, flat', () => {
  it('the Alt wheel: two hits a call, no build; an object allocated a call is caught (positive control)', async () => {
    const { t, api } = await open();
    // Below the rows' full height both ways, so the vertical overview stays (a change there rebuilds, by design).
    t.vp(S).setRowView(0.5, 0, 'api');
    const y = t.layout(S).scrolledTop + 9.5;
    let hits = 0;
    const layout = t.ctx.layout;
    const hitTestAt = layout.hitTestAt.bind(layout);
    layout.hitTestAt = (section, out) => {
      hits++;
      return hitTestAt(section, out);
    };
    t.pointer('move', S, 400.5, y);
    // The wheel input, written inline by each step: a double passed to a helper is boxed, in the test too.
    const w = api.Inputs.wheel();
    w.section = t.ctx.layout.sectionIndex(S);
    w.x = 400.5;
    w.y = y;
    w.deltaX = 0;
    w.modifiers = api.Modifier.alt;
    let sink = 0;
    const step = (i: number) => {
      w.deltaY = i & 1 ? -4.5 : 4.5;
      sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
    };
    prime(step);
    const [builds, rescales] = [t.layouts.stats.builds, t.layouts.stats.rescales];
    hits = 0;
    const calls = 4000;
    for (let i = 0; i < calls; i++) step(i);
    // The wheel's own hit, and the pointer's again after the rescale: two a call, not three.
    expect(hits / calls).toBeLessThanOrEqual(2);
    expect(hits / calls).toBeGreaterThan(1);
    const r = heapScaling(N, step);
    expect(t.layouts.stats.builds).toBe(builds);
    expect(t.layouts.stats.rescales - rescales).toBeGreaterThan(calls + N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(
      controlCaught(() => {
        const leaky = leaking(step);
        prime(leaky);
        for (let i = 0; i < calls; i++) leaky(i);
        return heapScaling(N, leaky);
      }),
    ).toBe(true);
    expect(sink).toBeGreaterThanOrEqual(0);
  }, 180_000);

  it('the Alt wheel over a dense matrix: its track rescaled each call, no build, flat', async () => {
    const { t, api } = await openMatrix();
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
    const heights = new Set<number>();
    const step = (i: number) => {
      w.deltaY = i & 1 ? -4.5 : 4.5;
      sink = (sink + (t.ctx.input.wheel(w) ? 1 : 0)) & 0xffff;
    };
    for (let i = 0; i < 4; i++) {
      step(i);
      heights.add(t.layout(M).entryHeight[t.layout(M).pinnedCount]);
    }
    // The matrix's track (the first scrolled entry) follows the zoom.
    expect(heights.size).toBe(2);
    prime(step);
    const [builds, rescales] = [t.layouts.stats.builds, t.layouts.stats.rescales];
    const r = heapScaling(N, step);
    expect(t.layouts.stats.builds).toBe(builds);
    expect(t.layouts.stats.rescales - rescales).toBeGreaterThan(N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sink).toBeGreaterThanOrEqual(0);
  }, 180_000);

  it("the vertical overview's bottom grip: rescaled each move, no build, flat", async () => {
    const { t, api, WindowKind } = await open();
    t.vp(S).zoomTo(0, 0, 22);
    const o = t.frame(S, WindowKind.overview);
    const ox = (o[8] + o[9]) / 2;
    t.pointer('down', S, ox, o[11] - 1);
    // A pressed pointer's moves: one `PointerInput` reused, as an engine's.
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
    expect(t.ctx.ui.drag.value === null).toBe(true);
    step(0);
    step(1);
    expect(t.ctx.ui.drag.value?.owner).toBe('verticalOverview');
    prime(step);
    const [builds, rescales] = [t.layouts.stats.builds, t.layouts.stats.rescales];
    const r = heapScaling(N, step);
    expect(t.layouts.stats.builds).toBe(builds);
    expect(t.layouts.stats.rescales - rescales).toBeGreaterThan(N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    t.pointer('up', S, ox, o[11]);
  }, 180_000);
});
