/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The canvas half of the retention test (IN1; CS2): a context with a mounted surface, used (drawn, panned, hovered,
// selected), then disposed with its surface's host dropped, frees the context, its SeqSet and resolved view, the
// surface, its canvases and, with the context, the per-context shared images, scratch and measurer. Nothing in the page
// (the document's listeners, a MutationObserver, the font set, a ResizeObserver, the matchMedia listener, the
// scheduler) keeps any of them. Run with the browser heap config (`gc` exposed). The positive control: the same, with
// the host keeping the context, collected and counted alike, is caught.

import { expect, it } from 'vitest';
import { mountSection } from '../../src/canvas/index.js';
import { Inputs, PointerType } from '../../src/index.js';
import { syntheticFamily } from '../families.js';
import { controlCaught } from './heap-browser.js';
import { opened } from './open.js';

const gc = (globalThis as { gc?: () => void }).gc as () => void;

/**
 * A viewer used through a mounted surface and a second one, disposed; weak references to what must go. `kept`: the
 * host keeps the context there (the control).
 */
async function used(kept?: unknown[]) {
  const t = await opened(syntheticFamily(60, 8), { width: 600, height: 260, dpr: 2 });
  const { ctx, surface, host, section } = t;
  const host2 = document.createElement('div');
  host2.style.cssText = 'width:300px;height:200px';
  document.body.appendChild(host2);
  const second = mountSection(ctx, section, host2, { dpr: 1 });
  const vp = ctx.viewports.get(section);
  vp?.setView(12, 0, 'api');
  const p = Inputs.pointer();
  p.type = PointerType.move;
  p.section = 0;
  for (let k = 0; k < 20; k++) {
    p.x = 150 + k * 7;
    p.y = 90 + k;
    ctx.input.pointer(p);
    vp?.panBy(9, 0, 'pointer');
    t.advance();
    surface.renderNow();
    second.renderNow();
  }
  const refs = {
    ctx: new WeakRef(ctx),
    set: new WeakRef(ctx.graph.assembled.value?.set as object),
    view: new WeakRef(ctx.graph.assembled.value?.view as object),
    surface: new WeakRef(surface),
    second: new WeakRef(second),
    canvas: new WeakRef(surface.element.querySelector('canvas') as object),
  };
  // The completing context disposes both surfaces; the host removes their parents.
  ctx.dispose();
  host.remove();
  host2.remove();
  kept?.push(ctx);
  return refs;
}

/** What is still reachable of `refs` after collecting, a task apart, up to six times. */
async function collected(refs: Record<string, WeakRef<object>>) {
  const live = () =>
    Object.entries(refs)
      .filter(([, r]) => r.deref() !== undefined)
      .map(([k]) => k);
  for (let k = 0; k < 6 && live().length > 0; k++) {
    await new Promise((r) => setTimeout(r, 0));
    gc();
  }
  return live();
}

it('a disposed context with mounted surfaces, used, frees all of it once the host drops it', async () => {
  expect(await collected(await used())).toEqual([]);
});

it('the same with the host keeping the context is caught (positive control)', async () => {
  const kept: unknown[] = [];
  const live = await collected(await used(kept));
  expect(controlCaught(() => live.length, 1)).toBe(true);
  kept.length = 0;
});
