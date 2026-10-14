/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dispatcher's half of "the hover follows the content under a resting pointer" (seqspec:interaction#input-events):
// a moved view is hit again once, at the next frame or the engine's draw before it (`flush`), with no frame while the
// pointer is elsewhere, pressed or the keyboard holds the hover; many moves in a frame are one re-hit; a hover
// subscriber that moves the
// view is no loop. The browser half (both UIs, real events, the tooltip) is canvas/resting-hover.test.tsx.

import { describe, expect, it } from "vitest";
import { Hit, HoverKind } from "../src/index.ts";
import type { InputDispatcher } from "../src/input.ts";
import { FAMILY_SECTION as S, syntheticFamily } from "./families.ts";
import { driven } from "./input-fixtures.ts";

const open = async () => {
  const t = await driven(syntheticFamily(120, 8), { size: { [S]: [420, 260] } });
  // Hit tests counted: the dispatcher's own and the re-hits.
  const counter = { hits: 0 };
  const layout = t.ctx.layout;
  const hitTestAt = layout.hitTestAt.bind(layout);
  layout.hitTestAt = (section, out) => {
    counter.hits++;
    return hitTestAt(section, out);
  };
  const rh = t.layout(S).rows.height;
  const y = t.layout(S).scrolledTop + 1.5 * rh;
  const x = t.xOf(S, 0, 4) + 0.25;
  return { t, counter, rh, x, y, input: t.ctx.input as InputDispatcher };
};
type Open = Awaited<ReturnType<typeof open>>;
/** The hovered row and column, and what is under the pointer now. */
function cells({ t }: Open) {
  const h = t.ctx.hover.live.current;
  const p = t.ctx.pointer.current;
  const g = Hit.create();
  t.ctx.layout.hitTest(S, p.x, p.y, g);
  return {
    hover: h.kind === HoverKind.aligned ? [h.aligned.row, h.aligned.column] : [h.kind],
    under: g.kind === HoverKind.aligned ? [g.aligned.row, g.aligned.column] : [g.kind],
  };
}

describe("a moved view under a resting pointer", () => {
  it("is hit again at the next frame, once for many moves in it; the hover is what is under the pointer", async () => {
    const o = await open();
    const { t, counter, rh, x, y } = o;
    t.pointer("move", S, x, y);
    t.advance(16);
    const at = cells(o).hover;
    expect(at.length).toBe(2);
    counter.hits = 0;
    const v = t.ctx.hover.live.version.value;
    for (let k = 0; k < 4; k++) t.wheel(S, x, y, 0, rh / 2);
    // Each wheel hits once for itself; the hover waits for the frame.
    expect([counter.hits, cells(o).hover]).toEqual([4, at]);
    t.advance(16);
    const c = cells(o);
    expect(c.hover).toEqual(c.under);
    expect(c.hover).toEqual([(at[0] as number) + 2, at[1]]);
    expect([counter.hits, t.ctx.hover.live.version.value - v]).toEqual([5, 1]);
    // Nothing moved: no frame, no hit.
    t.advance(16);
    expect(counter.hits).toBe(5);
  });

  it("a view that moved within the same cell: no hover bump (no refresh, unlike a rebuild)", async () => {
    const o = await open();
    const { t, x, y } = o;
    t.pointer("move", S, x, y);
    const v = t.ctx.hover.live.version.value;
    t.vp(S).panBy(0, 1, "api");
    t.advance(16);
    expect(t.ctx.hover.live.version.value).toBe(v);
  });

  it("`flush` runs the waiting re-hit now (the engine draws the hover with the view); its frame then finds nothing", async () => {
    const o = await open();
    const { t, counter, rh, x, y, input } = o;
    t.pointer("move", S, x, y);
    const at = cells(o).hover;
    t.vp(S).panBy(0, 3 * rh, "api");
    counter.hits = 0;
    input.flush();
    expect([counter.hits, cells(o).hover]).toEqual([1, [(at[0] as number) + 3, at[1]]]);
    // Another move before the frame: no second request, the same frame.
    const frames = t.scheduler.pendingFrames;
    t.vp(S).panBy(0, rh, "api");
    input.flush();
    expect([counter.hits, t.scheduler.pendingFrames]).toEqual([2, frames]);
    input.flush();
    t.advance(16);
    expect(counter.hits).toBe(2);
  });

  it("nothing while the pointer is elsewhere, after a leave, during a press, or under the keyboard's hover", async () => {
    const o = await open();
    const { t, counter, rh, x, y, input } = o;
    const moved = () => {
      counter.hits = 0;
      t.vp(S).panBy(0, rh, "api");
      t.advance(16);
      return counter.hits;
    };
    expect(moved()).toBe(0); // never over the surface
    t.pointer("move", S, x, y);
    input.leave(t.ctx.layout.sectionIndex(S));
    expect([moved(), t.ctx.hover.live.current.kind]).toEqual([0, HoverKind.none]);
    t.pointer("move", S, x, y);
    t.pointer("down", S, x, y);
    expect(moved()).toBe(0);
    t.pointer("up", S, x, y);
    // The keyboard's cursor holds the hover: it stays.
    const h = Hit.create();
    t.ctx.layout.hitTest(S, x, y, h);
    t.ctx.hover.setFromHit(h, "keyboard");
    const v = t.ctx.hover.live.version.value;
    expect([moved(), t.ctx.hover.live.version.value, t.ctx.hover.live.current.owner]).toEqual([0, v, "keyboard"]);
  });

  it("no loop: a hover subscriber that moves the view is not hit again for it", async () => {
    const o = await open();
    const { t, counter, rh, x, y } = o;
    t.pointer("move", S, x, y);
    const sub = t.ctx.hover.live.version.subscribe(() => t.vp(S).panBy(0, rh, "api"));
    counter.hits = 0;
    t.vp(S).panBy(0, rh, "api");
    for (let k = 0; k < 5; k++) t.advance(16);
    sub.unsubscribe();
    expect(counter.hits).toBe(1);
  });

  it("an animated move: the hover follows each frame of it, as the engine draws it (`flush` before the draw)", async () => {
    const o = await open();
    const { t, x, y, input } = o;
    t.vp(S).zoomTo(0, 0, 8);
    t.advance(1000);
    t.pointer("move", S, x, y);
    const seen = new Set<number>();
    t.vp(S).zoomTo(0, 20, 28, { animate: true, durationMs: 200 });
    for (let k = 0; k < 14; k++) {
      t.advance(16);
      // The animation's frame moved the view after the re-hit's frame in this one: the draw's flush hits it again.
      input.flush();
      const c = cells(o);
      expect(c.hover).toEqual(c.under);
      seen.add(c.hover[1] as number);
    }
    expect(seen.size).toBeGreaterThan(3);
  });
});
