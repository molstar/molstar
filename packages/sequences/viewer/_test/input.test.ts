/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The input dispatcher and base events (seqspec:interaction#input-events): claims by priority, the slop, clicks and
// double clicks on the manual clock, a release outside the surface, keys while typing; then the findings it closes: one
// pointer, cancel, a double click by button and section, and a re-hit that follows a rebuild.

import { describe, expect, it } from "vitest";
import { case2, case3 } from "../../core/_test/format/cases.ts";
import { Button, type DragEvt, DragPhase, Hit, HitRegion, Inputs, PointerType } from "../src/index.ts";
import type { InputDispatcher } from "../src/input.ts";
import { driven } from "./input-fixtures.ts";

const S = "family.main";
const open = () => driven(case3(), { size: { [S]: [1000, 150] } });
/** On the members' rows (data), in the fourth column. */
const rowY = (t: Awaited<ReturnType<typeof open>>, d: number) => t.layout(S).scrolledTop + d * 18 + 9;

describe("base events", () => {
  it("pointer writes ctx.pointer first; hover goes out only when the element changes", async () => {
    const t = await open();
    const seen: number[] = [];
    let hovers = 0;
    t.ctx.interaction.pointer.subscribe((e) => seen.push(t.ctx.pointer.current.x, e.hit.region));
    t.ctx.interaction.hover.subscribe(() => hovers++);
    const x = t.xOf(S, 0, 4);
    t.pointer("move", S, x, rowY(t, 0));
    t.pointer("move", S, x + 1, rowY(t, 0)); // the same cell
    t.pointer("move", S, x + 1, rowY(t, 1));
    expect(seen).toEqual([x, HitRegion.data, x + 1, HitRegion.data, x + 1, HitRegion.data]);
    expect(hovers).toBe(2);
    t.ctx.input.leave(0);
    expect([t.ctx.pointer.current.inside, t.ctx.pointer.current.section]).toEqual([false, -1]);
  });

  it("claims: the higher priority first; later subscribers see the owner; ui.drag names it", async () => {
    const t = await open();
    // Without the defaults, to see the claims alone.
    for (const b of ["select", "dragPan"]) t.ctx.behaviors.get(b)?.remove();
    const order: string[] = [];
    const phases: number[] = [];
    t.ctx.interaction.drag.subscribe(
      (e) => {
        order.push(`low:${e.owner}`);
        e.claim("low");
      },
      { priority: 1 },
    );
    t.ctx.interaction.drag.subscribe(
      (e) => {
        if (e.phase === DragPhase.start) order.push(`high:${e.owner}`);
        e.claim("high");
        phases.push(e.phase);
      },
      { priority: 5 },
    );
    const owners: (string | null)[] = [];
    t.ctx.ui.drag.subscribe((d) => owners.push(d?.owner ?? null));
    t.drag(S, 400, rowY(t, 0), 440, rowY(t, 0));
    expect(order.slice(0, 2)).toEqual(["high:", "low:high"]);
    expect(phases).toEqual([DragPhase.start, 1, 1, 1, DragPhase.end]);
    expect(owners).toEqual([null, "high", null]);
  });

  it("an unclaimed drag does nothing: no moves or end, no click", async () => {
    const t = await open();
    for (const b of ["select", "dragPan"]) t.ctx.behaviors.get(b)?.remove();
    const phases: number[] = [];
    let clicks = 0;
    t.ctx.interaction.drag.subscribe((e) => phases.push(e.phase));
    t.ctx.interaction.click.subscribe(() => clicks++);
    t.drag(S, 400, rowY(t, 0), 460, rowY(t, 0));
    expect([phases, clicks, t.ctx.ui.drag.value]).toEqual([[DragPhase.start], 0, null]);
  });

  it("a press claimed on down skips the slop", async () => {
    const t = await open();
    const starts: number[] = [];
    t.ctx.interaction.pointer.subscribe((e) => {
      if (e.input.type === PointerType.down) e.claim("lab:grab");
    });
    t.ctx.interaction.drag.subscribe((e) => starts.push(e.phase, e.totalX));
    t.pointer("down", S, 400, rowY(t, 0));
    t.pointer("move", S, 401, rowY(t, 0)); // within the 3 px slop
    expect(starts).toEqual([DragPhase.start, 1]);
    expect(t.ctx.ui.drag.value?.owner).toBe("lab:grab");
  });

  it("clicks within the slop; a double and a triple click within DoubleClickMs on the manual clock", async () => {
    const t = await open();
    const counts: number[] = [];
    t.ctx.interaction.click.subscribe((e) => counts.push(e.count));
    const [x, y] = [t.xOf(S, 0, 4), rowY(t, 0)];
    t.click(S, x, y);
    t.advance(200);
    t.click(S, x + 1, y); // a double click
    t.advance(200);
    t.click(S, x, y); // a triple click
    t.advance(200);
    t.click(S, x, y); // starts the count again
    t.advance(400);
    t.click(S, x, y); // too late
    t.pointer("down", S, x, y);
    t.pointer("move", S, x + 2, y);
    t.pointer("up", S, x + 2, y); // still within the slop
    expect(counts).toEqual([1, 2, 3, 1, 1, 2]);
  });

  it("a claimed press that wobbles within the slop is a drag that ends where it began, and a click", async () => {
    const t = await open();
    const seen: string[] = [];
    t.ctx.interaction.pointer.subscribe((e) => {
      if (e.input.type === PointerType.down) e.claim("lab:grab");
    });
    t.ctx.interaction.drag.subscribe((e) => seen.push(`drag ${e.phase} ${e.totalX}`));
    t.ctx.interaction.click.subscribe((e) => seen.push(`click ${e.count}`));
    const y = rowY(t, 0);
    const wobble = (x: number, dx: number) => {
      t.pointer("down", S, x, y);
      t.pointer("move", S, x + dx, y);
      t.pointer("up", S, x + dx, y);
    };
    wobble(400, 1);
    t.advance(50);
    wobble(400, -1);
    expect(seen).toEqual([
      `drag ${DragPhase.start} 1`,
      `drag ${DragPhase.end} 1`,
      "click 1",
      `drag ${DragPhase.start} -1`,
      `drag ${DragPhase.end} -1`,
      "click 2",
    ]);
    seen.length = 0;
    t.advance(1000);
    // Out past the slop and back: a drag only; held past the click time: a drag only.
    t.pointer("down", S, 400, y);
    t.pointer("move", S, 410, y);
    t.pointer("up", S, 400, y);
    t.pointer("down", S, 400, y);
    t.pointer("move", S, 401, y);
    t.advance(400);
    t.pointer("up", S, 401, y);
    expect(seen.filter((s) => s.startsWith("click"))).toEqual([]);
  });

  it("on the navigator: a double click with 1 px wobbles still toggles; a drag pans and is no click", async () => {
    const F = "fv.main";
    const t = await driven(case2(), { size: { [F]: [1000, 200] } });
    t.ctx.ui.setStyle(F, "tracks");
    await t.ctx.graph.update();
    const g = t.frame(F, 0);
    const s = t.vp(F).live.current;
    const counts: number[] = [];
    t.ctx.interaction.click.subscribe((e) => counts.push(e.count));
    const x = g[0] + 0.4 * g[1];
    const press = (dx: number) => {
      t.pointer("down", F, x, 10);
      t.pointer("move", F, x + dx, 10);
      t.pointer("up", F, x + dx, 10);
    };
    press(1);
    t.advance(50);
    press(1);
    t.advance(1000);
    const target = Math.max(t.vp(F).props.minLetterCell, (s.width - (t.vp(F).axis?.insetX ?? 0)) / 60);
    expect([counts, s.cell]).toEqual([[1, 2], target]);

    const sx = s.scrollX;
    t.drag(F, x, 10, x + 40, 10);
    expect(s.scrollX).toBeGreaterThan(sx);
    expect(counts).toEqual([1, 2]);

    press(-1);
    t.advance(50);
    press(1);
    t.advance(1000);
    expect([counts, s.cell, s.scrollX]).toEqual([[1, 2, 1, 2], s.minCell, 0]);
  });

  it("a release outside the surface ends a drag", async () => {
    const t = await open();
    t.vp(S).zoomTo(0, 0, 10);
    const ends: DragEvt["hit"]["region"][] = [];
    t.ctx.interaction.drag.subscribe((e) => e.phase === DragPhase.end && ends.push(e.hit.region), { priority: -1 });
    t.pointer("down", S, 400, rowY(t, 0));
    t.pointer("move", S, 300, rowY(t, 0));
    t.pointer("up", S, -40, 500, { onDocument: true });
    expect([ends, t.ctx.ui.drag.value]).toEqual([[HitRegion.none], null]);
  });

  it("keys are ignored while typing", async () => {
    const t = await open();
    let keys = 0;
    t.ctx.interaction.key.subscribe(() => keys++);
    t.key(S, "Escape", { typing: true });
    t.key(S, "Escape");
    expect(keys).toBe(1);
  });
});

describe("the findings the dispatcher closes (IN2's hot-path design)", () => {
  it("a drag ends as a cancel by cancel() and its surface's dispose; a release never seen ends it at the last point", async () => {
    const t = await open();
    const seen: string[] = [];
    t.ctx.interaction.drag.subscribe((e) => seen.push(`${e.phase}:${e.owner}`));
    const y = rowY(t, 0);
    const start = () => {
      t.pointer("down", S, 500, y);
      t.pointer("move", S, 480, y);
      expect(t.ctx.ui.drag.value?.owner).toBe("dragPan");
    };
    const C = `${DragPhase.cancel}:dragPan`;
    start();
    t.ctx.input.cancel();
    expect([t.ctx.ui.drag.value, seen.at(-1)]).toEqual([null, C]);
    t.pointer("move", S, 470, y); // a hover now: nothing drags
    expect(t.ctx.ui.drag.value).toBeNull();

    // A move with no button down: the release was lost; the drag ends (not a cancel) where it last was, and the move
    // is a hover after it.
    const E = `${DragPhase.end}:dragPan`;
    start();
    const sx = t.vp(S).live.current.scrollX;
    const lost = t.ctx.pointer.version.value;
    const ends: number[] = [];
    const sub = t.ctx.interaction.drag.subscribe((e) => e.phase === DragPhase.end && ends.push(e.input.x));
    t.ctx.input.pointer(Object.assign(sinkRecord(), { type: PointerType.move, section: 0, x: 460, y, buttons: 0 }));
    sub.unsubscribe();
    expect([t.ctx.ui.drag.value, seen.at(-1), t.vp(S).live.current.scrollX === sx, ends]).toEqual([
      null,
      E,
      true,
      [480],
    ]);
    expect(t.ctx.pointer.current.x).toBe(460);
    expect(t.ctx.pointer.version.value).toBeGreaterThan(lost);

    start();
    t.ctx.input.cancel(5); // another surface's dispose: nothing
    expect(t.ctx.ui.drag.value?.owner).toBe("dragPan");
    t.ctx.input.cancel(0);
    expect([t.ctx.ui.drag.value, seen.at(-1)]).toEqual([null, C]);
    // A press never released (a down again) ends first, released at its last point; then the new press.
    start();
    t.pointer("down", S, 500, y);
    expect([t.ctx.ui.drag.value, seen.at(-1)]).toEqual([null, E]);
    t.pointer("up", S, 500, y);
  });

  it("another pointer's moves and up are ignored while a press lasts; its down ends a press whose release was lost (VS-05)", async () => {
    const t = await open();
    const y = rowY(t, 0);
    const vp = t.vp(S);
    vp.zoomTo(0, 5, 15);
    const s = vp.live.current;
    t.pointer("down", S, 500, y, { pointerId: 1 });
    t.pointer("move", S, 480, y, { pointerId: 1 });
    const at = s.scrollX;
    t.pointer("move", S, 90, y, { pointerId: 2 });
    t.pointer("up", S, 90, y, { pointerId: 2 });
    expect([t.ctx.ui.drag.value?.owner, s.scrollX]).toEqual(["dragPan", at]);
    t.pointer("move", S, 470, y, { pointerId: 1 });
    expect(s.scrollX - at).toBe(10);
    // A touch whose release the page never delivered: the next touch (a new id) ends it, then is a press of its own.
    const clicks: number[] = [];
    t.ctx.interaction.click.subscribe((e) => clicks.push(e.count));
    t.pointer("down", S, 300, y, { pointerId: 7 });
    expect(t.ctx.ui.drag.value).toBeNull();
    t.pointer("up", S, 300, y, { pointerId: 7 });
    expect(clicks).toEqual([1]);
  });

  it("right then left counts 1, 1; the same button on another section is no double click", async () => {
    const t = await open();
    const counts: number[] = [];
    t.ctx.interaction.click.subscribe((e) => counts.push(e.count));
    const [x, y] = [t.xOf(S, 0, 4), rowY(t, 0)];
    t.click(S, x, y, { button: Button.secondary });
    t.advance(100);
    t.click(S, x, y, { button: Button.primary });
    t.advance(100);
    t.click(S, x, y, { button: Button.primary });
    expect(counts).toEqual([1, 1, 2]);
  });

  it("a rebuild inside the dispatcher's own hit test is not hit again; one during a click is, once, after it", async () => {
    const t = await open();
    const [x, y] = [t.xOf(S, 0, 4), rowY(t, 0)];
    let hits = 0;
    const layout = t.ctx.layout;
    const hitTestAt = layout.hitTestAt.bind(layout);
    layout.hitTestAt = (section, out) => {
      hits++;
      return hitTestAt(section, out);
    };
    t.pointer("move", S, x, y);
    // Dirty: the next read rebuilds, inside the next move's own hit test, which reads the new layout (one hit, and
    // the hover bumped for the new element only).
    const builds = t.layouts.stats.builds;
    const bumps = t.ctx.hover.live.version.value;
    hits = 0;
    t.ctx.ui.toggleClade("sh2", 1);
    t.pointer("move", S, x, y + 18);
    expect([t.layouts.stats.builds, hits, t.ctx.hover.live.version.value - bumps]).toEqual([builds + 1, 1, 1]);
    // A click's handler that rebuilds it: the pointer is hit again after the dispatch (down 2, up 1, the re-hit 1), and
    // the same element refreshes the hover (a clade's "click expand").
    const sub = t.ctx.interaction.click.subscribe(() => {
      t.ctx.ui.toggleClade("sh2", 1);
      t.layouts.axis(S);
    });
    const before = t.ctx.hover.live.version.value;
    hits = 0;
    t.click(S, x, y + 18);
    sub.unsubscribe();
    expect([hits, t.ctx.hover.live.version.value > before]).toEqual([4, true]);
    // Two sections rebuilt during one dispatch, the pointer's first: the pointer's is the one hit again.
    const input = t.ctx.input as InputDispatcher;
    const both = t.ctx.interaction.click.subscribe(() => {
      input.rehit(0);
      input.rehit(5);
    });
    const was = t.ctx.hover.live.version.value;
    t.click(S, x, y + 18);
    both.unsubscribe();
    expect(t.ctx.hover.live.version.value > was).toBe(true);
  });

  it("a hit test on a section with no layout resets the hit, not the last element", async () => {
    const t = await open();
    const h = Hit.create();
    expect(t.ctx.layout.hitTestAt(0, Object.assign(h, { x: t.xOf(S, 0, 4), y: rowY(t, 0) }))).toBe(true);
    expect(h.region).toBe(HitRegion.data);
    expect(t.ctx.layout.hitTestAt(7, h)).toBe(false);
    expect([h.region, h.block, h.shown, h.section, h.system.kind]).toEqual([HitRegion.none, -1, -1, -1, -1]);
  });
});

const sinkRecord = () => Inputs.pointer();
