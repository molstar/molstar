/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import type { Origin } from "../src/index.js";

import { event, live, state } from "../src/state.js";

describe("state", () => {
  it("replays the current value, emits only on change, shallow for plain records", () => {
    const s = state({ a: 1, b: "x" });
    const seen: unknown[] = [];
    s.subscribe((v) => seen.push(v));
    s.set({ a: 1, b: "x" });
    s.update({ a: 2 });
    expect(seen).toEqual([
      { a: 1, b: "x" },
      { a: 2, b: "x" },
    ]);
    expect(s.value).toEqual({ a: 2, b: "x" });
  });

  it("a set during delivery arrives after, in order, to every subscriber", () => {
    const s = state(0);
    const first: [number, Origin | undefined][] = [];
    const second: number[] = [];
    s.subscribe((v) => {
      first.push([v, s.origin]);
      if (v === 1) {
        s.set(2, "behavior");
        s.set(3, "api");
        expect(s.value).toBe(3);
      }
    });
    s.subscribe((v) => second.push(v));
    s.set(1, "pointer");
    expect(first).toEqual([
      [0, undefined],
      [1, "pointer"],
      [2, "behavior"],
      [3, "api"],
    ]);
    // The second subscriber sees 1 before the queued values, never 3 before 1.
    expect(second).toEqual([0, 1, 2, 3]);
  });

  it("a subscription made during a delivery replays the value being delivered, then the queued ones", () => {
    const s = state(0);
    const late: number[] = [];
    s.subscribe((v) => {
      if (v !== 1) return;
      s.set(2);
      s.subscribe((x) => late.push(x));
    });
    s.set(1);
    expect(late).toEqual([1, 2]);
  });

  it("the queue grows past its ring and keeps order", () => {
    const s = state(0);
    const seen: number[] = [];
    s.subscribe((v) => {
      seen.push(v);
      if (v === 1) for (let i = 2; i <= 20; i++) s.set(i);
    });
    s.set(1);
    expect(seen).toEqual(Array.from({ length: 21 }, (_, i) => i));
  });

  it("works with rxjs operators and completes", () => {
    const s = state(1);
    let done = false;
    const seen: number[] = [];
    s.subscribe({ next: (v) => seen.push(v), complete: () => (done = true) });
    s.complete();
    s.set(5);
    expect(seen).toEqual([1]);
    expect(done).toBe(true);
  });

  it("a live record bumps its version", () => {
    const l = live({ x: 0 });
    const versions: number[] = [];
    l.version.subscribe((v) => versions.push(v));
    l.current.x = 3;
    l.bump();
    expect(versions).toEqual([0, 1]);
  });
});

describe("events", () => {
  it("subscribers run by priority, higher first, then subscription order; nothing replays", () => {
    const e = event<{ n: number }>();
    const order: string[] = [];
    e.subscribe(() => order.push("a0"));
    e.subscribe(() => order.push("b100"), { priority: 100 });
    e.subscribe(() => order.push("c50"), { priority: 50 });
    const d = e.subscribe(() => order.push("d100"), { priority: 100 });
    e.subscribe(() => order.push("e0"));
    e.emit({ n: 1 });
    expect(order).toEqual(["b100", "d100", "c50", "a0", "e0"]);
    d.unsubscribe();
    order.length = 0;
    e.emit({ n: 2 });
    expect(order).toEqual(["b100", "c50", "a0", "e0"]);
  });

  it("an unsubscribe during emission skips that subscriber; a throw does not stop the others", () => {
    const errors: unknown[] = [];
    const e = event<number>((err) => errors.push(err));
    const seen: string[] = [];
    let b: { unsubscribe(): void } | null = null;
    e.subscribe(
      () => {
        seen.push("a");
        b?.unsubscribe();
        throw new Error("boom");
      },
      { priority: 1 },
    );
    b = e.subscribe(() => seen.push("b"));
    e.subscribe(() => seen.push("c"));
    e.emit(1);
    expect(seen).toEqual(["a", "c"]);
    expect(errors).toHaveLength(1);
  });
});
