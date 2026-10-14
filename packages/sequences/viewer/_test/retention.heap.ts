/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Retention (IN1; seqspec:interaction#channels: nothing of the projector is module state): a disposed viewer frees its
// SeqSet and resolved view once the host drops it, after its projectors (the viewports' edge-rate one, a surface's)
// have projected, closed and focused through them. The positive control: the same, with the host keeping the set,
// collected and counted alike, is caught.

import type { Loci } from "@molstar/seq-core/model";
import { expect, it } from "vitest";
import { case14 } from "../../core/_test/format/cases.js";
import { SectionMarks } from "../src/marks.js";
import { controlCaught } from "./heap.js";
import { shown } from "./projection-fixtures.js";

const gc = (globalThis as { gc?: () => void }).gc as () => void;

/** A viewer used through its projectors, disposed; weak references to its set and view. `kept`: the host keeps the set. */
async function used(kept?: unknown[]) {
  const t = await shown(case14());
  const domain = t.ctx.refs.toLoci("domain") as Loci;
  expect(t.ctx.focus.fromLoci(domain, { transitionMs: 0 })).toBe(true);
  const marks = new SectionMarks();
  t.p.project(t.extents(domain, "containing"), 0, marks.hover);
  t.p.closeLoci(domain, "containing");
  const refs = { set: new WeakRef(t.set), view: new WeakRef(t.view) };
  t.ctx.dispose();
  kept?.push(t.set);
  return refs;
}

/** How many of the set and the view are still reachable after collecting, a task apart, up to four times. */
async function collected(refs: { set: WeakRef<object>; view: WeakRef<object> }) {
  for (let k = 0; k < 4 && refs.set.deref(); k++) {
    await new Promise((r) => setTimeout(r, 0));
    gc();
  }
  return (refs.set.deref() ? 1 : 0) + (refs.view.deref() ? 1 : 0);
}

it("a disposed viewer frees its SeqSet and resolved view", async () => {
  expect(await collected(await used())).toBe(0);
});

it("the same with the host keeping the set is caught (positive control)", async () => {
  const kept: unknown[] = [];
  const live = await collected(await used(kept));
  expect(controlCaught(() => live, 1)).toBe(true);
  kept.length = 0;
});
