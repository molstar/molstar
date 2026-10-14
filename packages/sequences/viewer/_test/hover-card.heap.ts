/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's hover handler (seqspec:interaction#the-hover-card) at pointer rate: the hover's version bumped with
// the pointer on an item whose card is open (what each move within the item does to the handler) allocates nothing, on
// RB0b's stepped guard. The positive control: the same bumps allocating an object each (`leaking`), measured alike,
// are caught.

import { expect, it } from "vitest";
import { case13 } from "../../core/_test/format/cases.ts";
import { createHoverCard } from "../src/hover-card.ts";
import { type ActionEntry, Hit, type ItemTarget, SeqViewerConfig } from "../src/index.ts";
import { state } from "../src/state.ts";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "./heap.ts";
import { opened } from "./layout-fixtures.ts";

const SECTION = "features.main";
const N = 20_000;
const Show: ActionEntry<ItemTarget> = { name: "test:show", title: "Show", when: () => true, run: () => {} };
const t = await opened(case13(), { size: { [SECTION]: [1000, 900] }, registry: [{ actions: [Show] }] });
t.ctx.config.set(SeqViewerConfig.UI.ItemActions, ["test:show"]);
// The card's own view of the hover: the context's hover record, a version only it listens to (the viewer's other
// hover subscribers, guarded elsewhere, stay out of this measurement).
const version = state(0);
const card = createHoverCard(
  Object.assign(Object.create(t.ctx), { hover: { live: { current: t.ctx.hover.live.current, version } } }),
);
const a = Hit.create();
t.ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, 150), t.yOf(SECTION, "features.domain"), a);
t.ctx.hover.setFromHit(a, "pointer");
version.set(1);
t.scheduler.advance(500);
const move = () => version.set(version.value + 1);

it("a hover change on the open card's item allocates nothing", () => {
  expect(card.current.value).not.toBe(null);
  const r = heapScaling(N, move);
  expect(isFlat(r), describeScaling(r, N)).toBe(true);
});

it("the same bumps allocating an object each are caught (positive control)", () => {
  expect(controlCaught(() => heapScaling(N, leaking(move)))).toBe(true);
});
