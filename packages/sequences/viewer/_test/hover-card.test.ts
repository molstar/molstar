/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover card's controller (seqspec:interaction#the-hover-card): inert without `UI.ItemActions`; after the delay on
// an item whose actions apply, open with the hover's label and those actions; kept across the gap into the card and
// while the pointer is in it; closed after the grace off both, by a zoom, and by running an action (which gets the item
// as its target).

import { describe, expect, it } from "vitest";
import { case13 } from "../../core/_test/format/cases.js";
import { trackMenu } from "../src/chrome.js";
import { createHoverCard, HOVER_CARD_GRACE_MS } from "../src/hover-card.js";
import { type ActionEntry, Hit, HoverKind, type ItemTarget, SeqViewerConfig } from "../src/index.js";
import { opened } from "./layout-fixtures.js";

const SECTION = "features.main";
const ran: ItemTarget[] = [];
const ShowIt: ActionEntry<ItemTarget> = {
  name: "test:show",
  title: "Show it",
  when: (_ctx, t) => t.annotation === "p53-features",
  run: (_ctx, t) => void ran.push(t),
};

async function viewer(actions: string[]) {
  const t = await opened(case13(), { size: { [SECTION]: [1000, 900] }, registry: [{ actions: [ShowIt] }] });
  t.ctx.config.set(SeqViewerConfig.UI.ItemActions, actions);
  const card = createHoverCard(t.ctx);
  const hit = Hit.create();
  /** The pointer on the Domain track at a residue (an item there), or off every item. */
  const hover = (index: number | null) => {
    const y = t.yOf(SECTION, "features.domain");
    t.ctx.layout.hitTest(SECTION, index === null ? 2 : t.xOf(SECTION, 0, index), index === null ? 2 : y, hit);
    t.ctx.hover.setFromHit(hit, "pointer");
    return hit.kind;
  };
  return { ...t, card, hover };
}

describe("the hover card", () => {
  it("opens after the delay on an item with actions, with the label and the actions that apply", async () => {
    const t = await viewer(["test:show", "test:missing"]);
    expect(t.hover(150)).toBe(HoverKind.item);
    t.scheduler.advance(399);
    expect(t.card.current.value).toBe(null);
    t.scheduler.advance(1);
    const c = t.card.current.value;
    expect(c?.actions).toEqual([{ name: "test:show", title: "Show it" }]);
    expect(c?.target).toMatchObject({ section: SECTION, track: "features.domain", annotation: "p53-features" });
    expect(c?.label?.title).toBeTruthy();
    t.card.dispose();
  });

  it("stays across the gap into the card and while the pointer is in it; closes after the grace off both", async () => {
    const t = await viewer(["test:show"]);
    t.hover(150);
    t.scheduler.advance(400);
    t.hover(null);
    t.scheduler.advance(HOVER_CARD_GRACE_MS - 50);
    t.card.enter();
    t.scheduler.advance(HOVER_CARD_GRACE_MS * 3);
    expect(t.card.current.value).not.toBe(null);
    t.card.leave();
    t.scheduler.advance(HOVER_CARD_GRACE_MS);
    expect(t.card.current.value).toBe(null);
    t.card.dispose();
  });

  it("closes on a zoom; running an action closes it and hands the item to the action", async () => {
    const t = await viewer(["test:show"]);
    t.hover(150);
    t.scheduler.advance(400);
    t.vp(SECTION).zoomAt(2, 500, "pointer");
    t.scheduler.advance(20);
    expect(t.card.current.value).toBe(null);
    t.hover(null);
    t.hover(150);
    t.scheduler.advance(400);
    const item = t.card.current.value?.target.item;
    t.card.run("test:show");
    expect(t.card.current.value).toBe(null);
    expect(ran.at(-1)).toMatchObject({ annotation: "p53-features", item });
    t.card.dispose();
  });

  it("stays across a reassembly of the same document (a skin's dark); a new document closes it (VB-07)", async () => {
    const t = await viewer(["test:show"]);
    t.hover(150);
    t.scheduler.advance(400);
    const before = t.ctx.graph.assembled.value;
    t.ctx.setSkin({ ...t.ctx.skin.value, dark: !t.ctx.skin.value.dark });
    await t.ctx.graph.update();
    expect(t.ctx.graph.assembled.value).not.toBe(before);
    expect(t.card.current.value).not.toBe(null);
    await t.ctx.setDocument(case13());
    expect(t.card.current.value).toBe(null);
    t.card.dispose();
  });

  it("is inert without UI.ItemActions, and on items no action applies to", async () => {
    const none = await viewer([]);
    none.hover(150);
    none.scheduler.advance(1000);
    expect(none.card.current.value).toBe(null);
    const other = await viewer(["test:show"]);
    expect(other.hover(null)).not.toBe(HoverKind.item);
    other.scheduler.advance(1000);
    expect(other.card.current.value).toBe(null);
  });

  it("moving from one item (its card open) to another opens the other's card after the delay", async () => {
    const t = await viewer(["test:show"]);
    t.hover(150);
    t.scheduler.advance(400);
    const first = t.card.current.value?.target.item;
    // Another item of the same track: a motif or region a residue range away (the first one found).
    let other = -1;
    for (const i of [310, 340, 20, 50, 360]) {
      t.hover(i);
      const now = t.ctx.hover.live.current.item.item;
      if (now !== first && now >= 0) {
        other = now;
        break;
      }
    }
    expect(other).toBeGreaterThanOrEqual(0);
    t.scheduler.advance(HOVER_CARD_GRACE_MS);
    t.scheduler.advance(400);
    expect(t.card.current.value?.target.item).toBe(other);
    t.card.dispose();
  });
});

describe("item actions from the keyboard (VB-08)", () => {
  it("the track's menu lists the cursor's item's actions first, and runs them on the item", async () => {
    const t = await viewer(["test:show"]);
    t.hover(150);
    const item = t.ctx.hover.live.current.item.item;
    // The keyboard's cursor on the same item: the hover is the keyboard's.
    const h = Hit.create();
    t.ctx.layout.hitTest(SECTION, t.xOf(SECTION, 0, 150), t.yOf(SECTION, "features.domain"), h);
    t.ctx.hover.setFromHit(h, "keyboard");
    const menu = trackMenu(t.ctx, SECTION, "features.domain");
    expect(menu[0]).toMatchObject({ key: `item:test:show:p53-features:${item}:-1`, label: "Show it" });
    menu[0].run();
    expect(ran.at(-1)).toMatchObject({ annotation: "p53-features", item, track: "features.domain" });
    // The pointer's hover is the card's: the menu lists no item actions then; nor on another track.
    t.hover(150);
    expect(trackMenu(t.ctx, SECTION, "features.domain").some((m) => m.key.startsWith("item:"))).toBe(false);
    t.ctx.hover.setFromHit(h, "keyboard");
    const view = t.ctx.graph.assembled.value?.view;
    const other = Array.from(view?.resolved.tracks.ref ?? [], (r) => view?.set.strings.get(r) ?? "").find(
      (r) => r.startsWith("features.") && r !== "features.domain",
    );
    expect(other).toBeTruthy();
    expect(trackMenu(t.ctx, SECTION, other as string).some((m) => m.key.startsWith("item:"))).toBe(false);
    t.card.dispose();
  });
});
