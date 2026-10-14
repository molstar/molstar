/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Menus and keyboard focus in the React UI (UX3; UI2's cases, as _test/dom/focus.test.ts has them for the vanilla UI;
// seqspec:interaction#slots-and-chrome, seqspec:interaction#accessibility): a menu opens on its trigger or on the
// keyboard cursor's Shift+F10, focuses its first item, moves with the arrows, closes on Escape, on a choice and on an
// outside press and gives focus back to its trigger; a host's Menu gets `open` and answers through `onOpenChange`; a
// node React removes lands focus on the next of its group, the previous, the section's surface; the issue list opens
// with the same contract.

import { act, cloneElement, type ReactElement, useState } from "react";
import { describe, expect, it, vi } from "vitest";
import "../../src/skin.css";
import { withDownloadExport } from "../../src/dom/index.ts";
import type { SeqViewerContext } from "../../src/index.ts";
import { createSeqViewer, manualScheduler } from "../../src/index.ts";
import { Menu, type MenuProps, SeqViewerLayout, SeqViewerProvider, type SlotName } from "../../src/react/index.ts";
import { fixtureMeasurer } from "../../src/testing/index.ts";
import { patch } from "../graph-fixtures.ts";
import { at, harness, headlessSpec, keys, pointer, settle, stubSurface, uiDoc, xOf } from "./fixtures.tsx";

type Trigger = ReactElement<{ "aria-expanded"?: boolean; onClick?: () => void }>;

const { opened, cleanups, parent } = harness();
const S = "chains.main";
const active = () => document.activeElement as HTMLElement | null;
const press = async (el: Element, key: string, init: KeyboardEventInit = {}) => {
  const e = new KeyboardEvent("keydown", { key, bubbles: true, cancelable: true, composed: true, ...init });
  await act(async () => void el.dispatchEvent(e));
  return e;
};
const surfaceOf = (c: HTMLElement) => c.querySelector(".seqv-surface [tabindex]") as HTMLElement;
const trigger = (c: HTMLElement) => c.querySelector(".seqv-track-actions button") as HTMLButtonElement;
const items = (c: HTMLElement) => [...c.querySelectorAll<HTMLButtonElement>("[role^=menuitem]")];
const click = (el: Element) => act(async () => (el as HTMLElement).click());

/** The track's box at the Domains header, the pointer on it. */
async function atDomains(ui: Parameters<typeof opened>[0] = {}, strict = false) {
  const t = await opened(ui, strict);
  const h = at(t.ctx, S, "t.dom");
  await act(async () => pointer(t.ctx)("move", S, h.header, h.y));
  return t;
}
/** The keyboard cursor on the Domains track. */
async function cursorOnDomains(ctx: SeqViewerContext) {
  const key = keys(ctx);
  await act(async () => ctx.input.focus(ctx.layout.sectionIndex(S), true));
  for (let i = 0; i < 8 && ctx.ui.keyCursor.value?.track !== "t.dom"; i++) await act(async () => key(S, "ArrowDown"));
  expect(ctx.ui.keyCursor.value?.track).toBe("t.dom");
  return key;
}

describe("the menu", () => {
  it("opens on its trigger with the first item focused, the arrows move between items, Home and End go to the ends", async () => {
    const { container } = await atDomains();
    const t = trigger(container);
    expect([t.getAttribute("aria-haspopup"), t.getAttribute("aria-expanded")]).toEqual(["menu", "false"]);
    await click(t);
    const all = items(container);
    expect(all.length).toBe(4);
    expect([t.getAttribute("aria-expanded"), t.getAttribute("aria-controls")]).toEqual([
      "true",
      container.querySelector("[role=menu]")?.id,
    ]);
    expect(active()).toBe(all[0]);
    await press(all[0], "ArrowDown");
    expect(active()).toBe(all[1]);
    await press(all[1], "ArrowDown");
    await press(all[2], "ArrowDown");
    expect(active()).toBe(all[3]);
    await press(all[3], "ArrowDown");
    expect(active()).toBe(all[0]);
    await press(all[0], "ArrowUp");
    expect(active()).toBe(all[3]);
    await press(all[3], "Home");
    expect(active()).toBe(all[0]);
    await press(all[0], "End");
    expect(active()).toBe(all[3]);
    // Closed, ArrowDown on the focused trigger opens it at the first item.
    await press(active() as HTMLElement, "Escape");
    expect(active()).toBe(t);
    await press(t, "ArrowDown");
    expect([t.getAttribute("aria-expanded"), active()]).toEqual(["true", items(container)[0]]);
  });

  it("Escape closes it and focus goes back to the trigger; so does choosing an item, which still runs", async () => {
    const { container, ctx } = await atDomains();
    const t = trigger(container);
    await click(t);
    const pressed = await press(active() as HTMLElement, "Escape");
    expect(pressed.defaultPrevented).toBe(true);
    expect([container.querySelector("[role=menu]"), active(), t.getAttribute("aria-expanded")]).toEqual([
      null,
      t,
      "false",
    ]);
    // A choice that leaves the trigger in place: "Focus" moves the viewport.
    await click(t);
    const focus = items(container)[1];
    expect(focus.textContent).toBe("Focus");
    focus.focus();
    await click(focus);
    await settle();
    expect([container.querySelector("[role=menu]"), active()]).toEqual([null, t]);
    expect(ctx.focus.current.value.size).toBeGreaterThan(0);
  });

  it("an open menu reads its items again on the streams' version: a failed region adds Retry, focus stays on its item", async () => {
    const { container, ctx } = await atDomains();
    await click(trigger(container));
    const labels = () => items(container).map((i) => i.textContent);
    expect(labels()).not.toContain("Retry");
    const first = active()?.dataset.seqvKey;
    vi.spyOn(ctx.streams, "hasFailed").mockReturnValue(true);
    await act(async () =>
      (ctx.streams.version as unknown as { set(v: number): void }).set(ctx.streams.version.value + 1),
    );
    expect(labels().at(-1)).toBe("Retry");
    expect(active()?.dataset.seqvKey).toBe(first);
  });

  it("Tab closes it without trapping focus", async () => {
    const { container } = await atDomains();
    const t = trigger(container);
    await click(t);
    await press(active() as HTMLElement, "Tab");
    expect([container.querySelector("[role=menu]"), active()]).toEqual([null, t]);
  });

  it("focus moving out of it closes it and stays where it went", async () => {
    const { container } = await atDomains();
    await click(trigger(container));
    expect(container.querySelector("[role=menu]")).not.toBeNull();
    const surface = surfaceOf(container);
    await act(async () => surface.focus());
    expect([container.querySelector("[role=menu]"), active()]).toEqual([null, surface]);
  });

  it("a press outside closes it (a capture-phase pointerdown, which Safari needs), a press inside does not", async () => {
    const { container } = await atDomains();
    const t = trigger(container);
    await click(t);
    await act(async () =>
      items(container)[0].dispatchEvent(new PointerEvent("pointerdown", { bubbles: true, composed: true })),
    );
    expect(container.querySelector("[role=menu]")).not.toBeNull();
    expect(active()).toBe(items(container)[0]);
    await act(async () => document.body.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true })));
    expect([container.querySelector("[role=menu]"), t.getAttribute("aria-expanded")]).toEqual([null, "false"]);
    // Focus is the press's to move, not taken back to the trigger.
    expect(active()).not.toBe(t);
    await act(async () => document.body.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true })));
    expect(container.querySelector("[role=menu]")).toBeNull();
  });

  it("Shift+F10 and ContextMenu on the keyboard cursor's track open its menu, the first item focused", async () => {
    const { container, ctx } = await opened();
    const key = await cursorOnDomains(ctx);
    expect(container.querySelector("[role=menu]")).toBeNull();
    // F10 alone is not the key.
    let handled: boolean[] = [];
    await act(async () => void handled.push(key(S, "F10")));
    expect([handled, container.querySelector("[role=menu]")]).toEqual([[false], null]);
    handled = [];
    await act(async () => void handled.push(key(S, "F10", true)));
    expect(handled).toEqual([true]);
    expect(items(container).map((i) => i.textContent)).toEqual([
      "Show only what this track covers",
      "Focus",
      "Colors ▸",
      "Hide track",
    ]);
    expect(active()).toBe(items(container)[0]);
    await press(active() as HTMLElement, "Escape");
    expect([container.querySelector("[role=menu]"), active()?.getAttribute("aria-label")]).toEqual([
      null,
      "Actions for Domains",
    ]);
    await act(async () => void key(S, "ContextMenu"));
    expect(container.querySelector("[role=menu]")).not.toBeNull();
  });

  it("a host's Menu receives `open` on Shift+F10, and its `onOpenChange` closes it", async () => {
    const seen: boolean[] = [];
    let last: MenuProps | null = null;
    const HostMenu = (p: MenuProps) => {
      [last] = [p];
      seen.push(p.open);
      return (
        <span className="host-menu">
          {cloneElement(p.trigger as Trigger, {
            "aria-expanded": p.open,
            onClick: () => p.onOpenChange(!p.open),
          })}
        </span>
      );
    };
    const { container, ctx } = await opened({ components: { Menu: HostMenu } });
    const key = await cursorOnDomains(ctx);
    expect(container.querySelector(".host-menu")).not.toBeNull();
    expect(seen.at(-1)).toBe(false);
    await act(async () => void key(S, "F10", true));
    expect(seen.at(-1)).toBe(true);
    await act(async () => (last as MenuProps | null)?.onOpenChange(false));
    expect(seen.at(-1)).toBe(false);
  });

  it("an open menu keeps its box while the pointer is elsewhere and the cursor is gone, and focus stays in it", async () => {
    const { container, ctx } = await atDomains();
    await click(trigger(container));
    const first = items(container)[0];
    expect(active()).toBe(first);
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex(S)));
    await settle();
    const res = at(ctx, S, "t.res");
    await act(async () => pointer(ctx)("move", S, xOf(ctx, S, 12), res.y));
    await settle();
    expect([container.querySelector("[role=menu]"), active()]).toEqual([expect.anything(), first]);
    // Closed and focus moved off the box: the box goes, and focus lands on the surface.
    await press(first, "Escape");
    await act(async () => surfaceOf(container).focus());
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex(S)));
    await settle();
    expect(container.querySelector(".seqv-track-actions")).toBeNull();
    expect(active()).toBe(surfaceOf(container));
  });
});

describe("a default Menu answers its owner's latest callback", () => {
  it("an inline `onOpenChange` that changes each render is the one called", async () => {
    const calls: number[] = [];
    const items = [{ key: "a", label: "A", run: () => {} }];
    const Owner = ({ n }: { n: number }) => {
      const [open, setOpen] = useState(false);
      return (
        <Menu
          trigger={<button type="button">open</button>}
          items={items}
          open={open}
          onOpenChange={(o) => {
            calls.push(n);
            setOpen(o);
          }}
        />
      );
    };
    const p = parent((i) => <Owner n={i} />);
    await p.render();
    await p.render();
    await p.render();
    await act(async () => (p.container.querySelector("button") as HTMLButtonElement).click());
    await press(active() as HTMLElement, "Escape");
    expect(calls).toEqual([2, 2]);
  });
});

describe("the track's box", () => {
  it("stays while its closed trigger holds focus and the pointer is elsewhere, and goes when focus leaves it", async () => {
    const { container, ctx } = await atDomains();
    const t = trigger(container);
    await click(t);
    await press(active() as HTMLElement, "Escape");
    expect(active()).toBe(t);
    // Nothing hovered (the pointer on the data area, no keyboard cursor): focus inside keeps it.
    await act(async () => pointer(ctx)("move", S, xOf(ctx, S, 12), at(ctx, S, "t.res").y));
    await settle();
    expect([trigger(container), active()]).toEqual([t, t]);
    await act(async () => surfaceOf(container).focus());
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex(S)));
    await settle();
    expect(container.querySelector(".seqv-track-actions")).toBeNull();
  });

  it("stays while a host's Menu is open (its trigger `aria-expanded`), with focus outside the box", async () => {
    let last: MenuProps | null = null;
    const HostMenu = (p: MenuProps) => {
      last = p;
      return cloneElement(p.trigger as Trigger, { "aria-expanded": p.open });
    };
    const { container, ctx } = await atDomains({ components: { Menu: HostMenu } });
    await act(async () => (last as MenuProps | null)?.onOpenChange(true));
    expect(trigger(container).getAttribute("aria-expanded")).toBe("true");
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex(S)));
    await settle();
    expect(container.querySelector(".seqv-track-actions")).not.toBeNull();
    await act(async () => (last as MenuProps | null)?.onOpenChange(false));
    await act(async () => ctx.input.leave(ctx.layout.sectionIndex(S)));
    await act(async () => pointer(ctx)("move", S, xOf(ctx, S, 12), at(ctx, S, "t.res").y));
    await settle();
    expect(container.querySelector(".seqv-track-actions")).toBeNull();
  });
});

describe("focus kept, and where it lands", () => {
  it("a switch pressed keeps focus on its button, and so does a hover over other tracks", async () => {
    const { container, ctx } = await opened();
    const sel = '[data-seqv-key="choice:section.rowSizing:natural"]';
    const before = container.querySelector(sel) as HTMLButtonElement;
    expect(before).not.toBeNull();
    before.focus();
    await click(before);
    await settle();
    const after = container.querySelector(sel) as HTMLButtonElement;
    expect([after === before, after.getAttribute("aria-pressed"), active()]).toEqual([true, "true", after]);
    for (const track of ["t.res", "t.dom", "t.res"])
      await act(async () => pointer(ctx)("move", S, xOf(ctx, S, 12), at(ctx, S, track).y));
    expect(active()).toBe(after);
  });

  it("a switch with a menu in the header opens and closes like the track's; a choice leaves focus on the trigger", async () => {
    const { container, ctx } = await opened();
    const sel = '[data-seqv-key="switch:section.style"]';
    const t = container.querySelector(sel) as HTMLButtonElement;
    expect(t.getAttribute("aria-haspopup")).toBe("menu");
    await click(t);
    await press(active() as HTMLElement, "Escape");
    expect([container.querySelector("[role=menu]"), active()]).toEqual([null, t]);
    await click(t);
    const preset = items(container)[1];
    preset.focus();
    await click(preset);
    await settle();
    expect(ctx.ui.style.value.get(S)).toBe(preset.textContent);
    const after = container.querySelector(sel) as HTMLButtonElement;
    expect([after === t, active(), container.querySelector("[role=menu]")]).toEqual([true, after, null]);
  });

  it('"Show X" moves focus to the next chip, then, with none left, to the surface', async () => {
    const { container, ctx } = await opened();
    await act(async () => {
      await ctx.actions.run("track.hide", "t.res");
      await ctx.actions.run("track.hide", "t.dom");
    });
    await settle();
    const chips = () => [...container.querySelectorAll<HTMLButtonElement>(".seqv-chips button")];
    expect(chips().map((c) => c.textContent)).toEqual(["Show Residues", "Show Domains"]);
    chips()[0].focus();
    await click(chips()[0]);
    await settle();
    expect(chips().map((c) => c.textContent)).toEqual(["Show Domains"]);
    expect(active()).toBe(chips()[0]);
    await click(chips()[0]);
    await settle();
    expect(container.querySelector(".seqv-chips")).toBeNull();
    expect(active()).toBe(surfaceOf(container));
  });

  it("a last chip's neighbor before it takes focus when none follows", async () => {
    const { container, ctx } = await opened();
    await act(async () => {
      await ctx.actions.run("track.hide", "t.res");
      await ctx.actions.run("track.hide", "t.dom");
    });
    await settle();
    const chips = () => [...container.querySelectorAll<HTMLButtonElement>(".seqv-chips button")];
    chips()[1].focus();
    await click(chips()[1]);
    await settle();
    expect([chips().map((c) => c.textContent), active()]).toEqual([["Show Residues"], chips()[0]]);
  });

  it("chips that go in one commit pass focus to none of themselves: it lands on the surface", async () => {
    const { container, ctx } = await opened();
    await act(async () => {
      await ctx.actions.run("track.hide", "t.res");
      await ctx.actions.run("track.hide", "t.dom");
    });
    await settle();
    const chips = () => [...container.querySelectorAll<HTMLButtonElement>(".seqv-chips button")];
    expect(chips().length).toBe(2);
    chips()[0].focus();
    // The document again: both tracks shown, both chips removed together (the first's landing would be the second).
    await act(async () => void ctx.setDocument(uiDoc()));
    await settle();
    expect([container.querySelector(".seqv-chips"), active()]).toEqual([null, surfaceOf(container)]);
  });

  it("Hide from the track's menu leaves focus on the section's surface, not on the page", async () => {
    const { container } = await atDomains();
    await click(trigger(container));
    const hide = items(container)[3];
    expect(hide.textContent).toBe("Hide track");
    hide.focus();
    await click(hide);
    await settle();
    expect(container.querySelector(".seqv-track-actions")).toBeNull();
    expect(container.querySelector(".seqv-chips")?.textContent).toBe("Show Domains");
    expect(active()).toBe(surfaceOf(container));
  });

  it("the player's Next at the last snapshot and Previous at the first stay focusable, say so, and do nothing", async () => {
    const { container, ctx } = await opened();
    const next = container.querySelector('[aria-label="Next snapshot"]') as HTMLButtonElement;
    const prev = container.querySelector('[aria-label="Previous snapshot"]') as HTMLButtonElement;
    expect([prev.getAttribute("aria-disabled"), prev.hasAttribute("disabled")]).toEqual(["true", false]);
    prev.focus();
    await click(prev);
    expect([active(), ctx.snapshots.current.value]).toEqual([prev, null]);
    next.focus();
    await click(next);
    await settle();
    expect([ctx.snapshots.current.value, active()]).toEqual(["s1", next]);
    await click(next);
    await settle();
    expect([ctx.snapshots.current.value, next.getAttribute("aria-disabled"), next.hasAttribute("disabled")]).toEqual([
      "s2",
      "true",
      false,
    ]);
    await click(next);
    await settle();
    expect([ctx.snapshots.current.value, active()]).toEqual(["s2", next]);
  });

  it("an invalid snapshot's dot says so, stays focusable and applies nothing", async () => {
    const { container, ctx } = await opened();
    // s1 hides t.dom: a base patch removing it leaves s1 no longer applying.
    await act(async () => void ctx.applyPatch(patch([{ op: "remove", ref: "t.dom" }])));
    await settle();
    const dot = container.querySelector('[data-seqv-key="dot:s1"]') as HTMLButtonElement;
    expect([dot.getAttribute("data-invalid"), dot.getAttribute("aria-disabled"), dot.hasAttribute("disabled")]).toEqual(
      ["", "true", false],
    );
    dot.focus();
    const issues = ctx.issues.current.value.length;
    await click(dot);
    // Nothing was asked of the player: it would have reported the snapshot it was asked to apply as no longer applying.
    expect([active(), ctx.snapshots.current.value, ctx.issues.current.value.length]).toEqual([dot, null, issues]);
    const ok = container.querySelector('[data-seqv-key="dot:s2"]') as HTMLButtonElement;
    expect(ok.getAttribute("aria-disabled")).toBe("false");
    await click(ok);
    await settle();
    expect(ctx.snapshots.current.value).toBe("s2");
  });

  it("focus inside a section that leaves the view lands on the viewer's root", async () => {
    const { container, ctx } = await opened();
    surfaceOf(container).focus();
    expect(active()).toBe(surfaceOf(container));
    await act(async () => ctx.actions.run("view.switch", "msa"));
    await settle();
    // The chains section's surface is gone with its frame; the alignment's has not taken focus.
    expect(active()).toBe(container.querySelector(".seqv-root"));
  });

  it.each([
    ["ViewTabs", '[data-seqv-key="tab:chains"]', ".seqv-root"],
    ["SnapshotPlayer", '[data-seqv-key="dot:s2"]', ".seqv-root"],
    ["Legend", '[data-seqv-key="legend"]', ".seqv-root"],
    ["Export", '[data-seqv-key="export"]', ".seqv-root"],
    ["Issues", '[data-seqv-key="issues"]', ".seqv-root"],
    ["SectionHeader", '[data-seqv-key="choice:section.rowSizing:natural"]', ".seqv-surface [tabindex]"],
  ] as [SlotName, string, string][])(
    "focus inside %s when React removes it lands on the surface or the root",
    async (name, inside, lands) => {
      // Every part shown: a context with the export action, and an issue for the badge.
      const ctx = createSeqViewer(withDownloadExport(headlessSpec()), {
        scheduler: manualScheduler(),
        text: fixtureMeasurer(),
      });
      cleanups.push(() => ctx.dispose());
      await ctx.setDocument(uiDoc());
      ctx.report("seqv.test", "a reported problem");
      const p = parent((i) => (
        <SeqViewerProvider context={ctx} ui={{ surface: stubSurface, components: i === 0 ? {} : { [name]: "none" } }}>
          <SeqViewerLayout />
        </SeqViewerProvider>
      ));
      await p.render();
      const el = p.container.querySelector<HTMLElement>(inside);
      expect(el).not.toBeNull();
      el?.focus();
      expect(active()).toBe(el);
      await p.render();
      expect(p.container.querySelector(inside)).toBeNull();
      expect(active()).toBe(p.container.querySelector(lands));
    },
  );

  it("the tabs: arrows move between them, the focus stays on the selected one's group", async () => {
    const { container, ctx } = await opened();
    const tabs = () => [...container.querySelectorAll<HTMLButtonElement>("[role=tab]")];
    tabs()[0].focus();
    await press(tabs()[0], "ArrowRight");
    await settle();
    expect([ctx.views.current.value, active()]).toEqual(["msa", tabs()[1]]);
  });
});

describe("the issue list opens like a menu, without roving", () => {
  it("Escape and a press outside close it and focus goes back to the badge; its issues vanishing close it", async () => {
    const { container, ctx } = await opened();
    await act(async () => ctx.report("seqv.test", "a reported problem"));
    const badge = () => container.querySelector(".seqv-badge") as HTMLButtonElement;
    const list = () => container.querySelector(".seqv-issue-list");
    expect(badge().getAttribute("aria-expanded")).toBe("false");
    badge().focus();
    await click(badge());
    expect([badge().getAttribute("aria-expanded"), list()?.textContent]).toEqual([
      "true",
      "Viewerseqv.test a reported problem",
    ]);
    await press(badge(), "Escape");
    expect([list(), badge().getAttribute("aria-expanded"), active()]).toEqual([null, "false", badge()]);
    await click(badge());
    await act(async () => document.body.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true })));
    expect([list(), badge().getAttribute("aria-expanded")]).toEqual([null, "false"]);
    // Open, its issues cleared: the badge goes with the list, and the next issue finds it closed.
    await click(badge());
    expect(list()).not.toBeNull();
    await act(async () => ctx.setDocument(uiDoc()));
    await settle();
    expect([badge, list].map((f) => f())).toEqual([null, null]);
    await act(async () => ctx.report("seqv.test", "another"));
    expect([list(), badge().getAttribute("aria-expanded")]).toEqual([null, "false"]);
  });
});

describe("the document keeps no listener once the UI is gone", () => {
  it("an open menu's pointerdown listener goes with its unmount", async () => {
    const { container, ctx, unmount } = await atDomains();
    await click(trigger(container));
    const add = vi.spyOn(document, "addEventListener");
    const remove = vi.spyOn(document, "removeEventListener");
    await unmount();
    cleanups.splice(0);
    expect(remove.mock.calls.filter(([t]) => t === "pointerdown").length).toBeGreaterThan(0);
    expect(add.mock.calls.filter(([t]) => t === "pointerdown").length).toBe(0);
    ctx.dispose();
    vi.restoreAllMocks();
  });
});

describe("under StrictMode (effects and refs mounted twice)", () => {
  it("the track's menu opens, moves, closes and gives focus back; Hide lands focus on the surface", async () => {
    const { container } = await atDomains({}, true);
    const t = trigger(container);
    await click(t);
    expect([items(container).length, active(), t.getAttribute("aria-expanded")]).toEqual([
      4,
      items(container)[0],
      "true",
    ]);
    await press(active() as HTMLElement, "ArrowDown");
    expect(active()).toBe(items(container)[1]);
    await press(active() as HTMLElement, "Escape");
    expect([container.querySelector("[role=menu]"), active(), t.getAttribute("aria-expanded")]).toEqual([
      null,
      t,
      "false",
    ]);
    await click(t);
    const hide = items(container)[3];
    hide.focus();
    await click(hide);
    await settle();
    expect([container.querySelector(".seqv-track-actions"), active()]).toEqual([null, surfaceOf(container)]);
  });

  it("Shift+F10 opens the cursor's menu; the issue list opens and closes; one pointerdown listener while open", async () => {
    const { container, ctx } = await opened({}, true);
    const key = await cursorOnDomains(ctx);
    const add = vi.spyOn(document, "addEventListener");
    await act(async () => void key(S, "F10", true));
    expect([items(container).length, active()]).toEqual([4, items(container)[0]]);
    expect(add.mock.calls.filter(([t]) => t === "pointerdown").length).toBe(1);
    await press(active() as HTMLElement, "Escape");
    expect(active()?.getAttribute("aria-label")).toBe("Actions for Domains");
    await act(async () => ctx.report("seqv.test", "a reported problem"));
    const badge = container.querySelector(".seqv-badge") as HTMLButtonElement;
    badge.focus();
    await click(badge);
    expect(container.querySelector(".seqv-issue-list")?.textContent).toBe("Viewerseqv.test a reported problem");
    await press(badge, "Escape");
    expect([container.querySelector(".seqv-issue-list"), active()]).toEqual([null, badge]);
    vi.restoreAllMocks();
  });
});

describe("a default Menu whose trigger React replaces", () => {
  it("opens from the new trigger, which carries the ARIA state; the old one does nothing", async () => {
    const items = [{ key: "a", label: "A", run: () => {} }];
    const Owner = ({ n }: { n: number }) => {
      const [open, setOpen] = useState(false);
      const trigger = n === 0 ? <button type="button">open</button> : <a href="#menu">open</a>;
      return <Menu trigger={trigger} items={items} open={open} onOpenChange={setOpen} />;
    };
    const p = parent((i) => <Owner n={i} />, true);
    await p.render();
    const old = p.container.querySelector("button") as HTMLButtonElement;
    await p.render();
    const link = p.container.querySelector("a") as HTMLAnchorElement;
    expect([link.getAttribute("aria-haspopup"), link.getAttribute("aria-expanded")]).toEqual(["menu", "false"]);
    await act(async () => old.click());
    expect(p.container.querySelector("[role=menu]")).toBeNull();
    await act(async () => link.click());
    expect([p.container.querySelector("[role=menu]")?.textContent, link.getAttribute("aria-expanded")]).toEqual([
      "A",
      "true",
    ]);
    await press(active() as HTMLElement, "Escape");
    expect([p.container.querySelector("[role=menu]"), active()]).toEqual([null, link]);
  });
});
