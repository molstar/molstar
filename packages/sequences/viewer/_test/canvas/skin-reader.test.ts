/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The skin reader and the measurer a surface joins (CS2; seqspec:viewer#the-render-backend,
// seqspec:viewer#text-measuring, seqspec:interaction#the-skin), in Chromium: the tokens read from the host's computed
// style over the light or dark defaults (the dark ones complete, with their contrast), a theme toggle republished, a
// host's `setSkin` kept through it and `setSkin(null)` showing the skin last read, the measurer opt-out and a font load
// that only invalidates, the first surface setting the canvas measurer and the last one restoring the context's.

import { afterEach, describe, expect, it, vi } from "vitest";
import { refreshSkin } from "../../src/canvas/index.js";
import { DefaultSkin, type Skin } from "../../src/index.js";
import { syntheticFamily } from "../families.js";
import { canvasesOf, close, diff, type Opened, opened, pixelsOf } from "./open.js";

let t: Opened | null = null;
afterEach(() => {
  if (t) close(t);
  t = null;
  delete document.documentElement.dataset.theme;
  delete document.documentElement.dataset.mspTheme;
  document.body.style.removeProperty("--seq-surface");
});
const doc = () => syntheticFamily(10, 6);
/** A task later: the MutationObserver's callback has run, then a frame: the reader's read. */
const later = async () => {
  await new Promise((r) => setTimeout(r, 0));
  t?.advance();
};

describe("the skin reader", () => {
  it("reads the tokens from the host's computed style over the light defaults, any CSS color as rgb()", async () => {
    t = await opened(doc(), {
      style: "--seq-surface:#fff;--seq-ink:oklch(0.3 0 0);--seq-accent:rgb(1 2 3);--seq-font-ui:Arial",
    });
    const s = t.ctx.skin.value;
    expect(s.dark).toBe(false);
    expect(s.surface).toBe("rgb(255 255 255)");
    expect(s.accent).toBe("rgb(1 2 3)");
    expect(s.ink).toMatch(/^rgb\(\d+ \d+ \d+\)$/);
    expect(s.inkQuiet).toBe(DefaultSkin.inkQuiet);
    expect(s.fontUi).toBe("Arial");
    expect(t.ctx.text.fontString(t.ctx.text.font("ui", 12))).toContain("Arial");
  });

  it("a dark surface reads dark, over the complete dark defaults: highlight, conflict and the focus ring too", async () => {
    t = await opened(doc(), { style: "--seq-surface:rgb(20 22 26)" });
    const s = t.ctx.skin.value;
    expect(s.dark).toBe(true);
    expect([s.highlight, s.conflict, s.focusRing]).toEqual([
      "rgb(250 204 21 / 0.28)",
      "rgb(255 123 114)",
      "rgb(110 160 255)",
    ]);
    expect(s.surface).toBe("rgb(20 22 26)");
  });

  it("the root's data-theme and the color scheme republish the skin", async () => {
    t = await opened(doc());
    expect(t.ctx.skin.value.dark).toBe(false);
    document.documentElement.dataset.theme = "dark";
    await later();
    expect(t.ctx.skin.value.dark).toBe(true);
    document.documentElement.dataset.theme = "light";
    await later();
    expect(t.ctx.skin.value).toEqual(expect.objectContaining({ dark: false, surface: DefaultSkin.surface }));
  });

  it("dark is a surface below WCAG luminance 0.179 (Prism's canvas rule), not a channel average", async () => {
    // A mid green averages dark (0.21) but reads light: dark ink contrasts more on it.
    t = await opened(doc(), { style: "--seq-surface:rgb(0 160 0)" });
    expect(t.ctx.skin.value.dark).toBe(false);
    close(t);
    t = await opened(doc(), { style: "--seq-surface:rgb(0 110 0)" });
    expect(t.ctx.skin.value.dark).toBe(true);
  });

  it("Prism's `data-msp-theme` on the root republishes (blue is dark), as does a host's canvas color on an ancestor", async () => {
    t = await opened(doc());
    document.documentElement.dataset.mspTheme = "blue";
    await later();
    expect(t.ctx.skin.value.dark).toBe(true);
    document.documentElement.dataset.mspTheme = "light";
    await later();
    expect(t.ctx.skin.value.dark).toBe(false);
    // A canvas of the host's own color, set inline above the viewer after it mounted: unwatched (a host's drags write
    // inline properties), read on the host's `refreshSkin`.
    document.body.style.setProperty("--seq-surface", "rgb(20 22 26)");
    await later();
    expect(t.ctx.skin.value.dark).toBe(false);
    refreshSkin(t.ctx);
    await later();
    expect([t.ctx.skin.value.dark, t.ctx.skin.value.surface]).toEqual([true, "rgb(20 22 26)"]);
  });

  it("a class on the viewer's root (React's `ui.className`) or above it republishes; a transparent surface decides nothing", async () => {
    const sheet = new CSSStyleSheet();
    sheet.replaceSync(".night { --seq-surface: rgb(20 22 26) } .clear { --seq-surface: transparent }");
    document.adoptedStyleSheets = [...document.adoptedStyleSheets, sheet];
    const root = document.createElement("div");
    root.className = "seqv-root";
    document.body.append(root);
    try {
      t = await opened(doc(), { parent: root });
      expect(t.ctx.skin.value.dark).toBe(false);
      root.classList.add("night");
      await later();
      expect(t.ctx.skin.value.dark).toBe(true);
      root.classList.replace("night", "clear");
      document.documentElement.dataset.theme = "light";
      await later();
      expect(t.ctx.skin.value.dark).toBe(false);
    } finally {
      document.adoptedStyleSheets = document.adoptedStyleSheets.filter((s) => s !== sheet);
      if (t) close(t);
      t = null;
      root.remove();
    }
  });

  it("the mutation path does no work for what it does not watch, and one read a frame for what it does", async () => {
    t = await opened(doc());
    const spy = vi.spyOn(window, "getComputedStyle");
    try {
      // A host's drag: inline properties on an element above the viewer, every pointer move.
      for (let i = 0; i < 200; i++) {
        document.body.style.setProperty("--host-panel-width", `${300 + i}px`);
        await Promise.resolve();
      }
      await later();
      expect(spy).toHaveBeenCalledTimes(0);
      // Fifty class changes in one frame: one read (its computed style, and the tokens' when they changed).
      for (let i = 0; i < 50; i++) document.body.classList.toggle("busy");
      await later();
      expect(spy.mock.calls.length).toBeGreaterThan(0);
      expect(spy.mock.calls.length).toBeLessThanOrEqual(2);
    } finally {
      spy.mockRestore();
      document.body.style.removeProperty("--host-panel-width");
      document.body.classList.remove("busy");
    }
  });

  it("a host's setSkin survives a theme toggle; setSkin(null) shows the skin the reader last read", async () => {
    t = await opened(doc());
    const mine: Skin = { ...DefaultSkin, accent: "rgb(9 9 9)" };
    t.ctx.setSkin(mine);
    document.documentElement.dataset.theme = "dark";
    await later();
    expect(t.ctx.skin.value).toBe(mine);
    t.ctx.setSkin(null);
    expect(t.ctx.skin.value.dark).toBe(true);
    expect(t.ctx.skin.value.surface).toBe("rgb(24 26 30)");
    // Held again by a second override, and the reader still reading underneath it.
    t.ctx.setSkin(mine);
    document.documentElement.dataset.theme = "light";
    await later();
    expect(t.ctx.skin.value).toBe(mine);
    t.ctx.setSkin(null);
    expect(t.ctx.skin.value.dark).toBe(false);
  });

  it("a skin change redraws the passes: the main canvas takes the dark skin's ink", async () => {
    t = await opened(doc());
    const o = t;
    const [main] = canvasesOf(o);
    const light = pixelsOf(main);
    document.documentElement.dataset.theme = "dark";
    await later();
    o.advance();
    const dark = pixelsOf(main);
    expect(diff(light, dark, main.width).count).toBeGreaterThan(200);
    // The ruler's text is ink: the dark ink's color is among the pixels drawn.
    const inks = new Set<string>();
    for (let i = 0; i < dark.length; i += 4)
      if (dark[i + 3] === 255) inks.add(`${dark[i]} ${dark[i + 1]} ${dark[i + 2]}`);
    expect(inks.has("230 232 235")).toBe(true);
  });
});

describe("the measurer and font loads", () => {
  it("the first surface sets the canvas measurer, the last one disposed restores the context's", async () => {
    t = await opened(doc());
    const o = t;
    const canvas = o.ctx.text.port;
    // A second surface joins the same context.
    const host = document.createElement("div");
    host.style.cssText = "width:400px;height:200px";
    document.body.appendChild(host);
    const { mountSection } = await import("../../src/canvas/index.js");
    const second = mountSection(o.ctx, o.section, host, { dpr: 1 });
    expect(o.ctx.text.port).toBe(canvas);
    o.surface.dispose();
    expect(o.ctx.text.port).toBe(canvas);
    second.dispose();
    host.remove();
    expect(o.ctx.text.port).not.toBe(canvas);
    expect(o.ctx.text.port.measure("12px Arial", "WWWW")).toBeGreaterThan(0);
  });

  it("measurer: false keeps the host's measurer, and a font load only invalidates", async () => {
    t = await opened(doc(), { measurer: false });
    const o = t;
    const own = o.ctx.text.port;
    const version = o.ctx.text.version.value;
    document.fonts.dispatchEvent(new Event("loadingdone"));
    expect(o.ctx.text.version.value).toBe(version + 1);
    expect(o.ctx.text.port).toBe(own);
  });

  it("a font load with the canvas measurer invalidates and sets no measurer either", async () => {
    t = await opened(doc());
    const o = t;
    const canvas = o.ctx.text.port;
    const version = o.ctx.text.version.value;
    document.fonts.dispatchEvent(new Event("loadingdone"));
    expect(o.ctx.text.version.value).toBe(version + 1);
    expect(o.ctx.text.port).toBe(canvas);
  });

  it("a surface disposed stops listening for font loads", async () => {
    t = await opened(doc());
    const o = t;
    o.surface.dispose();
    const version = o.ctx.text.version.value;
    document.fonts.dispatchEvent(new Event("loadingdone"));
    expect(o.ctx.text.version.value).toBe(version);
  });
});
