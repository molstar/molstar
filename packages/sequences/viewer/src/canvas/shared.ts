/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a context's surfaces share (seqspec:viewer#the-render-backend, seqspec:viewer#text-measuring,
// seqspec:interaction#the-skin): the canvas text measurer and the font-load listener, set by the first mounted surface
// and restored by the last; the skin reader, one per context, which reads the `--seq-*` tokens from a surface's
// computed style at most once a frame after a change (a theme attribute or class on the viewer's root or an element
// above it, the color scheme, or a host's `refreshSkin`) and publishes the skin it read, which a
// host's `setSkin` holds back until `setSkin(null)`. A font load only invalidates the text caches (the measurer's em
// boxes too, through the text version): it never sets a measurer, so a host's own (`measurer: false`) stays. Also the
// browser scheduler and `canvasEmBox`, which `/dom`'s export re-exports.

import { SeqViewerConfig } from "../config.js";
import { DarkSkin, DefaultSkin, type SeqViewerContext } from "../context.js";
import type { Skin } from "../graph.js";
import { rafScheduler, type Scheduler } from "../scheduler.js";
import type { EmBox, MetricsMeasurer } from "../svg/metrics.js";
import type { TextMeasurer } from "../text.js";

/** `ports.scheduler` for a page: animation frames, timers and `performance.now()` of `window`. */
export function windowScheduler(): Scheduler {
  return rafScheduler({
    requestAnimationFrame: (cb) => window.requestAnimationFrame(cb),
    cancelAnimationFrame: (id) => window.cancelAnimationFrame(id),
    setTimeout: (cb, ms) => window.setTimeout(cb, ms),
    clearTimeout: (id) => window.clearTimeout(id),
    now: () => performance.now(),
  });
}

/**
 * Canvas `measureText` on one offscreen context, its font set only when it changes; `metrics`, the em box `canvasEmBox`
 * reads, so `renderSvg` over a context it measures writes baselines as `dy` (read again when `version` moves).
 */
export function canvasMeasurer(version?: () => number): MetricsMeasurer {
  const c = new OffscreenCanvas(1, 1).getContext("2d") as OffscreenCanvasRenderingContext2D;
  let font = "";
  return {
    measure(f, text) {
      if (f !== font) c.font = font = f;
      return c.measureText(text).width;
    },
    metrics: canvasEmBox(version),
  };
}

/**
 * The em box of a CSS font string as this browser's Canvas2D places `top`, `middle` and `bottom` baselines: where the foot of
 * an "H" lies under each, read at 100 px once per font and weight (hundredths of an em), so a page's export writes the
 * baselines as `dy` exactly where its canvas draws them, whatever fonts it has. With it the ink tops of "H" and "x" above
 * the alphabetic baseline (`cap`, `xHeight`), read at the font string's own size (to a hundredth of a px, from the top
 * row's coverage), once per font string: at a device size they are the hinted heights letters are centered on
 * (`TextManager.ink`). `version` is read on each call: when it moved, every box is read again, since a web font not
 * loaded yet when its box was first read gave its fallback's (the context's text version, which a font load bumps). The
 * probe's canvas is made on the first call.
 */
export function canvasEmBox(version?: () => number): (font: string) => EmBox | null {
  let g: OffscreenCanvasRenderingContext2D | null = null;
  const cache = new Map<string, EmBox | null>();
  const sized = new Map<string, EmBox | null>();
  let at = version ? version() : 0;
  const foot = (p: OffscreenCanvasRenderingContext2D, baseline: CanvasTextBaseline): number => {
    p.clearRect(0, 0, 160, 320);
    p.textBaseline = baseline;
    p.fillText("H", 10, 100);
    const d = p.getImageData(0, 0, 160, 320).data;
    for (let y = 319; y >= 0; y--) for (let x = 0; x < 160; x++) if (d[(y * 160 + x) * 4 + 3] > 127) return y + 1;
    return 0;
  };
  /**
   * The ink top of `ch` above an alphabetic baseline at y = 100 (px): its top row inks the share of a px its coverage is of
   * its column's most (a stem a px wide covers less than a px across, whatever its height).
   */
  const head = (p: OffscreenCanvasRenderingContext2D, ch: string): number => {
    p.clearRect(0, 0, 160, 320);
    p.textBaseline = "alphabetic";
    p.fillText(ch, 10, 100);
    const d = p.getImageData(0, 0, 160, 100).data;
    for (let y = 0; y < 100; y++) {
      let share = 0;
      for (let x = 0; x < 160; x++) {
        const a = d[(y * 160 + x) * 4 + 3];
        if (a === 0) continue;
        let most = a;
        for (let r = y + 1; r < 100; r++) most = Math.max(most, d[(r * 160 + x) * 4 + 3]);
        share = Math.max(share, a / most);
      }
      if (share > 0) return 99 - y + share;
    }
    return 0;
  };
  return (font) => {
    if (version && version() !== at) {
      at = version();
      cache.clear();
      sized.clear();
    }
    let full = sized.get(font);
    if (full !== undefined) return full;
    const probe = font.replace(/[\d.]+px/, "100px");
    let box = cache.get(probe);
    if (box === undefined) {
      g ??= new OffscreenCanvas(160, 320).getContext("2d", {
        willReadFrequently: true,
      }) as OffscreenCanvasRenderingContext2D;
      g.font = probe;
      const alphabetic = foot(g, "alphabetic");
      const [top, bottom] = [foot(g, "top"), foot(g, "bottom")];
      box = alphabetic > 0 ? { ascent: (top - alphabetic) / 100, descent: (alphabetic - bottom) / 100 } : null;
      cache.set(probe, box);
    }
    full = null;
    if (box) {
      // At the font's own size up to 100 px (the probe's room), else the probe's heights scaled.
      const px = Number(/([\d.]+)px/.exec(font)?.[1] ?? 0);
      const p = g as OffscreenCanvasRenderingContext2D;
      p.font = px > 0 && px <= 100 ? font : probe;
      const em = px > 0 && px <= 100 ? px : 100;
      full = { ...box, cap: head(p, "H") / em, xHeight: head(p, "x") / em };
    }
    sized.set(font, full);
    return full;
  };
}

const TOKENS: readonly (keyof Skin)[] = [
  "surface",
  "surfaceAlt",
  "ink",
  "inkQuiet",
  "inkFaint",
  "grid",
  "accent",
  "hover",
  "highlight",
  "selection",
  "selectionBorder",
  "unobserved",
  "pending",
  "grip",
  "gripHover",
  "window",
  "conflict",
  "focusRing",
];
const cssName = (k: string) => `--seq-${k.replace(/[A-Z]/g, (c) => `-${c.toLowerCase()}`)}`;
/** Every property the reader reads, named once. */
const NAMES: readonly string[] = [...TOKENS.map(cssName), "--seq-font-mono", "--seq-font-ui"];

/** WCAG 2's relative luminance of an sRGB channel (0..255). */
const linear = (c: number) => (c <= 10.31475 ? c / 3294.6 : ((c / 255 + 0.055) / 1.055) ** 2.4);
/**
 * Below it a surface reads dark (light ink): where light and dark text contrast equally on it, WCAG luminance 0.179, as
 * Mol* Prism's canvas theme decides.
 */
const DARK_LUMINANCE = 0.179;

/** Any CSS color (oklch, color-mix, var-substituted) as `[rgb(…), WCAG luminance 0..1]`, through a 1×1 canvas; edge-rate. */
function normalizer() {
  const c = new OffscreenCanvas(1, 1).getContext("2d", {
    willReadFrequently: true,
  }) as OffscreenCanvasRenderingContext2D;
  return (css: string): [string, number] => {
    c.clearRect(0, 0, 1, 1);
    c.fillStyle = "#000";
    c.fillStyle = css;
    c.fillRect(0, 0, 1, 1);
    const [r, g, b, a] = c.getImageData(0, 0, 1, 1).data;
    // A transparent color has no luminance of its own: NaN, so it decides nothing.
    const y = a === 0 ? Number.NaN : 0.2126 * linear(r) + 0.7152 * linear(g) + 0.0722 * linear(b);
    return [a === 255 ? `rgb(${r} ${g} ${b})` : `rgb(${r} ${g} ${b} / ${(a / 255).toFixed(3)})`, y];
  };
}

/** The page's own theme: the root's `data-theme`, else Prism's `data-msp-theme` (blue is dark), else the color scheme. */
function pageDark(): boolean {
  const { theme, mspTheme } = document.documentElement.dataset;
  if (theme === "dark" || theme === "light") return theme === "dark";
  if (mspTheme && mspTheme !== "auto") return mspTheme !== "light";
  return matchMedia("(prefers-color-scheme: dark)").matches;
}

/**
 * The skin at `el`: its tokens over the light or dark defaults (dark: a surface below `DARK_LUMINANCE`, or, with no
 * surface token, the page's theme).
 */
export function readSkin(el: Element, norm = normalizer()): Skin {
  const style = getComputedStyle(el);
  const surface = style.getPropertyValue("--seq-surface").trim();
  const y = surface ? norm(surface)[1] : Number.NaN;
  const dark = Number.isNaN(y) ? pageDark() : y < DARK_LUMINANCE;
  const out: Record<string, unknown> = { ...(dark ? DarkSkin : DefaultSkin), dark };
  for (const k of TOKENS) {
    const v = style.getPropertyValue(cssName(k)).trim();
    if (v) out[k] = norm(v)[0];
  }
  const mono = style.getPropertyValue("--seq-font-mono").trim();
  const ui = style.getPropertyValue("--seq-font-ui").trim();
  if (mono) out.fontMono = mono;
  if (ui) out.fontUi = ui;
  return out as unknown as Skin;
}

/**
 * The viewer's root (a host's `ui.className` on it) and the elements above it (crossing a shadow root to its host), up
 * to the document's root.
 */
function hostAncestors(el: Element): Element[] {
  const out: Element[] = [];
  let n: Node | null = el.closest(".seqv-root") ?? el.parentNode;
  while (n) {
    if (n instanceof Element) out.push(n);
    n = n instanceof ShadowRoot ? n.host : n.parentNode;
  }
  if (!out.includes(document.documentElement)) out.push(document.documentElement);
  return out;
}

const same = (a: Skin | null, b: Skin) => a !== null && (Object.keys(b) as (keyof Skin)[]).every((k) => a[k] === b[k]);

interface Shared {
  surfaces: Set<Element>;
  previous: TextMeasurer;
  measurer: TextMeasurer | null;
  /** A read on the next frame (one, however many asked). */
  schedule: () => void;
  stop: () => void;
}
const shared = new WeakMap<SeqViewerContext, Shared>();

/**
 * A surface joins its context's shared parts: the first sets the canvas measurer (unless `measurer: false`) and starts
 * the skin reader and the font listener; the returned leave stops them with the last, restoring the measurer.
 */
/**
 * Reads the skin again on the next frame (seqspec:interaction#the-skin): for a host that changes the canvas's tokens in
 * a way the reader does not watch (an inline `--seq-*` or color on an element above the viewer, a stylesheet swapped).
 * Nothing while no surface of the context is mounted.
 */
export function refreshSkin(ctx: SeqViewerContext): void {
  shared.get(ctx)?.schedule();
}

export function joinShared(ctx: SeqViewerContext, el: Element, opts: { measurer?: boolean } = {}): () => void {
  let s = shared.get(ctx);
  if (!s) {
    // Its em boxes read again after a font load (`fonts` below bumps the text version).
    const measurer = opts.measurer === false ? null : canvasMeasurer(() => ctx.text.version.value);
    const norm = normalizer();
    let last: Skin | null = null;
    let raw = "";
    const read = () => {
      frame = 0;
      const first = (s as Shared).surfaces.values().next().value as Element | undefined;
      if (!first || ctx.disposed) return;
      // The tokens as written, and the page's theme: the same as last time, nothing is normalized or published.
      const style = getComputedStyle(first);
      let now = pageDark() ? "dark" : "light";
      for (const n of NAMES) now += `|${style.getPropertyValue(n)}`;
      if (now === raw) return;
      raw = now;
      const skin = readSkin(first, norm);
      if (!same(last, skin)) {
        last = skin;
        ctx.publishSkin(skin);
      }
      ctx.text.setFamilies(skin.fontMono, skin.fontUi);
    };
    let frame = 0;
    const schedule = () => {
      if (!frame && !ctx.disposed) frame = ctx.scheduler.requestFrame(read);
    };
    // The page's theme on its root, and a host's canvas between it and the viewer (its own theme or class on an
    // ancestor, `ui.className` on the root): read when it joins. Not `style`: a host writes inline properties on its
    // layout during drags (a panel's width), so a color set inline asks for a read with `refreshSkin`.
    const mo = new MutationObserver(schedule);
    const watch = { attributes: true, attributeFilter: ["data-theme", "data-msp-theme", "class"] };
    for (const a of hostAncestors(el)) mo.observe(a, watch);
    const media = matchMedia("(prefers-color-scheme: dark)");
    media.addEventListener("change", schedule);
    // A web font that loads after text was measured: every prepared text is stale (pretext's `clearCache`).
    const fonts = () => !ctx.disposed && ctx.text.invalidate();
    document.fonts.addEventListener("loadingdone", fonts);
    s = {
      surfaces: new Set(),
      previous: ctx.text.port,
      measurer,
      schedule,
      stop: () => {
        mo.disconnect();
        if (frame) ctx.scheduler.cancelFrame(frame);
        media.removeEventListener("change", schedule);
        document.fonts.removeEventListener("loadingdone", fonts);
      },
    };
    shared.set(ctx, s);
    s.surfaces.add(el);
    // `Viewport.ReducedMotion` null: the page's preference, read once (seqspec:interaction#accessibility).
    const motion = SeqViewerConfig.Viewport.ReducedMotion;
    if (ctx.config.get(motion) === null && matchMedia("(prefers-reduced-motion: reduce)").matches)
      ctx.config.set(motion, true);
    if (measurer) ctx.setTextMeasurer(measurer);
    read();
  } else s.surfaces.add(el);
  const held = s;
  return () => {
    held.surfaces.delete(el);
    if (held.surfaces.size > 0) return;
    held.stop();
    shared.delete(ctx);
    if (!ctx.disposed && held.measurer && ctx.text.port === held.measurer) ctx.setTextMeasurer(held.previous);
  };
}
