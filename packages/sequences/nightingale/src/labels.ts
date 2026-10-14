/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A track's view and feature labels (seqspec:libraries#nightingale): Nightingale's tracks draw no text, so `mount` lays
// its own labels over a track, a block's name on the block-label track and a range's label inside it, and hovers every
// row by position. Platform-free: the display range a track shows on its pixels both ways, where a label goes, and the
// ink that reads on its fill.

/** What a track shows: its display range (1-based, inclusive) on its viewport [left, right) in pixels. */
export interface NightingaleAxis {
  start: number;
  end: number;
  left: number;
  right: number;
}

/**
 * A track's view from Nightingale's own state, as its `xScale` maps it (withZoom): `display-end` -1 or unset is the
 * length; the viewport is the width less the margins. Written into `out` (updated every zoom frame).
 */
export function axisView(
  out: NightingaleAxis,
  displayStart: number | undefined,
  displayEnd: number | undefined,
  length: number,
  width: number,
  marginLeft: number,
  marginRight: number,
): NightingaleAxis {
  out.start = displayStart && displayStart > 0 ? displayStart : 1;
  out.end = Math.max(out.start, displayEnd === undefined || displayEnd === -1 ? length : displayEnd);
  out.left = marginLeft;
  out.right = Math.max(marginLeft, width - marginRight);
  return out;
}

/** The x of a position's left edge (1-based); `end + 1` is the right edge of `end`. */
export const xOf = (v: NightingaleAxis, position: number): number =>
  v.left + ((position - v.start) * (v.right - v.left)) / (v.end + 1 - v.start);

/** The axis position (1-based) at an x on the track, or -1 off the viewport (a gutter) or the shown range. */
export function positionAt(v: NightingaleAxis, x: number): number {
  if (!(x >= v.left && x < v.right)) return -1;
  const p = Math.floor(v.start + ((x - v.left) * (v.end + 1 - v.start)) / (v.right - v.left));
  return p >= v.start && p <= v.end ? p : -1;
}

/** A label's box: the visible span of its widest fragment, its text centered in it; `none`: hidden. */
export interface LabelBox {
  mode: "full" | "cut" | "none";
  left: number;
  width: number;
}
/** Room around the text (both sides together), the characters a cut label keeps at least, the lowest feature that
 * takes text. */
const LABEL_PAD = 6;
const CUT_CHARS = 3;
const MIN_HEIGHT = 9;

/** The shortest a label is cut to: its first `CUT_CHARS` characters and "…"; null when it is no longer than that. */
export const cutText = (text: string): string | null =>
  [...text].length > CUT_CHARS ? `${[...text].slice(0, CUT_CHARS).join("")}…` : null;

/**
 * Where a feature's label goes at a view: in the widest visible span of its fragments (1-based, inclusive). Written
 * into `out` (each label keeps its box: this runs on every zoom frame). The fit comes as numbers: `text`, the text's
 * width in pixels; `cut`, the narrowest cut text's (its first characters and an ellipsis; at least `text` where none
 * is); `height`, the feature's drawn height (below `MIN_HEIGHT` no label is shown). No object a call, and `xOf`'s
 * arithmetic inline, so no double crosses a helper call (NG-12).
 */
export function placeLabel(
  out: LabelBox,
  v: NightingaleAxis,
  fragments: readonly { start: number; end: number }[],
  text: number,
  cut: number,
  height: number,
): LabelBox {
  let left = 0;
  let width = -1;
  const scale = (v.right - v.left) / (v.end + 1 - v.start);
  for (let k = 0; k < fragments.length; k++) {
    const f = fragments[k];
    const a = Math.max(v.left, v.left + (f.start - v.start) * scale);
    const b = Math.min(v.right, v.left + (f.end + 1 - v.start) * scale);
    if (b - a > width) {
      left = a;
      width = b - a;
    }
  }
  out.left = left;
  out.width = Math.max(0, width);
  if (width <= 0 || height < MIN_HEIGHT) out.mode = "none";
  else if (text + LABEL_PAD <= width) out.mode = "full";
  else out.mode = cut + LABEL_PAD <= width ? "cut" : "none";
  return out;
}

/** sRGB relative luminance of `#rgb` or `#rrggbb`; null for anything else. */
function luminance(hex: string): number | null {
  const m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return null;
  const h = m[1].length === 3 ? [...m[1]].map((c) => c + c).join("") : m[1];
  const [r, g, b] = [0, 2, 4].map((k) => {
    const c = Number.parseInt(h.slice(k, k + 2), 16) / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

/** Dark and light ink for labels on Nightingale's white ground. */
export const INK = { dark: "#1f1f1f", light: "#ffffff" } as const;

/**
 * The ink that reads on a feature's fill (Nightingale's: `fill`, else `color`, else gray, at `opacity` over white),
 * whichever of dark and light has the higher contrast.
 */
export function inkOn(fill: string | undefined, opacity = 0.9): string {
  const l = luminance(fill ?? "") ?? luminance("#808080") ?? 0;
  const seen = opacity * l + (1 - opacity); // over white, in linear light
  const [dark, light] = [luminance(INK.dark) ?? 0, 1];
  return (seen + 0.05) / (dark + 0.05) >= (light + 0.05) / (seen + 0.05) ? INK.dark : INK.light;
}

/** The style properties of a laid label that `layLabels` writes (a `CSSStyleDeclaration` has them). */
export interface StyleLike {
  display: string;
  left: string;
  width: string;
  top: string;
  height: string;
  lineHeight: string;
}
/** A track's layout, as Nightingale's feature tracks keep it. */
export interface LayoutLike {
  getFeatureYPos(f: unknown): number;
  getFeatureHeight(f?: unknown): number;
}
/** One label laid over a feature: its text's width (measured once), its box, and what was last written to its style. */
export interface Laid {
  readonly feature: unknown;
  readonly text: string;
  readonly fragments: readonly { start: number; end: number }[];
  readonly style: StyleLike;
  /** Measured with `measure`, -1 until then; `cut` is the narrowest cut's. */
  width: number;
  cut: number;
  readonly box: LabelBox;
  /** The last written row (top, height) and visibility: a label is written again only when one of them, or its box, changed. */
  top: number;
  height: number;
  visible: boolean;
}

/**
 * Lays labels at a view (every zoom frame, and the manager's attribute storm on each highlight): measures a label once
 * (`measure` reads the computed font only then), places its box, and writes its style only where the box, its row or
 * its height changed. Returns the labels written.
 */
export function layLabels(
  labels: readonly Laid[],
  view: NightingaleAxis,
  layout: LayoutLike | undefined,
  measure: (text: string) => number,
): number {
  let writes = 0;
  for (let i = 0; i < labels.length; i++) {
    const l = labels[i];
    if (l.width < 0) {
      l.width = measure(l.text);
      const cut = cutText(l.text);
      l.cut = cut ? measure(cut) : l.width;
    }
    const height = layout ? layout.getFeatureHeight(l.feature) : 0;
    const left = l.box.left;
    const width = l.box.width;
    placeLabel(l.box, view, l.fragments, l.width, l.cut, height);
    const visible = l.box.mode !== "none";
    const top = visible && layout ? layout.getFeatureYPos(l.feature) : 0;
    const same = l.box.left === left && l.box.width === width && top === l.top && height === l.height;
    if (visible === l.visible && (!visible || same)) continue;
    writes++;
    l.visible = visible;
    const style = l.style;
    if (!visible) {
      style.display = "none";
      continue;
    }
    l.top = top;
    l.height = height;
    style.display = "";
    style.left = `${l.box.left}px`;
    style.width = `${l.box.width}px`;
    style.top = `${top}px`;
    style.height = style.lineHeight = `${height}px`;
  }
  return writes;
}

/** What `onDimensionsChange` is on a Nightingale element (withDimensions): called after its width or height changed. */
export interface DimensionsLike {
  onDimensionsChange?(): void;
}
/**
 * Calls `changed` whenever the element has taken a new width or height (its own `onDimensionsChange`, run after
 * `width` and `height` are set, by its ResizeObserver and by an update that changed them). The track's `width` is a
 * property, no attribute or event announces it, and it lags the layout (first drawn at the width `mount` gave it, then
 * set from its measured size in a later frame): a view read from it before is stale. Returns the way to stop.
 */
export function watchDimensions(node: DimensionsLike, changed: () => void): () => void {
  const own = Object.hasOwn(node, "onDimensionsChange");
  const was = node.onDimensionsChange;
  node.onDimensionsChange = function (this: DimensionsLike) {
    was?.call(this);
    changed();
  };
  return () => {
    if (own) node.onDimensionsChange = was;
    else delete node.onDimensionsChange;
  };
}
