/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `metricsMeasurer` (seqspec:viewer#the-render-backend, seqspec:viewer#text-measuring): a `TextMeasurer` over committed
// advance tables, pure data, so layout and `renderSvg` run in Node with widths close to a browser's and no font parser.
// It also answers `metrics(font)`, the em box that Canvas2D's `top`, `middle` and `bottom` baselines lie on: the SVG
// backend writes those baselines as an explicit `dy` from it (every SVG reader honors `dy`; Inkscape ignores
// `dominant-baseline`). Platform-free.

import type { TextMeasurer } from "../text.ts";

/** One family and weight: advances of ASCII 32..126 in `unitsPerEm`, and the em box about the alphabetic baseline. */
export interface FontMetrics {
  /** A CSS family name the font stacks start with (`ui-sans-serif`), or a generic one (`sans-serif`), matched without case. */
  readonly family: string;
  readonly weight: number;
  readonly unitsPerEm: number;
  readonly advances: Uint16Array;
  /** The advance of any other code point. */
  readonly fallback: number;
  /** The em box's top above the alphabetic baseline, as Canvas2D's `textBaseline = "top"` places it. */
  readonly ascent: number;
  /** The em box's bottom below it (`"bottom"`). */
  readonly descent: number;
  /** The ink top of a capital ("H") and of a lower-case letter ("x") above the baseline: letters are centered on them. */
  readonly cap?: number;
  readonly xHeight?: number;
}
/**
 * The em box in em: its top above the alphabetic baseline and its bottom below it; and the ink top of a capital and of a
 * lower-case letter above that baseline (`TextManager.ink`; absent: 0.72 and 0.54 em), at the font string's own size
 * where the measurer reads them from a raster (`canvasEmBox`: hinting included).
 */
export interface EmBox {
  readonly ascent: number;
  readonly descent: number;
  readonly cap?: number;
  readonly xHeight?: number;
}
/** A measurer that also knows the em box of a font string; `renderSvg` reads it off the context's measurer. */
export interface MetricsMeasurer extends TextMeasurer {
  metrics(font: string): EmBox | null;
}
export const hasMetrics = (m: TextMeasurer): m is MetricsMeasurer =>
  typeof (m as MetricsMeasurer).metrics === "function";

/** A CSS font string as the text manager makes them: weight, px size, families. */
export const CSS_FONT = /^(?:(?:italic|oblique|normal)\s+)?(\d+|bold|normal)\s+([\d.]+)px\s+(.*)$/;
const MONO = /^(?:ui-)?monospace$|^sf ?mono|^menlo$|^monaco$|^consolas$|^courier|mono$/;
const SERIF = /^(?:serif|times|times new roman|georgia|palatino|cambria)$/;
const SANS =
  /^(?:sans-serif|system-ui|ui-sans-serif|-apple-system|blinkmacsystemfont|helvetica|arial|segoe ui|roboto)$/;

/** The names a family token may be found under, most specific first: itself, the system font it names, its generic. */
function names(token: string, out: string[]): void {
  const t = token
    .trim()
    .replace(/^["']|["']$/g, "")
    .toLowerCase();
  out.push(t);
  if (t === "system-ui" || t === "-apple-system" || t === "blinkmacsystemfont") out.push("ui-sans-serif");
  if (t.startsWith("sfmono") || t === "sf mono" || t === "menlo") out.push("ui-monospace");
  out.push(MONO.test(t) ? "monospace" : SERIF.test(t) ? "serif" : SANS.test(t) ? "sans-serif" : t);
}

interface Chosen {
  readonly t: FontMetrics;
  readonly px: number;
  readonly box: EmBox;
}

/** A measurer over `tables`: the table of the font's first family a table has (the usual stacks' system and generic names, else
 * the first table's family), the nearest weight; widths are the advances' sum, no kerning, other code points `fallback`. */
export function metricsMeasurer(tables: readonly FontMetrics[]): MetricsMeasurer {
  if (tables.length === 0) throw new Error("metricsMeasurer: no tables.");
  const cache = new Map<string, Chosen | null>();
  const choose = (font: string): Chosen | null => {
    let c = cache.get(font);
    if (c !== undefined) return c;
    const m = CSS_FONT.exec(font);
    if (!m) c = null;
    else {
      const w = m[1] === "bold" ? 700 : m[1] === "normal" ? 400 : Number(m[1]);
      const list: string[] = [];
      for (const token of m[3].split(",")) names(token, list);
      // The first family a table has; a stack with none of them is read as the first table's.
      const fam = (
        list.find((n) => tables.some((t) => t.family.toLowerCase() === n)) ?? tables[0].family
      ).toLowerCase();
      let best: FontMetrics | null = null;
      for (const t of tables)
        if (t.family.toLowerCase() === fam && (!best || Math.abs(t.weight - w) < Math.abs(best.weight - w))) best = t;
      c = best
        ? {
            t: best,
            px: Number(m[2]),
            box: {
              ascent: best.ascent / best.unitsPerEm,
              descent: best.descent / best.unitsPerEm,
              cap: best.cap === undefined ? undefined : best.cap / best.unitsPerEm,
              xHeight: best.xHeight === undefined ? undefined : best.xHeight / best.unitsPerEm,
            },
          }
        : null;
    }
    cache.set(font, c);
    return c;
  };
  return {
    measure(font, text) {
      const c = choose(font);
      if (!c) return 0;
      const { advances, fallback } = c.t;
      let sum = 0;
      for (let i = 0; i < text.length; i++) {
        const code = text.charCodeAt(i);
        if ((code & 0xfc00) === 0xdc00) continue; // a low surrogate: its pair was counted once
        sum += code >= 32 && code < 127 ? advances[code - 32] : fallback;
      }
      return (sum / c.t.unitsPerEm) * c.px;
    },
    metrics(font) {
      return choose(font)?.box ?? null;
    },
  };
}
