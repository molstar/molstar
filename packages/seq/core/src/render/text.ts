/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TextMeasure } from '@molstar/core/util/text-measure';

// The text metrics port (_spec/README.md "Render", Text): `measure(font, text)` in CSS px, one `TextMeasure` per font
// behind it. A font is CSS shorthand naming a family role (`12px sans`, `bold 10px mono`). The browser adapter
// measures with a canvas (the viewer's); Node measures with the bundled table: an Arial-compatible sans and a
// fixed-advance mono. Letters are centered in their cells and never measured.

export interface TextMetrics {
  /** The width of `text` set in `font`, CSS px. */
  measure(font: string, text: string): number;
  /** Drops every cached width (a web font arrived): the caller lays out again. */
  clear(): void;
  /** Bumped by `clear`: part of the draw lists' key, so a frame after it measures again. */
  readonly version: number;
}

/** A metrics port over a width function per font (a canvas's `measureText` in the browser). */
export function createTextMetrics(measureIn: (font: string) => (text: string) => number): TextMetrics {
  const byFont = new Map<string, TextMeasure>();
  let version = 0;
  return {
    get version() {
      return version;
    },
    measure(font, text) {
      let m = byFont.get(font);
      if (!m) byFont.set(font, (m = new TextMeasure(measureIn(font))));
      return m.width(text);
    },
    clear() {
      for (const m of byFont.values()) m.clear();
      version++;
    },
  };
}

/** A font's size in px (12 when it names none) and whether its family is the mono role. */
export function parseFont(font: string): { readonly size: number; readonly mono: boolean } {
  const m = /(\d+(?:\.\d+)?)px/.exec(font);
  return { size: m ? Number(m[1]) : 12, mono: /\bmono(space)?\b/i.test(font) };
}

// Advances per 1000 em of an Arial-compatible sans (the metric-compatible Liberation Sans) for ASCII 32 to 126.
const SansAdvances = [
  278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278, 556, 556, 556, 556, 556, 556, 556,
  556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833,
  722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 333, 556, 556, 500, 556,
  556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, 334,
  260, 334, 584,
];

/** A character outside ASCII: an average advance. */
const SansDefault = 556;
const MonoAdvance = 600;

/** The bundled table's width of `text` in `font`: no kerning, a fixed advance for the mono role. */
export function bundledWidth(font: string, text: string): number {
  const { size, mono } = parseFont(font);
  let units = 0;
  for (let i = 0; i < text.length; i++) {
    const c = text.charCodeAt(i);
    // a low surrogate continues the character its high one began
    if (c >= 0xdc00 && c <= 0xdfff) continue;
    units += mono ? MonoAdvance : c >= 32 && c <= 126 ? SansAdvances[c - 32] : SansDefault;
  }
  return (units * size) / 1000;
}

/** Node's port: the bundled metrics. */
export function bundledTextMetrics(): TextMetrics {
  return createTextMetrics((font) => (text) => bundledWidth(font, text));
}
