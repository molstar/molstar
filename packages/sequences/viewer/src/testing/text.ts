/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `/testing`: the conformance kit a renderer runs (seqspec:viewer#extension-points) and the fixtures tests and
// benchmarks share. The fixture text measurer (LY1): fixed widths per character, so layouts are deterministic in Node
// (seqspec:viewer#text-measuring).

import { fontPx, type TextMeasurer } from "../text.ts";

/** Em widths: a monospace font's every character; else by class. Combining marks and joiners are zero-width. */
export function fixtureWidth(code: number, mono: boolean): number {
  if ((code >= 0x300 && code < 0x370) || code === 0x200d || (code >= 0xfe00 && code < 0xfe10)) return 0;
  if (mono) return 0.6;
  if (code >= 0x1f000) return 1.2;
  if (code === 32) return 0.25;
  if (code >= 48 && code <= 57) return 0.55;
  if (code >= 65 && code <= 90) return 0.65;
  if (code >= 97 && code <= 122) return 0.5;
  return code < 128 ? 0.35 : 0.8;
}

/** Fixed widths per code point, times the font's px; counts its calls. */
export function fixtureMeasurer(): TextMeasurer & { calls: number } {
  const m = {
    calls: 0,
    measure(font: string, text: string): number {
      m.calls++;
      const px = fontPx(font);
      const mono = font.includes("mono");
      let w = 0;
      for (const ch of text) w += fixtureWidth(ch.codePointAt(0) as number, mono);
      return w * px;
    },
  };
  return m;
}
