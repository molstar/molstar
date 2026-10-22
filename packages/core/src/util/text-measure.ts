/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Text measured for glyph-by-glyph layout: a text split into graphemes (the characters a reader sees) and where each
// grapheme's pen starts, from a measure function for one font (canvas `measureText` in the apps), cached per font.
// The font atlas lays labels out with it (graphics `text-builder.ts`).
//
// Adapted from pretext (https://github.com/chenglou/pretext, commit 0d102a5cd933611e86a7e60bda7b4baa5995c6b0): its
// `pair-context` fit, where a grapheme's advance after another is the two measured together less the first, so the
// kerning between them is in the advance; widths cached per font, dropped whole past a bound (its
// MAX_CACHED_SEGMENTS); graphemes from `Intl.Segmenter`, ASCII a code unit each. Not taken: line breaking, whitespace
// modes, bidi, emoji width correction and the per-engine profiles.
//
// pretext is MIT licensed: Copyright (c) 2026 Pretext contributors. Permission is hereby granted, free of charge,
// to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright notice and this permission notice
// shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

export { TextMeasure, splitGraphemes };

/** The widths a `TextMeasure` keeps; past it they are all dropped and measured again as asked (residue ids, numbers). */
const MaxCachedWidths = 1 << 16;

let segmenter: Intl.Segmenter | null | undefined;

function isAscii(text: string) {
  for (let i = 0, il = text.length; i < il; ++i) {
    if (text.charCodeAt(i) > 0x7f) return false;
  }
  return true;
}

/**
 * The graphemes of `text`, in order: an emoji with its modifiers, a letter with its combining marks or a flag is one.
 * Code points where `Intl.Segmenter` is missing.
 */
function splitGraphemes(text: string): string[] {
  if (isAscii(text)) return text.split('');
  if (segmenter === undefined) {
    segmenter =
      typeof Intl !== 'undefined' && 'Segmenter' in Intl
        ? new Intl.Segmenter(undefined, { granularity: 'grapheme' })
        : null;
  }
  if (!segmenter) return Array.from(text);
  const graphemes: string[] = [];
  for (const s of segmenter.segment(text)) graphemes.push(s.segment);
  return graphemes;
}

namespace TextMeasure {
  export interface Layout {
    readonly graphemes: readonly string[];
    /** Where each grapheme's pen starts, the first at 0. */
    readonly offsets: Float64Array;
    /** Where the last grapheme ends: the text measured whole, up to rounding. */
    readonly width: number;
  }
}

/** Widths of text in one font, from `measure` (a width in any unit), cached. */
class TextMeasure {
  private readonly widths = new Map<string, number>();

  /** The width of `text` alone. */
  width(text: string): number {
    let w = this.widths.get(text);
    if (w === undefined) {
      if (this.widths.size >= MaxCachedWidths) this.widths.clear();
      w = this.measure(text);
      this.widths.set(text, w);
    }
    return w;
  }

  /** How far `grapheme` moves the pen after `previous`: the kerning between the two included (pretext's pair-context). */
  advance(grapheme: string, previous?: string): number {
    if (!previous) return this.width(grapheme);
    return this.width(previous + grapheme) - this.width(previous);
  }

  /**
   * The graphemes of `text` and where each starts. A grapheme starts where the pair it closes ends less its own width,
   * so the pair's kerning moves it, not the grapheme before it.
   */
  layout(text: string): TextMeasure.Layout {
    const graphemes = splitGraphemes(text);
    const offsets = new Float64Array(graphemes.length);
    let end = 0;
    for (let i = 0, il = graphemes.length; i < il; ++i) {
      const g = graphemes[i];
      end += this.advance(g, i > 0 ? graphemes[i - 1] : undefined);
      offsets[i] = end - this.width(g);
    }
    return { graphemes, offsets, width: end };
  }

  /** Drops every width, e.g. when the font behind `measure` changed (a web font arrived). */
  clear() {
    this.widths.clear();
  }

  constructor(private readonly measure: (text: string) => number) {}
}
