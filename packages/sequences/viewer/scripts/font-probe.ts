/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What `metricsMeasurer`'s tables hold, read from a browser (seqspec:viewer#the-render-backend,
// seqspec:viewer#text-measuring): per family and weight, the advance of each ASCII character (32..126) and of "é" (the
// `fallback`), in thousandths of an em, and the em box Canvas2D lays its `top`, `middle` and `bottom` baselines on: its
// ascent and descent above and below the alphabetic baseline. The box is read from where Canvas2D puts the flat foot of
// an "H" under each baseline, at 1,000 px, so it is whatever the browser's text stack does, not what a font file says
// (the font's own ascent differs: Chrome takes em-box edges from other metrics). Also the ink tops of "H" and "x" above
// the alphabetic baseline (`cap`, `xHeight`), which letters are centered on.
//
// `probeFonts` runs in the page: it is self-contained (the generator serializes it into Chrome; the browser test calls it).

/** The families the committed tables are made of, by the name `metricsMeasurer` finds them under, and their CSS. */
export const FONT_FAMILIES: Record<string, string> = {
  "ui-sans-serif": "ui-sans-serif, system-ui, sans-serif",
  "ui-monospace": "ui-monospace, Menlo, monospace",
  "sans-serif": "sans-serif",
  monospace: "monospace",
  serif: "serif",
};
export const FONT_WEIGHTS = [400, 600];

export interface ProbedFont {
  family: string;
  weight: number;
  advances: number[];
  fallback: number;
  ascent: number;
  descent: number;
  cap: number;
  xHeight: number;
}

export function probeFonts({
  families,
  weights,
}: {
  families: Record<string, string>;
  weights: number[];
}): ProbedFont[] {
  const SIZE = 1000;
  const WIDTHS_AT = 12;
  const out: ProbedFont[] = [];
  const c = new OffscreenCanvas(1400, 3200);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  /** The lowest row with ink in the "H" drawn with `baseline` at y = 1,600. */
  const foot = (baseline: CanvasTextBaseline): number => {
    g.clearRect(0, 0, 1400, 3200);
    g.textBaseline = baseline;
    g.fillStyle = "#000";
    g.fillText("H", 100, 1600);
    const d = g.getImageData(100, 0, 600, 3200).data;
    for (let y = 3199; y >= 0; y--) for (let x = 0; x < 600; x++) if (d[(y * 600 + x) * 4 + 3] > 127) return y + 1;
    return 0;
  };
  /** The highest row with ink in `ch` drawn on the alphabetic baseline at y = 1,600. */
  const head = (ch: string): number => {
    g.clearRect(0, 0, 1400, 3200);
    g.textBaseline = "alphabetic";
    g.fillStyle = "#000";
    g.fillText(ch, 100, 1600);
    const d = g.getImageData(100, 0, 600, 3200).data;
    for (let y = 0; y < 3200; y++) for (let x = 0; x < 600; x++) if (d[(y * 600 + x) * 4 + 3] > 127) return y;
    return 1600;
  };
  for (const [family, css] of Object.entries(families))
    for (const weight of weights) {
      g.textBaseline = "alphabetic";
      g.font = `${weight} ${WIDTHS_AT}px ${css}`;
      const em = (s: string) => Math.round((g.measureText(s).width / WIDTHS_AT) * 1000);
      const advances: number[] = [];
      for (let ch = 32; ch < 127; ch++) advances.push(em(String.fromCharCode(ch)));
      const fallback = em("é");
      g.font = `${weight} ${SIZE}px ${css}`;
      const alphabetic = foot("alphabetic");
      out.push({
        family,
        weight,
        advances,
        fallback,
        ascent: foot("top") - alphabetic,
        descent: alphabetic - foot("bottom"),
        cap: alphabetic - head("H"),
        xHeight: alphabetic - head("x"),
      });
    }
  return out;
}
