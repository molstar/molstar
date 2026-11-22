/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Core's SVG read back as what it drew where (the test's side of the field test, never the consumer's): letters,
// cells, range pills, labels, bars and heads, by the attributes the SVG backend writes. The consumer imports nothing
// of this file.

export interface Box {
  readonly x0: number;
  readonly x1: number;
  readonly fill: string;
  readonly y: number;
  readonly h: number;
  /** The fill's opacity (`fill-opacity`; 1 when absent). */
  readonly alpha: number;
}

export interface CoreDrawn {
  /** Each track's plot box, by its accessible name (`label (system)`). */
  readonly tracks: ReadonlyMap<string, { x: number; y: number; width: number; height: number }>;
  /** Rects (cells, bars) and rounded shapes (range pills), in draw order. */
  readonly rects: readonly Box[];
  readonly paths: readonly Box[];
  /** Letters: their x (cell centre) and character, and the fill the run is drawn in. */
  readonly letters: readonly { x: number; c: string; fill: string }[];
  /** Text drawn at 10px (labels inside marks) and 12px (track labels). */
  readonly texts: readonly { x: number; y: number; size: number; text: string }[];
  /** The accessible names of the items drawn (`label, system from–to`), the ranges' labels in full. */
  readonly items: readonly string[];
  /** Everything else that holds a fill, by tag: circles. */
  readonly circles: readonly { cx: number; cy: number; r: number; fill: string }[];
}

const attr = (tag: string, name: string) => new RegExp(`\\s${name}="([^"]*)"`).exec(tag)?.[1];

/** `#rrggbb` of a fill that may carry an alpha. */
const rgb = (c: string | undefined) => (c ?? '').slice(0, 7);

export function readCoreSvg(svg: string): CoreDrawn {
  const tracks = new Map<string, { x: number; y: number; width: number; height: number }>();
  for (const m of svg.matchAll(
    /<g data-seq-track="" role="group" aria-label="([^"]*)"><title>[^<]*<\/title><rect ([^>]*)\/>/g,
  )) {
    const r = m[2];
    tracks.set(m[1], {
      x: Number(attr(` ${r}`, 'x')),
      y: Number(attr(` ${r}`, 'y')),
      width: Number(attr(` ${r}`, 'width')),
      height: Number(attr(` ${r}`, 'height')),
    });
  }
  const data = svg.slice(svg.indexOf('data-pass="data"'), svg.indexOf('data-seq-accessible'));
  const rects: Box[] = [];
  for (const m of data.matchAll(/<rect ([^>]*)\/>/g)) {
    const t = ` ${m[1]}`;
    const x = Number(attr(t, 'x'));
    rects.push({
      x0: x,
      x1: x + Number(attr(t, 'width')),
      fill: rgb(attr(t, 'fill')),
      y: Number(attr(t, 'y')),
      h: Number(attr(t, 'height')),
      alpha: Number(attr(t, 'fill-opacity') ?? 1),
    });
  }
  const paths: Box[] = [];
  for (const m of data.matchAll(/<path d="([^"]*)" fill="([^"]*)"(?: fill-opacity="([^"]*)")?\/>/g)) {
    const nums = m[1].match(/-?\d+(?:\.\d+)?/g)?.map(Number) ?? [];
    // an even index is x and an odd one y: the commands carry coordinates only
    const xs = nums.filter((_, i) => i % 2 === 0);
    const ys = nums.filter((_, i) => i % 2 === 1);
    paths.push({
      x0: Math.min(...xs),
      x1: Math.max(...xs),
      fill: rgb(m[2]),
      y: Math.min(...ys),
      h: Math.max(...ys) - Math.min(...ys),
      alpha: Number(m[3] ?? 1),
    });
  }
  const letters: { x: number; c: string; fill: string }[] = [];
  const texts: { x: number; y: number; size: number; text: string }[] = [];
  for (const m of svg.matchAll(/<text ([^>]*)>(.*?)<\/text>/g)) {
    const t = ` ${m[1]}`;
    if (m[2].startsWith('<tspan')) {
      // a run of letters: a tspan each, at its cell's centre, in its ink
      for (const s of m[2].matchAll(/<tspan x="([^"]*)" fill="([^"]*)">(.)<\/tspan>/g))
        letters.push({ x: Number(s[1]), c: s[3], fill: rgb(s[2]) });
    } else
      texts.push({ x: Number(attr(t, 'x')), y: Number(attr(t, 'y')), size: Number(attr(t, 'font-size')), text: m[2] });
  }
  const items = [...svg.matchAll(/<g data-seq-item="" role="img" aria-label="([^"]*)">/g)].map((m) => m[1]);
  const circles: { cx: number; cy: number; r: number; fill: string }[] = [];
  for (const m of data.matchAll(/<circle ([^>]*)\/>/g)) {
    const t = ` ${m[1]}`;
    circles.push({
      cx: Number(attr(t, 'cx')),
      cy: Number(attr(t, 'cy')),
      r: Number(attr(t, 'r')),
      fill: rgb(attr(t, 'fill')),
    });
  }
  return { tracks, rects, paths, letters, texts, items, circles };
}
