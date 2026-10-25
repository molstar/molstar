/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dialects' color names as RGB (`0xRRGGBB`): PyMOL's named colors (the common ones of its color table, as
// fractions there), `greyNN`/`grayNN`, and VMD's ColorIDs 0–32 with their default colors. One table serves the
// commands, the VMD look's preset and completion, so they cannot disagree (_spec §5.3).

const rgb = (r: number, g: number, b: number) =>
  (Math.round(r * 255) << 16) | (Math.round(g * 255) << 8) | Math.round(b * 255);

/** PyMOL's named colors, as its color table gives them (fractions). */
const PymolTable: Readonly<Record<string, readonly [number, number, number]>> = {
  white: [1, 1, 1],
  black: [0, 0, 0],
  blue: [0, 0, 1],
  green: [0, 1, 0],
  red: [1, 0, 0],
  cyan: [0, 1, 1],
  yellow: [1, 1, 0],
  magenta: [1, 0, 1],
  orange: [1, 0.5, 0],
  pink: [1, 0.65, 0.85],
  grey: [0.5, 0.5, 0.5],
  gray: [0.5, 0.5, 0.5],
  wheat: [0.99, 0.82, 0.65],
  salmon: [1, 0.6, 0.6],
  slate: [0.5, 0.5, 1],
  lightblue: [0.75, 0.75, 1],
  marine: [0, 0.5, 1],
  purple: [0.75, 0, 0.75],
  violet: [1, 0.5, 1],
  lime: [0.5, 1, 0.5],
  forest: [0.2, 0.6, 0.2],
  teal: [0, 0.75, 0.75],
  deepteal: [0.1, 0.6, 0.6],
  hotpink: [1, 0, 0.5],
  skyblue: [0.2, 0.5, 0.8],
  palegreen: [0.65, 0.9, 0.65],
  lightpink: [1, 0.75, 0.87],
  olive: [0.77, 0.7, 0],
  brown: [0.65, 0.32, 0.17],
  firebrick: [0.698, 0.13, 0.13],
  raspberry: [0.7, 0.3, 0.4],
  ruby: [0.6, 0.2, 0.2],
  chocolate: [0.555, 0.222, 0.111],
  sand: [0.72, 0.55, 0.3],
  tv_red: [1, 0.2, 0.2],
  tv_green: [0.2, 1, 0.2],
  tv_blue: [0.3, 0.3, 1],
  tv_yellow: [1, 1, 0.2],
  tv_orange: [1, 0.55, 0.15],
  lightorange: [1, 0.8, 0.5],
  yelloworange: [1, 0.87, 0.37],
  limon: [0.75, 1, 0.25],
  splitpea: [0.52, 0.75, 0],
  chartreuse: [0.5, 1, 0],
  smudge: [0.55, 0.7, 0.4],
  deepolive: [0.6, 0.6, 0.1],
  limegreen: [0, 1, 0.5],
  density: [0.1, 0.1, 0.6],
  deepblue: [0.25, 0.25, 0.65],
  bluewhite: [0.85, 0.85, 1],
  lightteal: [0.4, 0.7, 0.7],
  aquamarine: [0.5, 1, 1],
  palecyan: [0.8, 1, 1],
  greencyan: [0.25, 1, 0.75],
  deepsalmon: [1, 0.5, 0.5],
  warmpink: [0.85, 0.2, 0.5],
  dirtyviolet: [0.7, 0.5, 0.5],
  lightmagenta: [1, 0.2, 0.8],
  purpleblue: [0.5, 0.2, 1],
  deeppurple: [0.6, 0.1, 0.6],
  violetpurple: [0.55, 0.25, 0.6],
  paleyellow: [1, 1, 0.5],
  carbon: [0.2, 1, 0.2],
  nitrogen: [0.2, 0.2, 1],
  oxygen: [1, 0.3, 0.3],
  hydrogen: [0.9, 0.9, 0.9],
  sulfur: [0.9, 0.775, 0.25],
  silver: [0.75, 0.75, 0.75],
  gold: [1, 0.82, 0.14],
};

/** A PyMOL color: a name of its table, `greyNN`/`grayNN`, `0xRRGGBB`, `[r, g, b]` (0–1 or 0–255). */
export function pymolColor(text: string): number | undefined {
  const t = text.trim().toLowerCase();
  const named = PymolTable[t];
  if (named) return rgb(...named);
  const grey = t.match(/^gr[ae]y(\d\d)$/);
  if (grey) {
    const v = Number(grey[1]) / 99;
    return rgb(v, v, v);
  }
  const hex = t.match(/^(?:0x|#)([0-9a-f]{6})$/);
  if (hex) return parseInt(hex[1], 16);
  const list = t.match(/^\[\s*([\d.]+)\s*,\s*([\d.]+)\s*,\s*([\d.]+)\s*\]$/);
  if (list) {
    const xs = [Number(list[1]), Number(list[2]), Number(list[3])];
    if (xs.some((x) => !Number.isFinite(x) || x < 0)) return undefined;
    const scale = xs.some((x) => x > 1) ? 255 : 1;
    if (xs.some((x) => x > scale)) return undefined;
    return rgb(xs[0] / scale, xs[1] / scale, xs[2] / scale);
  }
  return undefined;
}

/** PyMOL's color names (for completion and "did you mean"). */
export const PymolColorNames: readonly string[] = Object.keys(PymolTable);

/** VMD's ColorIDs by name, 0–32 in order. */
export const VmdColorIds: Readonly<Record<string, number>> = {
  blue: 0,
  red: 1,
  gray: 2,
  orange: 3,
  yellow: 4,
  tan: 5,
  silver: 6,
  green: 7,
  white: 8,
  pink: 9,
  cyan: 10,
  purple: 11,
  lime: 12,
  mauve: 13,
  ochre: 14,
  iceblue: 15,
  black: 16,
  yellow2: 17,
  yellow3: 18,
  green2: 19,
  green3: 20,
  cyan2: 21,
  cyan3: 22,
  blue2: 23,
  blue3: 24,
  violet: 25,
  violet2: 26,
  magenta: 27,
  magenta2: 28,
  red2: 29,
  red3: 30,
  orange2: 31,
  orange3: 32,
};

/** VMD's default colors of its ColorIDs, in id order. */
const VmdTable: readonly (readonly [number, number, number])[] = [
  [0, 0, 1],
  [1, 0, 0],
  [0.35, 0.35, 0.35],
  [1, 0.5, 0],
  [1, 1, 0],
  [0.5, 0.5, 0.2],
  [0.6, 0.6, 0.6],
  [0, 1, 0],
  [1, 1, 1],
  [1, 0.6, 0.6],
  [0.25, 0.75, 0.75],
  [0.65, 0, 0.65],
  [0.5, 0.9, 0.4],
  [0.9, 0.4, 0.7],
  [0.5, 0.3, 0],
  [0.5, 0.5, 0.75],
  [0, 0, 0],
  [0.88, 0.97, 0.02],
  [0.55, 0.9, 0.02],
  [0, 0.9, 0.04],
  [0, 0.9, 0.5],
  [0, 0.88, 1],
  [0, 0.76, 1],
  [0.02, 0.38, 0.67],
  [0.01, 0.04, 0.93],
  [0.27, 0, 0.98],
  [0.45, 0, 0.9],
  [0.9, 0, 0.9],
  [1, 0, 0.66],
  [0.98, 0, 0.23],
  [0.81, 0, 0],
  [0.89, 0.35, 0],
  [0.96, 0.72, 0],
];

/** A VMD color by ColorID (0–32) or by its name. */
export function vmdColor(idOrName: string): number | undefined {
  const t = idOrName.trim().toLowerCase();
  const id = /^\d+$/.test(t) ? Number(t) : VmdColorIds[t];
  const c = id === undefined ? undefined : VmdTable[id];
  return c ? rgb(...c) : undefined;
}

/**
 * VMD's Name coloring (its default): the color of an element's first letter, by ColorID. Carbon is cyan, nitrogen
 * blue, oxygen red, sulfur yellow, hydrogen white, phosphorus tan; others are VMD's "other" (tan here).
 */
export const VmdNameColors: Readonly<Record<string, number>> = {
  C: vmdColor('cyan')!,
  N: vmdColor('blue')!,
  O: vmdColor('red')!,
  S: vmdColor('yellow')!,
  H: vmdColor('white')!,
  P: vmdColor('tan')!,
};
