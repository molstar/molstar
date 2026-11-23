/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';
import { contrastRgba, inkFor, parseColor } from '@molstar/seqx-core/theme/index';

// Colors without CSS (_spec/README.md "Render", Tokens): the token table resolves every data token (`seq-*`) and the
// chrome tokens (`msui-*`) the draw lists use to 32-bit RGBA (`0xRRGGBBAA`). In the browser the viewer reads it once
// from the computed style of its root; in Node it is a built-in table per theme. Its version is a hash of its
// content, part of the draw lists' key: a theme switch redraws, an equal table does not.

/** The tokens the draw lists use. */
export const SeqTokenNames = [
  'msui-bg',
  'msui-fg',
  'msui-text',
  'msui-text-muted',
  'msui-line',
  'msui-accent',
  'msui-accent-muted',
  'seq-track',
  'seq-sequence',
  'seq-ink-dark',
  'seq-ink-light',
  'seq-missing',
  'seq-genotype-0',
  'seq-genotype-1',
  'seq-genotype-2',
  'seq-genotype-3',
  'seq-pending',
  'seq-conflict',
  'seq-modified',
  'seq-break',
  'seq-range',
  'seq-value',
  'seq-hover',
  'seq-density-low',
  'seq-density-high',
  'seq-aa-hydrophobic',
  'seq-aa-aromatic',
  'seq-aa-polar',
  'seq-aa-positive',
  'seq-aa-negative',
  'seq-aa-special',
  'seq-nt-a',
  'seq-nt-c',
  'seq-nt-g',
  'seq-nt-t',
  'seq-nt-other',
  'seq-cat-0',
  'seq-cat-1',
  'seq-cat-2',
  'seq-cat-3',
  'seq-cat-4',
  'seq-cat-5',
  'seq-cat-6',
  'seq-cat-7',
  'seq-ramp-0',
  'seq-ramp-1',
  'seq-ramp-2',
  'seq-ramp-3',
  'seq-ramp-4',
  'seq-div-neg',
  'seq-div-mid',
  'seq-div-pos',
  'seq-ss-helix',
  'seq-ss-strand',
  'seq-ss-turn',
  'seq-variant',
  'seq-ruler',
  'seq-window',
  'seq-grip',
  'seq-wash',
  'seq-centromere',
  'seq-cell-edge',
] as const;

export type SeqTokenName = (typeof SeqTokenNames)[number];

/** The categorical palette's size. */
export const CategoryCount = 8;

/** Font stacks by family role: a draw list's font names a role (`12px sans`), a backend writes its stack. */
export interface FontStacks {
  readonly sans: string;
  readonly mono: string;
}

export interface SeqTokens {
  readonly theme: string;
  readonly colors: Readonly<Record<SeqTokenName, number>>;
  readonly fonts: FontStacks;
  /** `seq-logo-chroma`: the floor of a logo glyph's Oklab chroma (a number, not a color). */
  readonly logoChroma: number;
  /** A hash of the colors and fonts: equal tables have equal versions. */
  readonly version: number;
}

/** `@molstar/ui`'s font tokens (`--msui-font-sans`, `--msui-font-mono`). */
export const DefaultFonts: FontStacks = {
  sans: '"Helvetica Neue", "Segoe UI", Helvetica, "Source Sans Pro", Arial, sans-serif',
  mono: 'ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace',
};

/** `seq-logo-chroma` per theme: the logo glyphs' Oklab chroma floor (tuned on a globin family's logo in each theme). */
export const LightLogoChroma = 0.12;
export const DarkLogoChroma = 0.1;

/** A token table from colors as CSS hex strings; a token left out, or not a hex color, is an error. */
export function createTokens(
  theme: string,
  colors: Record<SeqTokenName, string>,
  fonts = DefaultFonts,
  logoChroma = theme === 'dark' ? DarkLogoChroma : LightLogoChroma,
): SeqTokens {
  const out = {} as Record<SeqTokenName, number>;
  let key = `${theme}|${fonts.sans}|${fonts.mono}|${logoChroma}`;
  for (const name of SeqTokenNames) {
    const c = parseColor(colors[name]);
    if (c === undefined) throw new Error(`The token ${name} is not a hex color: ${colors[name]}`);
    out[name] = c;
    key += `|${c}`;
  }
  return { theme, colors: out, fonts, logoChroma, version: hashString(key) >>> 0 };
}

// The built-in tables follow `@molstar/ui`'s light and dark themes (`--msui-bg`, `--msui-fg`, `--msui-accent`); the
// derived chrome tokens are their mixes, computed once here. Residue classes are pastel fills that keep dark ink at
// 4.5:1 or better in the light theme, and deep fills under light ink in the dark one. The categories are Okabe–Ito's
// seven hues (mixed a tenth toward the surface) and a red, at least 0.06 apart in Oklab under normal vision and the
// three simulated deficiencies in both tables (_spec/README.md "Render", Tokens); the ramp is monotone in lightness
// from near the track's color to strong; `seq-centromere` holds 3:1 against the surface. `seq-cell-edge` is a letter
// cell's hairline in dark (3:1 against the page: the dark cells sit at 1.4–1.7:1 and cannot go lighter without their
// ink failing 4.5:1), none (transparent) in light.

export const LightTokens: SeqTokens = createTokens('light', {
  'msui-bg': '#f1efea',
  'msui-fg': '#2d271e',
  'msui-text': '#2d271e',
  'msui-text-muted': '#6b655c',
  'msui-line': '#d3d0c9',
  'msui-accent': '#8f4c06',
  'msui-accent-muted': '#e4d8c9',
  'seq-track': '#e9e6e0',
  'seq-sequence': '#c4c9cf',
  'seq-ink-dark': '#2d271e',
  'seq-ink-light': '#ffffff',
  'seq-missing': '#c9c5bd',
  'seq-genotype-0': '#dfe7f0',
  'seq-genotype-1': '#7c9dc2',
  'seq-genotype-2': '#2c4f78',
  'seq-genotype-3': '#aec3da',
  'seq-pending': '#dedad3',
  'seq-conflict': '#c2410c',
  'seq-modified': '#a21caf',
  'seq-break': '#8a847a',
  'seq-range': '#7a8fa6',
  'seq-value': '#5b7c99',
  'seq-hover': '#2d271e',
  'seq-density-low': '#e3e8ee',
  'seq-density-high': '#36506b',
  'seq-aa-hydrophobic': '#cfe0f5',
  'seq-aa-aromatic': '#d8d2f2',
  'seq-aa-polar': '#cdebd6',
  'seq-aa-positive': '#f5cfcf',
  'seq-aa-negative': '#f6d9be',
  'seq-aa-special': '#efe8b8',
  'seq-nt-a': '#cdebd6',
  'seq-nt-c': '#cfe0f5',
  'seq-nt-g': '#f6e2b6',
  'seq-nt-t': '#f5cfcf',
  'seq-nt-other': '#dedad3',
  'seq-cat-0': '#187eb8',
  'seq-cat-1': '#e7a717',
  'seq-cat-2': '#18a67f',
  'seq-cat-3': '#d085ae',
  'seq-cat-4': '#66bae9',
  'seq-cat-5': '#d86c17',
  'seq-cat-6': '#f0e553',
  'seq-cat-7': '#d23c3c',
  'seq-ramp-0': '#dfe7f0',
  'seq-ramp-1': '#aec3da',
  'seq-ramp-2': '#7c9dc2',
  'seq-ramp-3': '#4f76a3',
  'seq-ramp-4': '#2c4f78',
  'seq-div-neg': '#2f6fb0',
  'seq-div-mid': '#e6e2dc',
  'seq-div-pos': '#b8452a',
  'seq-ss-helix': '#d9536f',
  'seq-ss-strand': '#e8a33d',
  'seq-ss-turn': '#5aa0c8',
  'seq-variant': '#7b4fa6',
  'seq-ruler': '#8a847a',
  'seq-window': '#8f4c06',
  'seq-grip': '#6b655c',
  'seq-wash': '#f1efea33',
  'seq-centromere': '#b4306a',
  'seq-cell-edge': '#00000000',
});

export const DarkTokens: SeqTokens = createTokens('dark', {
  'msui-bg': '#121419',
  'msui-fg': '#cfd6e1',
  'msui-text': '#cfd6e1',
  'msui-text-muted': '#8d95a1',
  'msui-line': '#2f333b',
  'msui-accent': '#5aa6f8',
  'msui-accent-muted': '#1d2f45',
  'seq-track': '#1a1d23',
  'seq-sequence': '#3b424c',
  'seq-ink-dark': '#121419',
  'seq-ink-light': '#e8edf4',
  'seq-missing': '#3a3f48',
  'seq-genotype-0': '#1f2a36',
  'seq-genotype-1': '#44709e',
  'seq-genotype-2': '#a9c8ea',
  'seq-genotype-3': '#2f4d6e',
  'seq-pending': '#262a31',
  'seq-conflict': '#f97316',
  'seq-modified': '#e879f9',
  'seq-break': '#7d8592',
  'seq-range': '#6f86a3',
  'seq-value': '#7ea6c9',
  'seq-hover': '#e8edf4',
  'seq-density-low': '#1f2731',
  'seq-density-high': '#9cc3e6',
  'seq-aa-hydrophobic': '#28456a',
  'seq-aa-aromatic': '#3d3466',
  'seq-aa-polar': '#24533a',
  'seq-aa-positive': '#6a2a2a',
  'seq-aa-negative': '#6b4220',
  'seq-aa-special': '#5c5520',
  'seq-nt-a': '#24533a',
  'seq-nt-c': '#28456a',
  'seq-nt-g': '#6b5520',
  'seq-nt-t': '#6a2a2a',
  'seq-nt-other': '#30343c',
  'seq-cat-0': '#0269a3',
  'seq-cat-1': '#d19102',
  'seq-cat-2': '#02906a',
  'seq-cat-3': '#b96f99',
  'seq-cat-4': '#4fa4d4',
  'seq-cat-5': '#c25702',
  'seq-cat-6': '#dacf3e',
  'seq-cat-7': '#ff3c69',
  'seq-ramp-0': '#1f2a36',
  'seq-ramp-1': '#2f4d6e',
  'seq-ramp-2': '#44709e',
  'seq-ramp-3': '#6f9bcb',
  'seq-ramp-4': '#a9c8ea',
  'seq-div-neg': '#5a9ad8',
  'seq-div-mid': '#2a2e36',
  'seq-div-pos': '#e0734f',
  'seq-ss-helix': '#e0708a',
  'seq-ss-strand': '#e9b25a',
  'seq-ss-turn': '#6fb2d8',
  'seq-variant': '#b08ad6',
  'seq-ruler': '#7d8592',
  'seq-window': '#5aa6f8',
  'seq-grip': '#8d95a1',
  'seq-wash': '#12141933',
  'seq-centromere': '#e86aa0',
  'seq-cell-edge': '#646b76',
});

/** The built-in table of a theme (Node, exports): `light` unless `dark`. */
export function builtInTokens(theme: string): SeqTokens {
  return theme === 'dark' ? DarkTokens : LightTokens;
}

// A spec's colors resolve through the theme entry (`resolveColor`: literals, `token:<name>`, `light-dark(a, b)`).
export { colorHex, parseColor, resolveColor } from '@molstar/seqx-core/theme/index';

/** The color with its alpha replaced (0 to 1). */
export function withAlpha(c: number, alpha: number): number {
  return ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, alpha)) * 255)) >>> 0;
}

/**
 * Ink follows its fill (4.5:1): the table's dark or light ink, whichever contrasts more with `fill` (Mol*'s
 * `Color.contrast`); pure black or white where neither reaches 4.5:1.
 */
export function inkOn(fill: number, tokens: SeqTokens): number {
  return inkFor(fill, tokens.colors['seq-ink-dark'], tokens.colors['seq-ink-light']);
}

/** Linear interpolation of two colors' channels, alpha included; `t` in 0 to 1. */
export function mixColors(a: number, b: number, t: number): number {
  let out = 0;
  for (let shift = 24; shift >= 0; shift -= 8) {
    const x = (a >>> shift) & 0xff;
    const y = (b >>> shift) & 0xff;
    out = out * 256 + Math.round(x + (y - x) * t);
  }
  return out >>> 0;
}

/** Muted text holds 4.5:1 against its well (_spec/README.md "Themes", Colors). */
export const MutedTextContrast = 4.5;

/**
 * The ink of a muted text (ruler tick labels, a ruler's note) on its wells: `muted` where it holds 4.5:1 on each (or is
 * the text color itself: forced colors map the muted ink to the system's text color, which is kept, never blended),
 * else the text color at 75% over the first well, opaque (`mixColors` with that well, so the ratio is the one drawn), if
 * that holds on every well, else the text color itself.
 */
export function mutedInk(tokens: SeqTokens, muted: number, wells: readonly number[]): number {
  const text = tokens.colors['msui-text'];
  const holds = (ink: number) => wells.every((w) => contrastRgba(ink, w) >= MutedTextContrast);
  if (muted === text || holds(muted)) return muted;
  const soft = (mixColors(wells[0], text, 0.75) | 0xff) >>> 0;
  return holds(soft) ? soft : text;
}
