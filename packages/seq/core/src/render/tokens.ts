/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashString } from '@molstar/core/data/util/hash-functions';

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
  'seq-pending',
  'seq-conflict',
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
  /** A hash of the colors and fonts: equal tables have equal versions. */
  readonly version: number;
}

/** `@molstar/ui`'s font tokens (`--msui-font-sans`, `--msui-font-mono`). */
export const DefaultFonts: FontStacks = {
  sans: '"Helvetica Neue", "Segoe UI", Helvetica, "Source Sans Pro", Arial, sans-serif',
  mono: 'ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace',
};

/** A token table from colors as CSS hex strings; a token left out, or not a hex color, is an error. */
export function createTokens(theme: string, colors: Record<SeqTokenName, string>, fonts = DefaultFonts): SeqTokens {
  const out = {} as Record<SeqTokenName, number>;
  let key = `${theme}|${fonts.sans}|${fonts.mono}`;
  for (const name of SeqTokenNames) {
    const c = parseColor(colors[name]);
    if (c === undefined) throw new Error(`The token ${name} is not a hex color: ${colors[name]}`);
    out[name] = c;
    key += `|${c}`;
  }
  return { theme, colors: out, fonts, version: hashString(key) >>> 0 };
}

// The built-in tables follow `@molstar/ui`'s light and dark themes (`--msui-bg`, `--msui-fg`, `--msui-accent`); the
// derived chrome tokens are their mixes, computed once here. Residue classes and the categorical palette are pastel
// fills that keep dark ink at 4.5:1 or better in the light theme, and deep fills under light ink in the dark one.

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
  'seq-pending': '#dedad3',
  'seq-conflict': '#c2410c',
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
  'seq-cat-0': '#7a8fa6',
  'seq-cat-1': '#a6857a',
  'seq-cat-2': '#7aa68a',
  'seq-cat-3': '#9a7aa6',
  'seq-cat-4': '#a69c7a',
  'seq-cat-5': '#7aa3a6',
  'seq-cat-6': '#a67a95',
  'seq-cat-7': '#8c8c8c',
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
  'seq-pending': '#262a31',
  'seq-conflict': '#f97316',
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
  'seq-cat-0': '#6f86a3',
  'seq-cat-1': '#a3806f',
  'seq-cat-2': '#6fa384',
  'seq-cat-3': '#946fa3',
  'seq-cat-4': '#a3986f',
  'seq-cat-5': '#6f9fa3',
  'seq-cat-6': '#a36f8f',
  'seq-cat-7': '#8a8f98',
});

/** The built-in table of a theme (Node, exports): `light` unless `dark`. */
export function builtInTokens(theme: string): SeqTokens {
  return theme === 'dark' ? DarkTokens : LightTokens;
}

/** `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa` as `0xRRGGBBAA`; undefined for anything else. */
export function parseColor(s: string | undefined): number | undefined {
  if (typeof s !== 'string' || !/^#([0-9a-f]{3,4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(s)) return undefined;
  let h = s.slice(1);
  if (h.length <= 4) h = [...h].map((c) => c + c).join('');
  if (h.length === 6) h += 'ff';
  return Number.parseInt(h, 16) >>> 0;
}

/**
 * A spec's color: a literal (`#rrggbb`) or a token name (`token:accent` reads `seq-accent`, then `msui-accent`, then
 * `accent`). Undefined when it is neither: the caller reports it and uses its default.
 */
export function resolveColor(value: unknown, tokens: SeqTokens): number | undefined {
  if (typeof value !== 'string') return undefined;
  if (value.startsWith('token:')) {
    const name = value.slice(6);
    const colors = tokens.colors as Readonly<Record<string, number>>;
    return colors[`seq-${name}`] ?? colors[`msui-${name}`] ?? colors[name];
  }
  return parseColor(value);
}

/** `#rrggbbaa`. */
export function colorHex(c: number): string {
  return `#${(c >>> 0).toString(16).padStart(8, '0')}`;
}

/** The color with its alpha replaced (0 to 1). */
export function withAlpha(c: number, alpha: number): number {
  return ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, alpha)) * 255)) >>> 0;
}

/** Relative luminance (WCAG) of a color's RGB, 0 to 1. */
export function luminance(c: number): number {
  const ch = (v: number) => {
    const s = v / 255;
    return s <= 0.03928 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  };
  return 0.2126 * ch((c >>> 24) & 0xff) + 0.7152 * ch((c >>> 16) & 0xff) + 0.0722 * ch((c >>> 8) & 0xff);
}

/** Ink follows its fill: the table's dark or light ink, whichever contrasts more with `fill`. */
export function inkOn(fill: number, tokens: SeqTokens): number {
  const dark = tokens.colors['seq-ink-dark'];
  const light = tokens.colors['seq-ink-light'];
  const l = luminance(fill);
  const contrast = (a: number, b: number) => (Math.max(a, b) + 0.05) / (Math.min(a, b) + 0.05);
  return contrast(l, luminance(dark)) >= contrast(l, luminance(light)) ? dark : light;
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

const AminoClass: Readonly<Record<string, SeqTokenName>> = {};
for (const [letters, name] of [
  ['AVLIMC', 'seq-aa-hydrophobic'],
  ['FWY', 'seq-aa-aromatic'],
  ['STNQ', 'seq-aa-polar'],
  ['KRH', 'seq-aa-positive'],
  ['DE', 'seq-aa-negative'],
  ['GP', 'seq-aa-special'],
] as const)
  for (const l of letters) (AminoClass as Record<string, SeqTokenName>)[l] = name;

const NucleotideClass: Readonly<Record<string, SeqTokenName>> = {
  A: 'seq-nt-a',
  C: 'seq-nt-c',
  G: 'seq-nt-g',
  T: 'seq-nt-t',
  U: 'seq-nt-t',
};

/** A letter's cell fill by its class in the alphabet; missing for a gap or an unknown letter. */
export function letterToken(letter: string, alphabet: 'protein' | 'dna' | 'rna'): SeqTokenName {
  const u = letter.toUpperCase();
  if (u === '-' || u === '.') return 'seq-missing';
  if (alphabet === 'protein') return AminoClass[u] ?? 'seq-missing';
  return NucleotideClass[u] ?? 'seq-nt-other';
}

/** A category's palette entry: stable for a name (its hash), so a category keeps its color across tiles. */
export function categoryToken(category: string): SeqTokenName {
  return `seq-cat-${(hashString(category) >>> 0) % CategoryCount}` as SeqTokenName;
}
