/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type CanvasInk, parseColor } from '@molstar/prism-ui/core';
import { DarkSkin, DefaultSkin, type Skin } from '@molstar/seq-viewer';

// The sequence view's colors from the canvas ink (prism-sequence.md §5): Seq*'s skin, written from the ink's theme and
// background and Prism's palette for that theme (`theme.css`: the text and accent colors of the light and dark
// themes). Nothing reads the DOM or CSS: the same ink gives the same skin in a story, a shadow root and a test.

/** Prism's text and accent colors per canvas theme (styles/theme.css). */
const Palette = {
  light: { text: 0x332b1f, accent: 0xae5d04 },
  dark: { text: 0xccd4e0, accent: 0x51a2fb },
} as const;

const channel = (c: number, shift: number) => (c >> shift) & 0xff;
const rgb = (c: number, alpha?: number) =>
  `rgb(${channel(c, 16)} ${channel(c, 8)} ${channel(c, 0)}${alpha === undefined ? '' : ` / ${Math.round(alpha * 100) / 100}`})`;

/** `from` stepped `t` (0..1) toward `to`, per channel (sRGB). */
function mix(from: number, to: number, t: number): number {
  const at = (shift: number) => Math.round(channel(from, shift) + (channel(to, shift) - channel(from, shift)) * t);
  return (at(16) << 16) | (at(8) << 8) | at(0);
}

/** The fonts the view draws with: the host's (`--msp-font-mono`, `--msp-font-sans`), not part of the ink. */
export interface SequenceFonts {
  readonly mono: string;
  readonly ui: string;
  /** The panel text size, px (`--msp-font-size-sm`): the text version's letters are at most this large. */
  readonly size?: number;
}

/**
 * The skin for an ink: the surface is the canvas's color (the view draws on it), the ink and accent are Prism's for the
 * canvas theme, the quieter inks and the grid are steps from the surface toward them, the selection is the accent at
 * 25% (as the mock strip draws it). `fonts` replaces Seq*'s default stacks.
 */
export function skinFor(ink: CanvasInk, fonts?: SequenceFonts): Skin {
  const dark = ink.theme === 'dark';
  const base = dark ? DarkSkin : DefaultSkin;
  const { text, accent } = Palette[ink.theme];
  const surface = parseColor(ink.color) ?? (dark ? 0x111318 : 0xffffff);
  // Shades go toward the text color (white on dark, black on light), as `--msp-shade` does.
  const toward = (t: number) => mix(surface, text, t);
  return {
    ...base,
    dark,
    surface: rgb(surface),
    surfaceAlt: rgb(toward(0.05)),
    ink: rgb(text),
    inkQuiet: rgb(mix(text, surface, 0.33)),
    inkFaint: rgb(mix(text, surface, 0.7)),
    grid: rgb(toward(0.12)),
    accent: rgb(accent),
    hover: rgb(accent, 0.2),
    selection: rgb(accent, 0.25),
    selectionBorder: rgb(accent, 0.7),
    // Positions not modelled in the structure: a quiet text tone, readable as letters (Seq*'s `unobservedGray`).
    unobserved: rgb(mix(text, surface, 0.55)),
    pending: rgb(toward(0.1)),
    grip: rgb(mix(text, surface, 0.45)),
    gripHover: rgb(text),
    window: rgb(text, 0.1),
    focusRing: rgb(accent),
    ...(fonts ? { fontMono: fonts.mono, fontUi: fonts.ui } : {}),
  };
}
