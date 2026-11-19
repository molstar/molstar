/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  builtInTokens,
  colorHex,
  createTokens,
  type FontStacks,
  parseColor,
  type SeqTokenName,
  SeqTokenNames,
  type SeqTokens,
} from '@molstar/seqx-core/render';

// The browser's token table (core _spec "Render", Tokens; viewer _spec "Canvas2D backend", Tokens): the custom
// properties `--msui-*` (`@molstar/ui`'s theme) and `--seq-*` (a host's data colors) of a root's computed style, as
// the core's `SeqTokens`. A token the CSS does not define, or whose value is no color, takes the built-in table's of
// the theme the root shows.

export const ThemeAttribute = 'data-msui-theme';

/** The theme a root shows: `data-msui-theme` on its nearest ancestor (itself included) or the document element, else `light`. */
export function themeOf(root: Element): string {
  const own = root.closest?.(`[${ThemeAttribute}]`)?.getAttribute(ThemeAttribute);
  return own || root.ownerDocument?.documentElement.getAttribute(ThemeAttribute) || 'light';
}

/** The color of a CSS value as a hex string: hex, `rgb()`/`rgba()` and `color(srgb …)` (what a browser resolves to). */
export function cssColorToHex(value: string): string | undefined {
  const v = value.trim();
  if (parseColor(v) !== undefined) return v;
  const m = /^(rgba?|color)\(\s*(?:srgb\s+)?([^)]*)\)$/i.exec(v);
  if (!m) return undefined;
  const parts = m[2]!.split(/[\s,/]+/).filter(Boolean);
  if (parts.length < 3) return undefined;
  const unit = m[1]!.toLowerCase() === 'color';
  const ch = (s: string) => {
    const pct = s.endsWith('%');
    const n = Number.parseFloat(s);
    if (!Number.isFinite(n)) return undefined;
    const x = unit ? (pct ? n / 100 : n) * 255 : pct ? (n / 100) * 255 : n;
    return Math.max(0, Math.min(255, Math.round(x)));
  };
  const rgb = parts.slice(0, 3).map(ch);
  if (rgb.some((c) => c === undefined)) return undefined;
  let a = 1;
  if (parts[3] !== undefined) {
    const n = Number.parseFloat(parts[3]);
    if (!Number.isFinite(n)) return undefined;
    a = parts[3].endsWith('%') ? n / 100 : n;
  }
  const hex = (c: number) => c.toString(16).padStart(2, '0');
  return `#${rgb.map((c) => hex(c!)).join('')}${hex(Math.round(Math.max(0, Math.min(1, a)) * 255))}`;
}

/** The system colors forced colors name (`CanvasText` and the rest), as hex strings. */
export interface SystemColors {
  readonly Canvas: string;
  readonly CanvasText: string;
  readonly Highlight: string;
  readonly GrayText: string;
}

/** What a browser gives the system colors without a theme to say (its default light scheme), for one that cannot resolve them. */
export const DefaultSystemColors: SystemColors = {
  Canvas: '#ffffff',
  CanvasText: '#000000',
  Highlight: '#3399ff',
  GrayText: '#6d6d6d',
};

/** A color with an opacity (`#rrggbb` with `alpha` in 0..1). */
const withAlpha = (hex: string, alpha: number) =>
  `${hex.slice(0, 7)}${Math.round(alpha * 255)
    .toString(16)
    .padStart(2, '0')}`;

/**
 * The tokens under forced colors (viewer _spec "Canvas2D backend", Forced colors): the chrome and the emphasis (the
 * page, ink, lines, accent, selection and hover, the window and its grips, the cell hairline) take the system colors, so
 * the data's ink, the selection and the hover stay visible on whatever the reader's palette is. The data colors stay
 * (the data's meaning is in them), as `forced-color-adjust` leaves a legend's swatches.
 */
export function forcedColorTokens(tokens: SeqTokens, system: SystemColors = DefaultSystemColors): SeqTokens {
  const hex = (name: SeqTokenName) => colorHex(tokens.colors[name]);
  const colors = Object.fromEntries(SeqTokenNames.map((n) => [n, hex(n)])) as Record<SeqTokenName, string>;
  const { Canvas, CanvasText, Highlight, GrayText } = system;
  Object.assign(colors, {
    'msui-bg': Canvas,
    'msui-fg': CanvasText,
    'msui-text': CanvasText,
    'msui-text-muted': GrayText,
    'msui-line': CanvasText,
    'msui-accent': Highlight,
    'msui-accent-muted': withAlpha(Highlight, 0.4),
    'seq-track': Canvas,
    'seq-sequence': CanvasText,
    'seq-missing': GrayText,
    'seq-pending': GrayText,
    'seq-break': CanvasText,
    'seq-hover': Highlight,
    'seq-ruler': CanvasText,
    'seq-window': Highlight,
    'seq-grip': CanvasText,
    'seq-wash': withAlpha(Highlight, 0.2),
    'seq-cell-edge': CanvasText,
  } satisfies Partial<Record<SeqTokenName, string>>);
  return createTokens(tokens.theme, colors, tokens.fonts, tokens.logoChroma);
}

/** Whether the platform forces colors now. */
export function forcedColorsActive(win: Window | null | undefined): boolean {
  return !!win?.matchMedia?.('(forced-colors: active)')?.matches;
}

/**
 * The token table of a root's computed style. A value the CSS leaves as a function (`color-mix(…)`: a custom property
 * keeps its text) is resolved by the browser through a probe element's `color`.
 */
export function readCssTokens(root: Element): SeqTokens {
  const doc = root.ownerDocument;
  const win = doc.defaultView;
  const theme = themeOf(root);
  const style = win?.getComputedStyle(root);
  const fallback = builtInTokens(theme === 'dark' || style?.colorScheme === 'dark' ? 'dark' : 'light');
  if (!style || !win) return fallback;

  let probe: HTMLElement | undefined;
  const resolve = (value: string): string | undefined => {
    const direct = cssColorToHex(value);
    if (direct || !value.trim()) return direct;
    probe ??= root.appendChild(doc.createElement('span'));
    probe.style.color = '';
    probe.style.color = `color-mix(in srgb, ${value} 100%, transparent)`;
    if (!probe.style.color) return undefined;
    return cssColorToHex(win.getComputedStyle(probe).color);
  };

  const colors = {} as Record<SeqTokenName, string>;
  for (const name of SeqTokenNames)
    colors[name] = resolve(style.getPropertyValue(`--${name}`)) ?? colorHex(fallback.colors[name]);
  // under forced colors the system colors, resolved by the browser through the probe's `color`
  let system: SystemColors | undefined;
  if (forcedColorsActive(win)) {
    const sys = (name: keyof SystemColors) => {
      probe ??= root.appendChild(doc.createElement('span'));
      probe.style.color = '';
      probe.style.color = name;
      return cssColorToHex(win.getComputedStyle(probe).color) ?? DefaultSystemColors[name];
    };
    system = {
      Canvas: sys('Canvas'),
      CanvasText: sys('CanvasText'),
      Highlight: sys('Highlight'),
      GrayText: sys('GrayText'),
    };
  }
  probe?.remove();

  const font = (name: string, def: string) => style.getPropertyValue(`--msui-font-${name}`).trim() || def;
  const fonts: FontStacks = {
    sans: font('sans', fallback.fonts.sans),
    mono: font('mono', fallback.fonts.mono),
  };
  // `--seq-logo-chroma` is a number, not a color
  const chroma = Number.parseFloat(style.getPropertyValue('--seq-logo-chroma'));
  const tokens = createTokens(
    theme,
    colors,
    fonts,
    Number.isFinite(chroma) && chroma >= 0 ? chroma : fallback.logoChroma,
  );
  return system ? forcedColorTokens(tokens, system) : tokens;
}
