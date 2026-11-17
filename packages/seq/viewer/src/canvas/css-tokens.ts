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
  probe?.remove();

  const font = (name: string, def: string) => style.getPropertyValue(`--msui-font-${name}`).trim() || def;
  const fonts: FontStacks = {
    sans: font('sans', fallback.fonts.sans),
    mono: font('mono', fallback.fonts.mono),
  };
  return createTokens(theme, colors, fonts);
}
