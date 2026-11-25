/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A spec's colors as spelled (_spec/README.md "Themes", Colors): reading one needs no color space, so the spec reader's
// checks import this alone (Mol*'s `Color` and its spaces stay out of reading a spec).

/** What resolving a color reads: the table's theme (`light`, `dark`) and its colors by name. */
export interface TokenTable {
  readonly theme: string;
  readonly colors: Readonly<Record<string, number>>;
  /** Equal tables have equal versions. */
  readonly version: number;
}

/** `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa` as `0xRRGGBBAA`; undefined for anything else. */
export function parseColor(s: string | undefined): number | undefined {
  if (typeof s !== 'string' || !/^#([0-9a-f]{3,4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(s)) return undefined;
  let h = s.slice(1);
  if (h.length <= 4) h = [...h].map((c) => c + c).join('');
  if (h.length === 6) h += 'ff';
  return Number.parseInt(h, 16) >>> 0;
}

const LightDark = /^light-dark\(\s*([^,()\s]+)\s*,\s*([^,()\s]+)\s*\)$/;

/**
 * A spec's color: a literal (`#rrggbb`), a token (`token:accent` reads `seq-accent`, then `msui-accent`, then `accent`),
 * or `light-dark(a, b)` (`a` in a light table, `b` in a dark one; each a literal or a token). Undefined when it is none
 * of these, or a token the table lacks.
 */
export function resolveColor(value: unknown, tokens: TokenTable): number | undefined {
  if (typeof value !== 'string') return undefined;
  if (value.startsWith('token:')) {
    const name = value.slice(6);
    const colors = tokens.colors;
    return colors[`seq-${name}`] ?? colors[`msui-${name}`] ?? colors[name];
  }
  if (value.startsWith('light-dark(')) {
    const m = LightDark.exec(value);
    return m ? resolveColor(tokens.theme === 'dark' ? m[2] : m[1], tokens) : undefined;
  }
  return parseColor(value);
}

/** Whether a string is spelled as a color (a literal, a token name, `light-dark` of two): the check's rule. */
export function isColorSpelling(value: string): boolean {
  if (value.startsWith('token:')) return /^token:[a-z][a-z0-9-]*$/.test(value);
  const m = LightDark.exec(value);
  if (m) return isColorSpelling(m[1]) && isColorSpelling(m[2]) && !m[1].startsWith('light-dark');
  return parseColor(value) !== undefined;
}
