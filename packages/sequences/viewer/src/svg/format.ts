/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the SVG writer prints with (seqspec:viewer#the-render-backend): numbers as the shortest decimal at a fixed
// precision (device px over the DPR, in hundredths: never `toString` on a long double), XML escaping, and CSS colors as
// `#rrggbb` plus an opacity (SVG 1.1 readers, Inkscape and Illustrator, do not read `rgb(0 0 0 / 0.5)`). Platform-free.

const INTS: string[] = [];

/** A whole number of hundredths as a decimal: 1250 → "12.5", 300 → "3", -5 → "-0.05". */
export function dec(h: number): string {
  if (h % 100 !== 0) return String(h / 100);
  const v = h / 100;
  if (v < 0 || v >= 1024) return String(v);
  INTS[v] ??= String(v);
  return INTS[v];
}
/** `v` to three decimals (scales, opacities), shortest. */
export function dec3(v: number): string {
  return String(Math.round(v * 1000) / 1000);
}

/** XML text and attribute escaping; characters XML 1.0 cannot carry (controls but tab, LF, CR) are dropped. */
export function esc(s: string): string {
  let out = "";
  let from = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i);
    const r =
      c === 38
        ? "&amp;"
        : c === 60
          ? "&lt;"
          : c === 62
            ? "&gt;"
            : c === 34
              ? "&quot;"
              : c < 32 && c !== 9 && c !== 10 && c !== 13
                ? ""
                : null;
    if (r === null) continue;
    out += s.slice(from, i) + r;
    from = i + 1;
  }
  return from === 0 ? s : out + s.slice(from);
}

export const hex = (c: number): string => `#${(c & 0xffffff).toString(16).padStart(6, "0")}`;

const colors = new Map<string, readonly [string, number]>();
/** A CSS color as `#rrggbb` and its opacity; a keyword or function SVG 1.1 may not know is returned as given, opacity 1. */
export function cssColor(css: string): readonly [string, number] {
  let r = colors.get(css);
  if (!r) {
    if (colors.size > 512) colors.clear();
    r = parseColor(css);
    colors.set(css, r);
  }
  return r;
}
function parseColor(css: string): readonly [string, number] {
  const t = css.trim().toLowerCase();
  if (t[0] === "#") {
    const h = t.length === 4 || t.length === 5 ? [...t.slice(1)].map((c) => c + c).join("") : t.slice(1);
    if (/^[0-9a-f]{6}([0-9a-f]{2})?$/.test(h))
      return [`#${h.slice(0, 6)}`, h.length === 8 ? Number.parseInt(h.slice(6), 16) / 255 : 1];
  }
  const m = /^rgba?\(([^)]*)\)$/.exec(t);
  const p = m ? m[1].split(/[\s,/]+/).filter(Boolean) : [];
  if (p.length >= 3) {
    const ch = (v: string) =>
      Math.max(0, Math.min(255, Math.round(Number.parseFloat(v) * (v.endsWith("%") ? 2.55 : 1))));
    const a = p.length > 3 ? Number.parseFloat(p[3]) / (p[3].endsWith("%") ? 100 : 1) : 1;
    return [hex((ch(p[0]) << 16) | (ch(p[1]) << 8) | ch(p[2])), a >= 0 && a <= 1 ? a : 1];
  }
  return t === "transparent" ? ["#000000", 0] : [css, 1];
}
