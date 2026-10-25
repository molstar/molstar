/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The skin (seqspec:interaction#the-skin): `skin.css` declares the `--seq-*` tokens at zero specificity in three
// blocks, light, and two dark ones that are identical and define every token the light one does; the text tokens meet
// WCAG AA in both schemes, `--seq-conflict` (error text, the error badge, danger items) at least 4.5:1 on the surface.
// The chrome follows Prism's canvas-chrome rules: glass, no 1px borders, the picked item an accent wash, icons for the
// check and the chevron, and its optional tokens derived where they are used.

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { DarkSkin, DefaultSkin } from '../src/index.js';

const css = readFileSync(new URL('../src/skin.css', import.meta.url), 'utf8').replace(/\r\n/g, '\n');

/** The declarations of the first rule whose selector is `selector`, in order. */
function block(selector: string): [string, string][] {
  const at = css.indexOf(`${selector} {`);
  if (at < 0) throw new Error(`no rule ${selector}`);
  const body = css.slice(css.indexOf('{', at) + 1, css.indexOf('}', at));
  return [...body.matchAll(/(--seq-[\w-]+):\s*([^;]+);/g)].map((m) => [m[1], m[2].trim()]);
}
const light = block(':where(:root)');
const dark = block(
  ':where(:root[data-theme="dark"], :root:not([data-theme="light"]):is([data-msp-theme="dark"], [data-msp-theme="blue"]))',
);
const auto = block(':where(:root:not([data-theme="light"], [data-msp-theme="light"]))');

const token = (b: [string, string][], name: string) => b.find(([k]) => k === name)?.[1] as string;
const rgb = (v: string) => (v.match(/[\d.]+/g) ?? []).slice(0, 3).map(Number);
/** WCAG 2 relative luminance of an opaque `rgb(r g b)`. */
function luminance(v: string) {
  const [r, g, b] = rgb(v).map((c) => {
    const s = c / 255;
    return s <= 0.03928 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}
const contrast = (a: string, b: string) => {
  const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return (hi + 0.05) / (lo + 0.05);
};

describe('skin.css tokens', () => {
  it('the dark blocks are identical to each other and define every token the light block does', () => {
    expect(auto).toEqual(dark);
    expect(dark.map(([k]) => k).sort()).toEqual(light.map(([k]) => k).sort());
  });

  it('every token the skin names is declared, with the conflict token in both dark blocks', () => {
    const names = [
      'surface',
      'surface-alt',
      'ink',
      'ink-quiet',
      'ink-faint',
      'grid',
      'accent',
      'hover',
      'highlight',
      'selection',
      'selection-border',
      'unobserved',
      'conflict',
      'pending',
      'grip',
      'grip-hover',
      'window',
      'focus-ring',
      'font-mono',
      'font-ui',
    ].map((n) => `--seq-${n}`);
    for (const b of [light, dark, auto]) expect(names.filter((n) => token(b, n) === undefined)).toEqual([]);
    expect(token(dark, '--seq-conflict')).not.toBe(token(light, '--seq-conflict'));
  });

  it('is declared at zero specificity: `:where` on every token rule', () => {
    const selectors = [...css.replace(/\/\*[\s\S]*?\*\//g, '').matchAll(/([^{};]+)\{/g)].map((m) => m[1].trim());
    const tokenRules = selectors.filter((s) => s.includes(':root'));
    expect(tokenRules.length).toBe(3);
    for (const s of tokenRules) expect(s.startsWith(':where(')).toBe(true);
  });

  for (const [scheme, b] of [
    ['light', light],
    ['dark', dark],
  ] as const) {
    it(`${scheme}: text meets WCAG AA on the surface, conflict at least 4.5:1, marks at least 3:1`, () => {
      const on = (name: string, ground = '--seq-surface') => contrast(token(b, name), token(b, ground));
      expect(on('--seq-ink')).toBeGreaterThanOrEqual(7);
      expect(on('--seq-ink', '--seq-surface-alt')).toBeGreaterThanOrEqual(4.5);
      expect(on('--seq-ink-quiet')).toBeGreaterThanOrEqual(4.5);
      expect(on('--seq-conflict')).toBeGreaterThanOrEqual(4.5);
      expect(on('--seq-accent')).toBeGreaterThanOrEqual(4.5);
      expect(on('--seq-focus-ring')).toBeGreaterThanOrEqual(3);
      expect(on('--seq-grip')).toBeGreaterThanOrEqual(3);
      expect(on('--seq-grip-hover')).toBeGreaterThan(on('--seq-grip'));
    });
  }

  it("the reader's fallbacks are skin.css's blocks: DefaultSkin the light one, DarkSkin the dark one (UX1's routed items)", () => {
    const camel = (k: string) => k.replace('--seq-', '').replace(/-(\w)/g, (_, c: string) => c.toUpperCase());
    for (const [skin, b] of [
      [DefaultSkin, light],
      [DarkSkin, dark],
    ] as const) {
      const fromCss = Object.fromEntries(b.map(([k, v]) => [camel(k), v]));
      const colors = Object.entries(skin).filter(
        ([k, v]) => k !== 'dark' && !k.startsWith('font') && typeof v === 'string',
      );
      expect(Object.fromEntries(colors)).toEqual(
        Object.fromEntries(Object.entries(fromCss).filter(([k]) => !k.startsWith('font'))),
      );
    }
    expect([DefaultSkin.dark, DarkSkin.dark]).toEqual([false, true]);
  });

  it('the old dark skin did not meet it: the light conflict on the dark surface is below 4.5:1', () => {
    expect(contrast(token(light, '--seq-conflict'), token(dark, '--seq-surface'))).toBeLessThan(3);
  });

  /** The declarations of the rule whose selector list is exactly `selector`. */
  const rule = (selector: string) => {
    const at = css.indexOf(`\n${selector} {`);
    if (at < 0) throw new Error(`no rule ${selector}`);
    return css.slice(at, css.indexOf('}', at));
  };
  const rules = () => [...css.replace(/\/\*[\s\S]*?\*\//g, '').matchAll(/([^{}]+)\{([^{}]*)\}/g)];

  it('the chrome draws no 1px borders: a picked item is an accent wash with accent ink, a swatch a hairline shadow', () => {
    const borders = rules()
      .filter(([, , body]) => /border(?!-radius)(-\w+)?:\s*1px/.test(body))
      .map(([, sel]) => sel.trim());
    expect(borders).toEqual([]);
    const picked = rule(
      '.seqv-tab[aria-selected="true"],\n.seqv-choice[aria-pressed="true"],\n.seqv-dot[aria-current]',
    );
    expect(picked).toContain('var(--seqv-selected-ink)');
    expect(picked).toContain('var(--seqv-selected)');
    expect(rule('.seqv-swatch')).toContain('box-shadow: var(--seqv-hairline');
  });

  it('what floats over the canvas is glass with a hairline; only what opens for a moment adds a shadow', () => {
    const glass = rule('.seqv-tooltip,\n.seqv-hover-card');
    expect(glass).toContain('var(--seqv-glass-panel)');
    expect(glass).toContain('box-shadow: var(--seqv-hairline);');
    expect(glass).toContain('var(--seqv-radius)');
    expect(rule('.seqv-issue-list,\n.seqv-menu-list')).toContain('var(--seqv-hairline), var(--seqv-menu-shadow)');
    // No literal drop shadow anywhere but the overlay's default.
    const shadows = [...css.matchAll(/box-shadow:[^;]*rgb\(0 0 0[^;]*;/g)];
    expect(shadows).toEqual([]);
  });

  it('no text glyph for a check or a chevron, no ellipsis: icons and fades', () => {
    expect(css).not.toMatch(/content:\s*"(✓|▾)/);
    expect(css).not.toContain('text-overflow: ellipsis');
    expect(css).toContain('animation-timeline: scroll(self inline)');
  });

  it('each optional chrome token has a default where it is used, derived from the tokens a host maps', () => {
    for (const t of [
      'chrome',
      'chrome-panel',
      'chrome-filter',
      'chrome-edge',
      'menu-shadow',
      'control',
      'selected',
      'selected-ink',
      'font-size',
      'font-size-small',
      'radius',
      'radius-control',
      'chrome-gap',
      'canvas-edge',
    ])
      expect([t, new RegExp(`var\\(--seq-${t}, [^)]`).test(css)]).toEqual([t, true]);
    // Never declared on the root: a value there would be computed from the root's tokens, not the host's.
    expect([...light, ...dark].map(([k]) => k).filter((k) => /--seq-(chrome|control|selected|radius)/.test(k))).toEqual(
      [],
    );
  });

  it("the chrome's text column is the canvas's: `--seqv-label-x` is where the composer draws a track's label", () => {
    const composer = readFileSync(new URL('../src/draw/composer-headers.ts', import.meta.url), 'utf8');
    const drawn = [...composer.matchAll(/const x = (\d+) \+ (?:inset|indent)/g)].map((m) => Number(m[1]));
    expect(drawn.length).toBeGreaterThan(1);
    const css = /--seqv-label-x: (\d+)px;/.exec(readFileSync(new URL('../src/skin.css', import.meta.url), 'utf8'));
    expect(new Set([...drawn, Number(css?.[1])])).toEqual(new Set([drawn[0]]));
  });
});
