/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorNames } from '@molstar/core/util/color/names';
import { activeElement } from '../../dom/events.js';
import { bind, effect, type MaybeReactive } from '../../dom/scope.js';
import { drag } from '../../dom/scoped.js';
import { observableState } from '../../reactive/observable-state.js';
import type { ReadonlyState } from '../../reactive/readonly-state.js';
import { derive, isState } from '../../reactive/state.js';
import { classes } from './_generated/color-picker.module.css.js';

/** The palette's hues (vivid; their lighter and darker weights are the row under them, `colorWeights`). */
export const ColorPalette: readonly number[] = [
  0xe5484d, 0xf76b15, 0xffc53d, 0x46a758, 0x12a594, 0x0090ff, 0x6e56cf, 0xd6409f,
];

/**
 * The weights of a color, lightest to darkest: eight, itself among them, placed by how light it is (black is the last,
 * white the first, a mid color about the middle), the rest spread over the room it has toward white and toward black.
 */
export function colorWeights(c: number): number[] {
  const r = (c >> 16) & 0xff;
  const g = (c >> 8) & 0xff;
  const b = c & 0xff;
  const lightness = (Math.max(r, g, b) + Math.min(r, g, b)) / 510;
  const lighter = Math.round(7 * (1 - lightness));
  const darker = 7 - lighter;
  const weights: number[] = [];
  for (let i = lighter; i > 0; i--) weights.push(shade(c, (0.9 * i) / lighter));
  weights.push(c);
  for (let i = 1; i <= darker; i++) weights.push(shade(c, (-0.85 * i) / darker));
  return weights;
}

/** For CSS: `#rrggbb`. */
export function colorToHex(c: number): string {
  return `#${((c >>> 0) & 0xffffff).toString(16).padStart(6, '0')}`;
}

/** For people: the channels in decimal, `51, 98, 178`. */
export function colorToRgb(c: number): string {
  return `${(c >> 16) & 0xff}, ${(c >> 8) & 0xff}, ${c & 0xff}`;
}

const rgb = (r: number, g: number, b: number) => (r << 16) | (g << 8) | b;
const byte = (v: number) => Math.max(0, Math.min(255, Math.round(v)));

/**
 * A color as people write it, to the number Mol* uses (`0xRRGGBB`); undefined when it is not one:
 * - hex: `#3362b2`, `3362b2`, `#36b`, `0x3362b2`;
 * - channels, 0–255 or percentages: `51, 98, 178`, `51 98 178`, `rgb(51 98 178)`, `rgba(51,98,178,1)`, `rgb(20%,…)`;
 * - normalized channels (all at most 1, one with a decimal point): `0.2, 0.38, 0.7`;
 * - `hsl(218, 55%, 45%)`;
 * - a CSS color name: `steelblue`, `Steel Blue`.
 */
export function parseColor(text: string): number | undefined {
  const t = text.trim().toLowerCase();
  if (!t) return undefined;

  const hex = /^(?:#|0x)?([0-9a-f]{3}|[0-9a-f]{6})$/.exec(t);
  // Three bare digits are channels' worth of hex only without separators ("123" is #112233, "1, 2, 3" is not).
  if (hex) {
    const h = hex[1].length === 3 ? hex[1].replace(/./g, (d) => d + d) : hex[1];
    return Number.parseInt(h, 16);
  }

  const fn = /^(rgba?|hsla?)?\s*\(?\s*([^()]*?)\s*\)?$/.exec(t);
  if (fn) {
    // Separated by commas, spaces, semicolons or a slash (the alpha in `rgb(r g b / a)`, ignored).
    const parts = fn[2].split(/\s*[,;/\s]\s*/).filter(Boolean);
    if (parts.length === 3 || parts.length === 4) {
      const values = parts.slice(0, 3);
      if (values.every((p) => /^[+]?(\d+\.?\d*|\.\d+)%?$/.test(p))) {
        if (fn[1]?.startsWith('hsl')) return hsl(values);
        const percent = values.some((p) => p.endsWith('%'));
        const nums = values.map((p) => Number.parseFloat(p));
        const normalized = !percent && nums.every((v) => v <= 1) && values.some((p) => p.includes('.'));
        const to = (v: number) => (percent ? (v * 255) / 100 : normalized ? v * 255 : v);
        if (!percent && !normalized && nums.some((v) => v > 255)) return undefined;
        return rgb(byte(to(nums[0])), byte(to(nums[1])), byte(to(nums[2])));
      }
    }
  }

  const name = t.replace(/[\s_-]/g, '') as keyof typeof ColorNames;
  return Object.hasOwn(ColorNames, name) ? (ColorNames[name] as number) : undefined;
}

function hsl(values: string[]): number | undefined {
  const h = (((Number.parseFloat(values[0]) % 360) + 360) % 360) / 60;
  const s = Math.min(1, Number.parseFloat(values[1]) / 100);
  const l = Math.min(1, Number.parseFloat(values[2]) / 100);
  const c = (1 - Math.abs(2 * l - 1)) * s;
  const x = c * (1 - Math.abs((h % 2) - 1));
  const m = l - c / 2;
  const [r, g, b] =
    h < 1 ? [c, x, 0] : h < 2 ? [x, c, 0] : h < 3 ? [0, c, x] : h < 4 ? [0, x, c] : h < 5 ? [x, 0, c] : [c, 0, x];
  return rgb(byte((r + m) * 255), byte((g + m) * 255), byte((b + m) * 255));
}

/** Mixes a color toward white (`amount` > 0) or black (< 0), per channel. */
export function shade(c: number, amount: number): number {
  const target = amount > 0 ? 255 : 0;
  const t = Math.abs(amount);
  const ch = (shift: number) => {
    const v = (c >> shift) & 0xff;
    return Math.round(v + (target - v) * t) << shift;
  };
  return ch(16) | ch(8) | ch(0);
}

// The smooth palette: hues left to right (as the CSS gradient draws them, piecewise linear in RGB), mixed toward white
// above the middle and toward black below it, at most this much (white and black are on the greys under it).
const FieldEdge = 0.9;

/** The hue at `x` (0–1) at full strength. */
function hueAt(x: number): number {
  const h = (x * 6) % 6;
  const f = h - Math.floor(h);
  const up = byte(f * 255);
  const down = 255 - up;
  switch (Math.floor(h)) {
    case 0:
      return rgb(255, up, 0);
    case 1:
      return rgb(down, 255, 0);
    case 2:
      return rgb(0, 255, up);
    case 3:
      return rgb(0, down, 255);
    case 4:
      return rgb(up, 0, 255);
    default:
      return rgb(255, 0, down);
  }
}

/** The smooth palette's color at (`x`, `y`), both 0–1: what the gradients show there. */
export function fieldColor(x: number, y: number): number {
  const hue = hueAt(Math.max(0, Math.min(1, x)));
  const t = Math.max(0, Math.min(1, y));
  return t < 0.5 ? shade(hue, FieldEdge * (1 - t / 0.5)) : shade(hue, -FieldEdge * ((t - 0.5) / 0.5));
}

/** Where a color sits: on the field (its hue across, its lightness down) or, without a hue to speak of, on the greys. */
function placeOf(c: number): { greys: boolean; x: number; y: number } {
  const r = (c >> 16) & 0xff;
  const g = (c >> 8) & 0xff;
  const b = c & 0xff;
  const max = Math.max(r, g, b);
  const min = Math.min(r, g, b);
  if (max - min < 16) return { greys: true, x: (r + g + b) / 765, y: 0.5 };
  const d = max - min;
  const h = max === r ? ((g - b) / d + 6) % 6 : max === g ? (b - r) / d + 2 : (r - g) / d + 4;
  const l = (max + min) / 510;
  const y = l >= 0.5 ? 0.5 * (1 - (2 * l - 1) / FieldEdge) : 0.5 + (0.5 * (1 - 2 * l)) / FieldEdge;
  return { greys: false, x: h / 6, y: Math.max(0, Math.min(1, y)) };
}

export interface ColorPickerProps {
  value: MaybeReactive<number>;
  onChange(color: number): void;
}

/**
 * A color picker that looks and works the same everywhere (no system picker): the smooth palette (every hue, lighter
 * up, darker down; greys under it) to click or drag across, the palette's hues, the current color marked; then, in two
 * steps, the weights of the color picked (lighter to darker; picking one keeps the row, so it can be walked); a field
 * that takes a color in most ways people write one (`parseColor`) and shows it as decimal channels.
 */
export function ColorPicker(props: ColorPickerProps) {
  const value = held(props.value);
  const place = derive(value, placeOf);
  let current = 0;
  let input: HTMLInputElement | undefined;
  // The color the weights are of: every pick but a weight's (and any change from outside) moves it.
  const base = observableState(0);
  let fromWeight = false;
  const weights = derive(base.asReadonly, colorWeights);
  bind(value, (v) => {
    current = v;
    if (!fromWeight) base.set(v);
    // Not under the caret: what is being typed stays as typed until it is left (then it is written out again).
    if (input && activeElement(input) !== input) input.value = colorToRgb(v);
  });

  // Click or drag across the field or the greys (`drag`, dom/scoped.ts: attached once, the pointer captured, so a drag
  // can leave the strip).
  const strip = (pick: (x: number, y: number) => number) => (el: HTMLElement) => {
    const at = (e: PointerEvent) => {
      const r = el.getBoundingClientRect();
      props.onChange(pick((e.clientX - r.left) / r.width, (e.clientY - r.top) / r.height));
    };
    drag(el, {
      start(e) {
        if (e.button !== 0) return undefined;
        at(e);
        return true;
      },
      move: at,
    });
  };
  // Snapping at the ends: black and white exactly.
  const grey = (x: number) => {
    const v = x < 0.03 ? 0 : x > 0.97 ? 255 : byte(x * 255);
    return rgb(v, v, v);
  };

  return (
    <div class={classes.picker}>
      <div class={classes.field} aria-hidden="true" ref={strip(fieldColor)}>
        <span
          class={classes.marker}
          hidden={derive(place, (p) => p.greys)}
          style={{ left: derive(place, (p) => `${p.x * 100}%`), top: derive(place, (p) => `${p.y * 100}%`) }}
        />
      </div>
      <div class={classes.greys} aria-hidden="true" ref={strip((x) => grey(x))}>
        <span
          class={classes.marker}
          hidden={derive(place, (p) => !p.greys)}
          style={{ left: derive(place, (p) => `${p.x * 100}%`), top: '50%' }}
        />
      </div>
      <div class={classes.palette} role="listbox" aria-label="Palette">
        {ColorPalette.map((c) => (
          <button
            type="button"
            role="option"
            class={classes.swatch}
            style={{ background: colorToHex(c) }}
            aria-label={colorToRgb(c)}
            aria-selected={derive(value, (v) => v === c)}
            onClick={() => props.onChange(c)}
          />
        ))}
      </div>
      <div class={classes.palette} role="listbox" aria-label="Weights">
        {[0, 1, 2, 3, 4, 5, 6, 7].map((i) => {
          const color = derive(weights, (w) => w[i]);
          return (
            <button
              type="button"
              role="option"
              class={classes.swatch}
              style={{ background: derive(color, colorToHex) }}
              aria-label={derive(color, colorToRgb)}
              aria-selected={derive([color, value], (c, v) => c === v)}
              onClick={() => {
                fromWeight = true;
                props.onChange(colorWeights(base.value)[i]);
                fromWeight = false;
              }}
            />
          );
        })}
      </div>
      <div class={classes.custom}>
        <span class={classes.preview} style={{ background: derive(value, colorToHex) }} />
        <input
          ref={(el: HTMLInputElement) => {
            input = el;
            el.value = colorToRgb(current);
          }}
          type="text"
          class={classes.text}
          aria-label="Color"
          placeholder="r, g, b · #hex · rgb() · name"
          data-tip="r, g, b; hex, rgb(), hsl() or a name"
          spellcheck={false}
          onInput={(e) => {
            const c = parseColor(e.currentTarget.value);
            if (c !== undefined) props.onChange(c);
          }}
          onKeyDown={(e) => {
            if (e.key === 'Enter') e.currentTarget.value = colorToRgb(current);
          }}
          onBlur={(e) => {
            e.currentTarget.value = colorToRgb(current);
          }}
        />
      </div>
    </div>
  );
}

/** The value as a state: a state as it is; a stream (or a plain value) held in one, with the current scope. */
function held(source: MaybeReactive<number>): ReadonlyState<number> {
  if (isState(source)) return source as ReadonlyState<number>;
  const state = observableState(0);
  effect(source, (v) => state.set(v));
  return state.asReadonly;
}
