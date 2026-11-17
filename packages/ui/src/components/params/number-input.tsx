/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { round } from '@molstar/core/util';
import { getPrecision } from '@molstar/core/util/number';
import { effect, frameState } from '../../dom/index.js';

/** A number param's range (PD's `Range`): any of the bounds, the step. */
export interface NumberRange {
  min?: number;
  max?: number;
  step?: number;
}

/** The decimals kept and shown for a number without a step. */
const MaxDigits = 4;

/** A typed number: "0.5", "0,5" (a comma as the decimal point), " -1 "; NaN when it is not one. */
export function parseNumber(text: string): number {
  const t = text.trim().replace(',', '.');
  return t === '' ? Number.NaN : Number(t);
}

/** The decimals of a step (0.25 → 2, 2.5 → 1, 1 → 0); undefined without one. */
export function stepDigits(step?: number): number | undefined {
  return step !== undefined && step > 0 ? getPrecision(step) : undefined;
}

/**
 * A number as the field keeps it: rounded to its step's decimals (as v5's number field: the step's precision, not snapped
 * to a multiple of it; at most 4 decimals without a step), then clamped to whichever bounds exist.
 */
export function holdNumber(v: number, range: NumberRange = {}): number {
  let x = round(v, stepDigits(range.step) ?? MaxDigits);
  if (range.min !== undefined && x < range.min) x = range.min;
  if (range.max !== undefined && x > range.max) x = range.max;
  return x;
}

/**
 * A number's text: the step's decimals at least (0.5 with a step of 0.25 reads "0.50"), more when the value has them (up
 * to 4 past none), so a value set elsewhere reads as it is held; without a step, the value as it is (at most 4 decimals).
 */
export function formatNumber(v: number | undefined, step?: number): string {
  if (v === undefined || !Number.isFinite(v)) return '';
  const digits = stepDigits(step);
  if (digits === undefined) return String(round(v, MaxDigits));
  return v.toFixed(Math.max(digits, Math.min(MaxDigits, getPrecision(round(v, MaxDigits)))));
}

/** An input's value ('' for none), written only when it differs (a write moves the caret). */
export function setValue(el: HTMLInputElement | HTMLTextAreaElement, value: unknown) {
  const v = value === null || value === undefined ? '' : String(value);
  if (el.value !== v) el.value = v;
}

/**
 * A typed number: text, not a native number field (that one shows decimals in the browser's locale, "0,5", while every
 * other number reads "0.5"; a comma typed is read as the decimal point). A commit (change: blur, Enter) is held
 * (`holdNumber`) and set when it differs; the field then says what holds, whether or not it changed (a clamp to the same
 * value, "1.50" for 1.5), and text that is no number goes back to the value. The caller writes the value
 * (`setValue(el, formatNumber(v, step))` in its watch).
 */
export function NumberInput(props: {
  current(): number | undefined;
  commit(v: number): void;
  /** Its bounds and step; a function when they move (an interval's end: up to the other end). */
  range?: NumberRange | (() => NumberRange);
  class?: string;
  'aria-label'?: string;
  'aria-labelledby'?: string;
  'data-tip'?: string;
}): HTMLInputElement {
  return (
    <input
      class={props.class}
      type="text"
      inputMode="decimal"
      aria-label={props['aria-label']}
      aria-labelledby={props['aria-labelledby']}
      data-tip={props['data-tip']}
      onChange={(e) => {
        const input = e.currentTarget;
        const range = typeof props.range === 'function' ? props.range() : props.range;
        const n = parseNumber(input.value);
        const current = props.current();
        const v = Number.isFinite(n) ? holdNumber(n, range) : current;
        if (v !== undefined && Number.isFinite(v) && !Object.is(v, current)) props.commit(v);
        input.value = formatNumber(v, range?.step);
      }}
    />
  ) as HTMLInputElement;
}

/**
 * Input published in the frame lane (dom §3: input events come several a frame; a set may update a structure's
 * representation): `offer` keeps the latest, the frame sets it unless the field already `holds` it. Only what was
 * offered since the last frame is published, so a drag back to the last dragged value after a change from elsewhere (a
 * typed value, a reset) is set again, and a resume sets nothing stale.
 */
export function frameInput<T>(set: (v: T) => void, holds: (v: T) => boolean): (v: T) => void {
  let pending = false;
  // Never equal: each frame's offer is published; whether it is new is the field's to say, not the last offer's.
  const offered = frameState<T | undefined>(undefined, () => false);
  effect(offered, (v) => {
    if (!pending) return;
    pending = false;
    if (!holds(v as T)) set(v as T);
  });
  return (v) => {
    pending = true;
    offered.offer(v);
  };
}
