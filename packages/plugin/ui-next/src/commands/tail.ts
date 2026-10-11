/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The words typed after a command's alias, read into its fields.
import { type FieldValue, type FormFieldView, fits, isStringField, readValue } from '../forms/flat-fields.js';
import type { CommandFormats } from './types.js';

export interface TailReading {
  /** Values the tail gave, by field id (each fits its field). */
  values: Record<string, FieldValue>;
  /** Words it did not place: typed into the first empty field. */
  rest: string[];
}

/**
 * Reads a tail into fields. `initial` holds values already set and `defaulted` which of them are only defaults (typed
 * words win over those). Three passes: lead words, named words, the rest.
 */
export function readTail<Ctx>(input: {
  fields: readonly FormFieldView[];
  tail: string;
  initial: Readonly<Record<string, FieldValue>>;
  defaulted: ReadonlySet<string>;
  formats: CommandFormats<Ctx>;
  ctx: Ctx;
}): TailReading {
  const { fields, initial, formats, ctx } = input;
  const defaulted = new Set(input.defaulted);
  const values: Record<string, FieldValue> = {};
  const open = (f: FormFieldView) => (initial[f.id] === undefined || defaulted.has(f.id)) && values[f.id] === undefined;

  // Words fill the first empty field that names them exactly: an option ("sticks"), or a value of the field's format
  // that it lists ("blue", "polymer", "by element"; "element symbol" reads as one). Format names count only when every
  // word is named ("spacefill blue"); in a phrase ("the P-loop and the gatekeeper") the words stay text.
  const named = (f: FormFieldView, phrase: string, withFormats: boolean): string | undefined => {
    const w = phrase.toLowerCase();
    const option = f.options?.find((o) => o.value.toLowerCase() === w || o.label.toLowerCase() === w);
    if (option || !withFormats || !f.format) return option?.value;
    const format = formats.get(f.format);
    const listed = format?.complete(phrase, ctx) ?? [];
    const exact = listed.find((c) => c.value.toLowerCase() === w || c.label.toLowerCase() === w);
    if (exact) return exact.value;
    const read = format?.check(phrase, ctx);
    const all = read?.ok ? (format?.complete('', ctx) ?? []) : [];
    return read?.ok && all.some((c) => c.value === read.value) ? read.value : undefined;
  };
  /** A phrase as a field's value: one of its options, or what its format reads. */
  const read = (f: FormFieldView, phrase: string): string | undefined => {
    const w = phrase.toLowerCase();
    const option = f.options?.find((o) => o.value.toLowerCase() === w || o.label.toLowerCase() === w);
    if (option) return option.value;
    if (f.options) return undefined;
    const r = f.format ? formats.get(f.format)?.check(phrase, ctx) : undefined;
    return r?.ok ? String(r.value) : f.format ? undefined : phrase;
  };
  const assign = (words: string[], withFormats: boolean) => {
    const taken: Record<string, FieldValue> = {};
    const left: string[] = [];
    const at: number[] = [];
    for (let i = 0; i < words.length; ) {
      let step = 0;
      for (let n = Math.min(4, words.length - i); n >= 1 && !step; n--) {
        const phrase = words.slice(i, i + n).join(' ');
        // Typed words win over defaults (the selection): a defaulted field can still be named.
        for (const f of fields) {
          if (!open(f) || taken[f.id] !== undefined || !isStringField(f)) continue;
          const v = named(f, phrase, withFormats);
          if (v === undefined || !fits(f, v)) continue;
          taken[f.id] = v;
          step = n;
          break;
        }
      }
      if (step) i += step;
      else {
        at.push(i);
        left.push(words[i++] as string);
      }
    }
    return { taken, left, at };
  };

  // 1. A word a field is written after ("onto", "by", "on") starts its value: the words up to the next such word,
  // taken when they read as an option or in the field's format ("struct.oxy onto struct.deoxy by tm-align on A, B");
  // otherwise they stay for what follows.
  const leads = new Map<string, FormFieldView>();
  for (const f of fields)
    if (f.before && !/\s/.test(f.before) && isStringField(f)) leads.set(f.before.toLowerCase(), f);
  const typed = input.tail.split(/\s+/).filter(Boolean);
  const words: string[] = [];
  for (let i = 0; i < typed.length; ) {
    const f = leads.get((typed[i] as string).toLowerCase());
    let end = i + 1;
    while (end < typed.length && !leads.has((typed[end] as string).toLowerCase())) end++;
    const phrase = typed.slice(i + 1, end).join(' ');
    const value = f && open(f) && phrase ? read(f, phrase) : undefined;
    if (f && value !== undefined && fits(f, value)) {
      values[f.id] = value;
      defaulted.delete(f.id);
      i = end;
    } else words.push(typed[i++] as string);
  }

  // 2. Named words.
  const all = assign(words, true);
  // 3. The words left over, when they are one run ("A:315" in "A:315 element symbol"), may read as the value of a
  // field still open (then everything is placed); otherwise they stay text.
  const run = all.at.every((n, k) => k === 0 || n === (all.at[k - 1] as number) + 1);
  const leftText = all.left.join(' ');
  let placed = false;
  if (all.left.length && run)
    for (const f of fields) {
      if (!f.format || all.taken[f.id] !== undefined || !open(f)) continue;
      const r = formats.get(f.format)?.check(leftText, ctx);
      if (r?.ok && fits(f, r.value)) {
        all.taken[f.id] = r.value;
        placed = true;
        break;
      }
    }
  const { taken, left } = !all.left.length || placed ? { ...all, left: [] } : assign(words, false);
  Object.assign(values, taken);
  // Every value as its field holds it (a boolean's "on" is `true`).
  for (const f of fields) {
    if (values[f.id] !== undefined) values[f.id] = readValue(f, values[f.id]) as FieldValue;
  }
  return { values, rest: left };
}
