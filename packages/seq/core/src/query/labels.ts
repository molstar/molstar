/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Labels (_spec/README.md "Queries", Numberings and labels): positions are 0-based half-open indices inside, labels
// are inclusive names of positions at the edge, converted once where a form is read and never mixed in one table.
// A label is an integer with an optional suffix: letters (an insertion code, `94A`) or `.n` (`100.1`); any other word
// is a name (a numbering of names). Labels have one order, so a label a numbering lacks is still placed. A `.n` label
// is written quoted in a form (`(at :label "100.1")`): bare, `100.1` is the number 100.1, not a label.

import { compareStrings } from './compare.js';

/** A label read: its integer, its suffix kind (none, letters, `.n`) and suffix; a name has no integer. */
export interface ParsedLabel {
  readonly int?: number;
  /** 0 none, 1 letters, 2 `.n`; names sort after every numbered label. */
  readonly kind: 0 | 1 | 2;
  readonly suffix: string;
}

const LabelPattern = /^([+-]?\d+)(?:([A-Za-z]+)|\.(\d+))?$/;

/** A label's integer, suffix kind and suffix (a name: no integer, the whole label as its suffix). */
export function parseLabel(label: string): ParsedLabel {
  const m = LabelPattern.exec(label);
  if (!m) return { kind: 0, suffix: label };
  if (m[2] !== undefined) return { int: Number(m[1]), kind: 1, suffix: m[2] };
  if (m[3] !== undefined) return { int: Number(m[1]), kind: 2, suffix: m[3] };
  return { int: Number(m[1]), kind: 0, suffix: '' };
}

/**
 * The order of labels: by integer; then no suffix, then letters (shorter first, then alphabetic ignoring case), then
 * `.n` by n; names after every numbered label, alphabetic ignoring case. Case and spelling break the last ties, so the
 * order is total: `99 < 100 < 100A < 100Z < 100AA < 100.1 < 100.2 < 101`.
 */
export function compareLabels(a: string, b: string): number {
  if (a === b) return 0;
  const x = parseLabel(a);
  const y = parseLabel(b);
  if (x.int === undefined || y.int === undefined) {
    if (x.int !== undefined) return -1;
    if (y.int !== undefined) return 1;
    return compareStrings(x.suffix.toUpperCase(), y.suffix.toUpperCase()) || compareStrings(a, b);
  }
  if (x.int !== y.int) return x.int - y.int;
  if (x.kind !== y.kind) return x.kind - y.kind;
  if (x.kind === 1) {
    if (x.suffix.length !== y.suffix.length) return x.suffix.length - y.suffix.length;
    return compareStrings(x.suffix.toUpperCase(), y.suffix.toUpperCase()) || compareStrings(a, b);
  }
  if (x.kind === 2) return Number(x.suffix) - Number(y.suffix) || compareStrings(a, b);
  return compareStrings(a, b);
}

/** The largest magnitude of a label's integer: an int32, so a label never wraps. */
const Int32Max = 2 ** 31 - 1;

/** What `readLabel` gives: the label's text, or why the value is not a label. */
export type LabelRead = { readonly label: string } | { readonly error: 'syntax' | 'int32' | 'fraction' };

/**
 * A form's label value: an integer (`94`) or a word (`94A`, `"100.1"`, a name). A number that is not an integer is
 * not a label (`fraction`: a `.n` label is written quoted), nor is an empty word; a label whose integer is outside int32
 * is refused, never wrapped. A written integer keeps its spelling as JavaScript prints it (`094` is read as the number
 * 94, so the label `94`); a quoted label keeps its own (`"100.10"` is not `"100.1"`).
 */
export function readLabel(value: unknown): LabelRead {
  if (typeof value === 'number') {
    if (!Number.isInteger(value)) return { error: Number.isFinite(value) ? 'fraction' : 'syntax' };
    return Math.abs(value) > Int32Max ? { error: 'int32' } : { label: `${value}` };
  }
  if (typeof value !== 'string' || !value.trim() || value !== value.trim()) return { error: 'syntax' };
  const p = parseLabel(value);
  if (p.int !== undefined && Math.abs(p.int) > Int32Max) return { error: 'int32' };
  return { label: value };
}

/** The spellings of "no insertion code" in an insertion-code column: `""`, `?` and `.`. */
export function isNoInsertionCode(code: string): boolean {
  return code === '' || code === '?' || code === '.';
}

/** Whether a label has a suffix (an insertion code or `.n`), so another code may not be added to it. */
export function hasLabelSuffix(label: string): boolean {
  const p = parseLabel(label);
  return p.int !== undefined && p.kind !== 0;
}
