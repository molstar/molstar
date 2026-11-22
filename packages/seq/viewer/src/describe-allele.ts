/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesList } from '@molstar/seqx-core/data';

// A variant allele as its card says it (_spec/README.md "Behaviors", hoverCard, A variant allele's card): the item
// card's words for a source row holding `ref` and `alt` (a VCF's).

/** The INFO fields an allele's card lists before "+k more" (the rest in the source's rows). */
export const AlleleInfoRows = 8;

/** A source's fixed fields (core _spec "Data", VCF, Fields): what is not an INFO key. */
const FixedFields = new Set(['label', 'name', 'ref', 'alt', 'qual', 'filter', 'category']);

/**
 * A variant allele: a row of a source's list holding `ref` and `alt` is one; its record's ID (`label`), its category and
 * its INFO fields (every other field with a value, in the source's order); undefined for any other row.
 */
export function alleleOf(
  t: RangesList,
  i: number,
): { readonly id: string; readonly category: string; readonly info: [string, string][] } | undefined {
  const fields = t.fields ?? [];
  if (!fields.some((f) => f.name === 'ref') || !fields.some((f) => f.name === 'alt')) return undefined;
  const text = (name: string): string => {
    const f = fields.find((x) => x.name === name);
    return f ? cellText(f.column, i) : '';
  };
  const info: [string, string][] = [];
  for (const f of fields) {
    if (FixedFields.has(f.name)) continue;
    const v = cellText(f.column, i);
    if (v) info.push([f.name, v]);
  }
  return { id: text('label'), category: text('category'), info };
}

/**
 * A column's cell as text: '' when missing (a null category, a NaN, int32's missing value); an integer as it is (an
 * INFO integer is often an id, ALLELEID: never grouped), a fraction to three significant decimals.
 */
function cellText(c: Column, i: number): string {
  if (c instanceof CategoryColumn) return c.get(i) ?? '';
  const v = c[i];
  if (v === undefined || Number.isNaN(v) || v === -2147483648) return '';
  if (Number.isInteger(v) || !Number.isFinite(v)) return String(v);
  const a = Math.abs(v);
  return Number(v.toFixed(a >= 100 ? 1 : a >= 1 ? 2 : 3)).toString();
}
