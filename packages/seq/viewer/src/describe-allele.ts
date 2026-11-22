/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesList } from '@molstar/seqx-core/data';
import { FixedFields as VcfFixedFields } from '@molstar/seqx-core/runtimes/vcf-parse';

// A variant allele as its card says it (_spec/README.md "Behaviors", hoverCard, A variant allele's card): the item
// card's words for a source row holding `ref` and `alt` (a VCF's).

/** The INFO fields an allele's card lists before "+k more" (the rest in the source's rows). */
export const AlleleInfoRows = 8;

/** The longest INFO value a card row shows (a long list of ids, a description): the rest as an ellipsis. */
export const AlleleInfoChars = 40;

/** A source's fixed fields (core's VCF reader's: core _spec "Data", VCF, Fields): what is not an INFO key. */
const FixedFields = new Set(VcfFixedFields.map((f) => f.name));

/** A variant allele's card words: its record's ID, category, QUAL and FILTER, and its INFO fields. */
export interface Allele {
  readonly id: string;
  readonly category: string;
  readonly qual: string;
  readonly filter: string;
  readonly info: readonly (readonly [string, string])[];
  /** Its card's rows after the category: ID, QUAL, FILTER, the first `AlleleInfoRows` INFO fields and "+k more". */
  readonly rows: readonly (readonly [string, string])[];
}

/** The last allele read of a list (its row): a hover held over one reads it once; a list let go takes it along. */
const Last = new WeakMap<RangesList, { i: number; allele: Allele | undefined }>();

/**
 * A variant allele: a row of a source's list holding `ref` and `alt` is one; its record's ID (`label`), its category,
 * QUAL and FILTER, and its INFO fields (every other field with a value, in the source's order, each at most
 * `AlleleInfoChars`); undefined for any other row. The last one read is kept (a hover again makes nothing).
 */
export function alleleOf(t: RangesList, i: number): Allele | undefined {
  const last = Last.get(t);
  if (last?.i === i) return last.allele;
  const allele = readAllele(t, i);
  Last.set(t, { i, allele });
  return allele;
}

function readAllele(t: RangesList, i: number): Allele | undefined {
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
    if (v) info.push([f.name, v.length > AlleleInfoChars ? `${v.slice(0, AlleleInfoChars - 1)}…` : v]);
  }
  const id = text('label');
  const qual = text('qual');
  // VCF's `.` is a missing FILTER (no row)
  const filtered = text('filter');
  const filter = filtered === '.' ? '' : filtered;
  const rows: [string, string][] = [];
  if (id) rows.push(['ID', id]);
  if (qual) rows.push(['QUAL', qual]);
  if (filter) rows.push(['FILTER', filter]);
  rows.push(...info.slice(0, AlleleInfoRows));
  if (info.length > AlleleInfoRows) rows.push(['INFO', `+${info.length - AlleleInfoRows} more`]);
  return { id, category: text('category'), qual, filter, info, rows };
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
