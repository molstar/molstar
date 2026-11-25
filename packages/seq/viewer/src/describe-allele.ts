/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import { type RangesList, VcfFixedFieldNames } from '@molstar/seqx-core/data';

// A variant allele as its card says it (_spec/README.md "Behaviors", hoverCard, A variant allele's card): the item
// card's words for a source row holding `ref` and `alt` (a VCF's).

/** The INFO fields an allele's card lists before "+k more" (the rest in the source's rows). */
export const AlleleInfoRows = 8;

/** The longest INFO value a card row shows (a long list of ids, a description): the rest as an ellipsis. */
export const AlleleInfoChars = 40;

/** A source's fixed fields (core's `VcfFixedFieldNames`: core _spec "Data", VCF, Fields): what is not an INFO key. */
const FixedFields = new Set<string>(VcfFixedFieldNames);

/** A variant allele's card words: its record's ID, category, QUAL and FILTER, and its INFO fields. */
export interface Allele {
  readonly id: string;
  readonly category: string;
  readonly qual: string;
  readonly filter: string;
  readonly info: readonly (readonly [string, string])[];
  /** Its card's rows after the category: ID, QUAL, FILTER, the first `AlleleInfoRows` INFO fields and "+k more". */
  readonly rows: readonly (readonly [string, string])[];
  /** Its REF and ALT as the source gives them ('' when not). */
  readonly ref: string;
  readonly alt: string;
}

/** The longest allele a title spells out: past it, its first letters and its length (the id holds a hash). */
export const AlleleTitleChars = 32;

/** An allele's text in a title: itself, or past `AlleleTitleChars` its first letters and its length (`GAGGCT… 469 bp`). */
const shortAllele = (a: string) => (a.length <= AlleleTitleChars ? a : `${a.slice(0, 6)}… ${a.length} bp`);

/**
 * An allele's card title (_spec/README.md "Behaviors", hoverCard, A variant allele's card): its id (`7674964:A>G`),
 * with a long REF or ALT, hashed in the id (`7674964:469:edaeae47>G`), as its first letters and its length
 * (`7674964:GAGGCT… 469 bp>G`): a reader reads a length, not a hash.
 */
export function alleleTitle(key: string | number, allele: Allele): string {
  const id = String(key);
  const m = /^(\d+):/.exec(id);
  if (!m || (allele.ref.length <= AlleleTitleChars && allele.alt.length <= AlleleTitleChars)) return id;
  return `${m[1]}:${shortAllele(allele.ref)}>${shortAllele(allele.alt)}`;
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
  return { id, category: text('category'), qual, filter, info, rows, ref: text('ref'), alt: text('alt') };
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

/**
 * What a source's row says of a transcript (core _spec "Data", GFF3 and GTF, Tables): its `canonical` field, `Canonical`
 * or `MANE Select` when the file's tags say so; '' when the row has none (the viewer never guesses one).
 */
export function canonicalOf(t: RangesList, i: number): string {
  const f = t.fields?.find((x) => x.name === 'canonical');
  return f ? cellText(f.column, i) : '';
}
