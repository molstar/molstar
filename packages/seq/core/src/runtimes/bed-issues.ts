/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, type Issue, issueOf } from '@molstar/core/repr/tree';

// The BED family's codes (_spec/README.md "Data", BED family): one issue a cause, with its count and its first ten
// places (EC-86), its severity the catalog's (EC-87).

export const BedIssues = defineIssues('seq.source.bed', {
  'seq.source.bed.line': {
    severity: 'warning',
    description:
      'A line that does not parse (a position not a whole number, an end before its start, blocks outside the item): left out',
  },
  'seq.source.bed.cells': {
    severity: 'info',
    description: 'A cell its field cannot hold read as missing, or columns past the first line’s not read',
  },
  'seq.source.bed.field': {
    severity: 'error',
    description: 'A field of `fields` names no column: the source is unavailable',
    hint: "Name a column by its 0-based number (`{ stain: '4' }`) or by its name in the format",
  },
  'seq.source.bed.transform': {
    severity: 'info',
    description: 'The format fixes its positions (from 0, ends exclusive): `base` and `ends` not read',
  },
  'seq.source.bed.whitespace': {
    severity: 'info',
    description: 'Lines without a tab are split on spaces',
  },
  'seq.source.bed.index': {
    severity: 'info',
    description: 'BEDPE and chrom.sizes are read whole: `index` not read',
  },
  'seq.source.bed.encoding': {
    severity: 'info',
    description: 'Replacement characters (U+FFFD): the text may not be UTF-8',
  },
});

export type BedCode = keyof typeof BedIssues.codes;

/** Places an issue keeps. */
const MaxPlaces = 10;

/** The causes a read met, one issue each with its count and first places ("line 12"). */
export class BedCauses {
  private readonly byKey = new Map<string, { code: BedCode; message: string; count: number; places: string[] }>();

  add(code: BedCode, message: string, place?: string) {
    const key = `${code}|${message}`;
    let c = this.byKey.get(key);
    if (!c) {
      c = { code, message, count: 0, places: [] };
      this.byKey.set(key, c);
    }
    c.count++;
    if (place !== undefined && c.places.length < MaxPlaces) c.places.push(place);
  }

  get size() {
    return this.byKey.size;
  }

  issues(address: string): Issue[] {
    return [...this.byKey.values()].map((c) => {
      const where = c.places.length ? ` (${c.places.join(', ')}${c.count > c.places.length ? ', …' : ''})` : '';
      return issueOf(BedIssues, c.code, `${c.message}: ${c.count}${where}`, {
        address,
        ...(BedIssues.codes[c.code].hint && { hint: BedIssues.codes[c.code].hint }),
      });
    });
  }
}
