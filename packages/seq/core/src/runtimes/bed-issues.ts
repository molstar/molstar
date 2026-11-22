/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, issueOf } from '@molstar/core/repr/tree';
import { Causes } from './internal/causes.js';

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

/** The causes a read met, one issue each with its count and first places ("line 12"). */
export class BedCauses extends Causes<BedCode> {
  constructor() {
    super((code, message, address) =>
      issueOf(BedIssues, code, message, {
        address,
        ...(BedIssues.codes[code].hint && { hint: BedIssues.codes[code].hint }),
      }),
    );
  }
}
