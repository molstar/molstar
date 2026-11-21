/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, type Issue, issueOf } from '@molstar/core/repr/tree';

// The GFF3 and GTF codes (_spec/README.md "Data", GFF3 and GTF): one issue a cause, with its count and its first ten
// places (EC-86), its severity the catalog's (EC-87). GFF3's under `seq.source.gff`, GTF's own under `seq.source.gtf`.

export const GffIssues = defineIssues('seq.source.gff', {
  'seq.source.gff.line': {
    severity: 'warning',
    description:
      'A line that does not parse (not nine columns, a position not a whole number from 1, an end before its start): left out',
  },
  'seq.source.gff.type': {
    severity: 'info',
    description: 'A type outside `transcriptTypes` with exon or CDS children read as a transcript',
    hint: 'Add the type to the runtime’s `transcriptTypes`',
  },
  'seq.source.gff.parent': {
    severity: 'warning',
    description:
      'A `Parent` chain deeper than 16, a cycle, or a parent no line has: those rows kept as items of their own',
  },
  'seq.source.gff.attributes': {
    severity: 'warning',
    description: 'A row with more than 256 attributes: the rest dropped',
  },
  'seq.source.gff.transform': {
    severity: 'info',
    description: 'The format fixes its positions (from 1, ends inclusive): `base` and `ends` not read',
  },
  'seq.source.gff.region': {
    severity: 'info',
    description: 'A `region`, `chromosome` or `supercontig` feature over a whole sequence describes it: not an item',
  },
  'seq.source.gff.fasta': {
    severity: 'info',
    description:
      'Everything after `##FASTA` is ignored (its bytes counted against the bound): letters are a FASTA source’s',
  },
  'seq.source.gff.index': {
    severity: 'info',
    description: 'GFF3 and GTF are read whole here: `index` not read',
  },
  'seq.source.gff.extent': {
    severity: 'info',
    description: 'A child outside its parent’s extent: the parent’s extent kept (malformed GFF3)',
  },
  'seq.source.gff.encoding': {
    severity: 'info',
    description: 'Replacement characters (U+FFFD): the text may not be UTF-8',
  },
});

export const GtfIssues = defineIssues('seq.source.gtf', {
  'seq.source.gtf.line': {
    severity: 'warning',
    description:
      'A line that does not parse (not nine columns, a position not a whole number from 1, an end before its start, attributes that are not `key "value";`): left out',
  },
  'seq.source.gtf.extent': {
    severity: 'warning',
    description: 'No transcript lines: a transcript’s extent is the hull of its rows, exact only when read whole',
    hint: 'Read it whole, or extents are partial by region',
  },
});

export type GffCode = keyof typeof GffIssues.codes | keyof typeof GtfIssues.codes;

/** Places an issue keeps. */
const MaxPlaces = 10;

/** The causes a read met, one issue each with its count and first places ("line 12"). */
export class GffCauses {
  private readonly byKey = new Map<string, { code: GffCode; message: string; count: number; places: string[] }>();

  add(code: GffCode, message: string, place?: string) {
    const key = `${code}|${message}`;
    let c = this.byKey.get(key);
    if (!c) {
      c = { code, message, count: 0, places: [] };
      this.byKey.set(key, c);
    }
    c.count++;
    if (place !== undefined && c.places.length < MaxPlaces) c.places.push(place);
  }

  /** Whether a cause with `code` was said. */
  has(code: GffCode): boolean {
    for (const c of this.byKey.values()) if (c.code === code) return true;
    return false;
  }

  issues(address: string): Issue[] {
    return [...this.byKey.values()].map((c) => {
      const where = c.places.length ? ` (${c.places.join(', ')}${c.count > c.places.length ? ', …' : ''})` : '';
      const message = `${c.message}: ${c.count}${where}`;
      if (c.code.startsWith('seq.source.gtf.')) {
        const code = c.code as keyof typeof GtfIssues.codes;
        const hint = (GtfIssues.codes[code] as { hint?: string }).hint;
        return issueOf(GtfIssues, code, message, { address, ...(hint && { hint }) });
      }
      const code = c.code as keyof typeof GffIssues.codes;
      const hint = (GffIssues.codes[code] as { hint?: string }).hint;
      return issueOf(GffIssues, code, message, { address, ...(hint && { hint }) });
    });
  }
}
