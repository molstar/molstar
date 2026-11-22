/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, issueOf } from '@molstar/core/repr/tree';
import { Causes } from './internal/causes.js';

// The GenBank codes (_spec/README.md "Data", GenBank): one issue a cause, with its count and its first ten places
// (EC-86), its severity the catalog's (EC-87). A feature is refused alone, never the file, but for the file's bounds.

export const GenbankIssues = defineIssues('seq.source.genbank', {
  'seq.source.genbank.location': {
    severity: 'warning',
    description:
      'A feature whose location is not read (a remote reference, a site between two bases `a^b`, a position past the record’s length, nested deeper than 16, more pieces than the bound, malformed): left out',
  },
  'seq.source.genbank.feature': {
    severity: 'warning',
    description:
      'A feature past a bound (more qualifiers than `maxQualifiers`, a qualifier’s text past `maxQualifierBytes`) or a feature table line that does not parse: left out',
  },
  'seq.source.genbank.approximate': {
    severity: 'info',
    description: 'A location `a.b` (one base unknown within): read as its hull',
  },
  'seq.source.genbank.translation': {
    severity: 'info',
    description: '`/translation` qualifiers are not read (the protein is the CDS’s letters translated)',
  },
  'seq.source.genbank.header': {
    severity: 'info',
    description: 'Header fields not read (REFERENCE, COMMENT, DBLINK, KEYWORDS, …)',
  },
  'seq.source.genbank.length': {
    severity: 'warning',
    description: 'ORIGIN’s letters disagree with LOCUS’s length: the letters’ length is the sequence’s',
  },
  'seq.source.genbank.record': {
    severity: 'warning',
    description:
      'A record not read: a CONTIG or WGS assembly of other records, a LOCUS line without a name or length, a name another record has, or the last record cut short (no `//`)',
  },
  'seq.source.genbank.names': {
    severity: 'warning',
    description: 'More than 256 qualifier names in the file: the rest not fields',
  },
  'seq.source.genbank.transform': {
    severity: 'info',
    description: 'The format fixes its positions (from 1, ends inclusive): `base` and `ends` not read',
  },
  'seq.source.genbank.encoding': {
    severity: 'info',
    description: 'Replacement characters (U+FFFD): the text may not be UTF-8',
  },
});

export type GenbankCode = keyof typeof GenbankIssues.codes;

/** The causes a read met, one issue each with its count and first places ("line 12"). */
export class GenbankCauses extends Causes<GenbankCode> {
  constructor() {
    super((code, message, address) => issueOf(GenbankIssues, code, message, { address }));
  }
}
