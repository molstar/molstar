/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues, issueOf } from '@molstar/core/repr/tree';
import { Causes } from './internal/causes.js';

// The VCF codes (_spec/README.md "Data", VCF): one issue a cause, with its count and its first ten places (EC-86), its
// severity the catalog's (EC-87).

export const VcfIssues = defineIssues('seq.source.vcf', {
  'seq.source.vcf.line': {
    severity: 'warning',
    description:
      'A line that does not parse (fewer than eight columns, a position not a whole number to 2^31 − 1, no REF): left out',
  },
  'seq.source.vcf.alleles': {
    severity: 'warning',
    description: 'A record with more ALT alleles than `maxAlleles`: left out',
    hint: 'Split the record (`bcftools norm -m-`), or raise the runtime’s `maxAlleles`',
  },
  'seq.source.vcf.allele': {
    severity: 'info',
    description: 'An ALT allele that is no variant (`.`, `*`, `<NON_REF>`, `<*>`): no item',
  },
  'seq.source.vcf.cell': {
    severity: 'info',
    description:
      'A cell its field cannot hold read as missing, or an INFO key named as a fixed field (`ref`, `qual` …) not read',
  },
  'seq.source.vcf.info': {
    severity: 'info',
    description: 'INFO keys past `maxInfoFields`, or keys the header does not declare: not fields',
    hint: 'Name the keys wanted in the source’s `fields` or `types`: they are read first',
  },
  'seq.source.vcf.header': {
    severity: 'info',
    description: 'A `##contig`, `##INFO` or `##FORMAT` line without an `ID` or with an unknown `Type`: not read',
  },
  'seq.source.vcf.range': {
    severity: 'warning',
    description: 'A record starts past its contig’s length (left out), or runs past it (clipped); both counted here',
  },
  'seq.source.vcf.transform': {
    severity: 'info',
    description:
      'The format fixes its positions (from 1, `POS − 1` to the end the index reads): `base` and `ends` not read',
  },
  'seq.source.vcf.samples': {
    severity: 'info',
    description: '`samples` names samples the file has not, or the file has none: those are not read',
  },
  'seq.source.vcf.encoding': {
    severity: 'info',
    description: 'Replacement characters (U+FFFD): the text may not be UTF-8',
  },
});

export type VcfCode = keyof typeof VcfIssues.codes;

/** The causes a read met, one issue each with its count and first places ("line 12", "17:7675088"). */
export class VcfCauses extends Causes<VcfCode> {
  constructor() {
    super((code, message, address) => {
      const hint = (VcfIssues.codes[code] as { hint?: string }).hint;
      return issueOf(VcfIssues, code, message, { address, ...(hint && { hint }) });
    });
  }
}
