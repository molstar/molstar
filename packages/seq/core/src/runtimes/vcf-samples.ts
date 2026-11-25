/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { FieldInfo } from '../data/protocol.js';
import type { VcfHeader } from './vcf-header.js';
import type { VcfCauses } from './vcf-issues.js';

// A VCF's samples held (_spec/README.md "Data", VCF: the matrix): what a source with samples lists at open (the matrix's
// states and value channels), apart from the matrix's reader (`vcf-genotypes.ts`, loaded on the first matrix request).

/** The samples a source holds and the FORMAT keys a matrix serves as values. */
export interface GenotypeSchema {
  /** The states: the samples held, in the header's order. */
  readonly states: readonly string[];
  /** Each state's sample column (0 the first sample), ascending. */
  readonly columns: Int32Array;
  /** The value channels: `GT` (copies), then the numeric FORMAT keys of one value a sample or an allele. */
  readonly fields: readonly FieldInfo[];
  /** A value channel's `Number` (`1` or `A`). */
  readonly number: ReadonlyMap<string, string>;
}

/**
 * The samples held: the header's, narrowed by `samples` (names the header lacks said once, as is a narrowing of a file
 * without samples); the matrix's value channels from the header's `##FORMAT` lines.
 */
export function genotypeSchema(
  header: VcfHeader,
  samples: readonly string[] | undefined,
  causes: VcfCauses,
): GenotypeSchema {
  let columns: number[];
  if (samples === undefined) columns = header.samples.map((_, k) => k);
  else {
    const wanted = new Set(samples);
    columns = [];
    for (let k = 0; k < header.samples.length; k++) if (wanted.has(header.samples[k])) columns.push(k);
    const have = new Set(header.samples);
    const missing = samples.filter((s) => !have.has(s));
    if (!header.samples.length) causes.add('seq.source.vcf.samples', 'The file has no samples: `samples` reads none');
    else if (missing.length)
      causes.add(
        'seq.source.vcf.samples',
        `Samples the file has not: ${missing.length} (${missing.slice(0, 10).join(', ')})`,
      );
  }
  const fields: FieldInfo[] = [{ name: 'GT', type: 'float32', doc: 'Copies of the allele (NaN uncalled)' }];
  const number = new Map<string, string>([['GT', '1']]);
  for (const d of header.format) {
    if (d.id === 'GT' || (d.number !== '1' && d.number !== 'A')) continue;
    if (d.type !== 'Integer' && d.type !== 'Float') continue;
    fields.push({ name: d.id, type: 'float32' });
    number.set(d.id, d.number);
  }
  return { states: columns.map((k) => header.samples[k]), columns: Int32Array.from(columns), fields, number };
}

/** What a matrix read may allocate: its tile's bytes (the store's `maxTileBytes`, 16 MiB). */
export interface MatrixBounds {
  readonly maxBytes: number;
}
