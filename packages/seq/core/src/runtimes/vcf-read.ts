/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { TabixLimits } from '@molstar/io/common/tabix';
import { SourceIssues } from '../data/issues.js';
import type { Coverage, OpenContext, OpenedSource, SourceEntry, SystemInfo } from '../data/protocol.js';
import type { VcfHeader } from './vcf-header.js';
import type { VcfCauses } from './vcf-issues.js';
import { fieldInfos, MaxPosition, type VcfSchema, vcfSchema } from './vcf-parse.js';

// What a VCF's whole and region reads share (_spec/README.md "Data", VCF): the runtime's bounds, the schema and the
// systems a header gives, and the coverage of them.

/** What the runtime may hold and read, whatever a file says (security.md: a file is untrusted input). */
export interface VcfOptions {
  /** One line's bytes (16 MiB: a line of 100,000 samples at about 160 bytes each). */
  readonly maxLineBytes: number;
  /** The header's bytes (16 MiB: gnomAD's is a few hundred kB). */
  readonly maxHeaderBytes: number;
  /** INFO keys read as fields (256: gnomAD v4 declares about 900); the source's `fields` and `types` first. */
  readonly maxInfoFields: number;
  /** Sample names on the `#CHROM` line (100,000). */
  readonly maxSamples: number;
  /** ALT alleles a record (256): past it the record is left out. */
  readonly maxAlleles: number;
  /** io's other tabix bounds of a read by region (`DefaultTabixLimits` for those not given). */
  readonly tabix?: Partial<TabixLimits>;
}

export const DefaultVcfOptions: VcfOptions = {
  maxLineBytes: 16 * 1024 * 1024,
  maxHeaderBytes: 16 * 1024 * 1024,
  maxInfoFields: 256,
  maxSamples: 100_000,
  maxAlleles: 256,
};

/** An opened VCF: a source, and its header (what a host shows of the file: its keys, contigs and samples). */
export interface VcfSource extends OpenedSource {
  readonly header: VcfHeader;
}

/** An open in progress: the source, its context, the runtime's options and the causes said. */
export interface VcfOpened {
  readonly source: SourceEntry;
  readonly ctx: OpenContext;
  readonly options: VcfOptions;
  readonly causes: VcfCauses;
}

export const IndexHint = 'Read it by region: bgzip it and give its tabix index';

/**
 * What a header gives a source: the schema of its INFO keys, and its systems' lengths (`##contig`, EC-188), a length
 * past 2^31 − 1 refused by name (the system not read, the others are).
 */
export function headerSchema(o: VcfOpened, header: VcfHeader) {
  const { source, options, causes } = o;
  if (header.unread)
    causes.add('seq.source.vcf.header', `${header.unread} header lines not read (no ID, or an unknown Type)`);
  const schema = vcfSchema(
    header,
    source.fields,
    source.types as Readonly<Record<string, 'int32' | 'float32' | 'float64' | 'category' | 'string'>> | undefined,
    options.maxInfoFields,
    causes,
  );
  const lengths = new Map<string, number | undefined>();
  const refused: string[] = [];
  for (const c of header.contigs) {
    if (c.length !== undefined && c.length > MaxPosition) refused.push(c.name);
    else lengths.set(c.name, c.length);
  }
  return { schema, lengths, refused };
}

/** The issue of contigs refused for their length. */
export function refusedIssues(address: string, refused: readonly string[]): Issue[] {
  if (!refused.length) return [];
  const message = `Longer than 2^31 − 1, not read: ${refused.slice(0, 10).join(', ')}`;
  return [issueOf(SourceIssues, 'seq.source.too_large', message, { address })];
}

/** A coverage of `names` (their lengths when known), each serving `ranges` with the schema's fields. */
export function coverageOf(
  names: readonly string[],
  lengths: ReadonlyMap<string, number | undefined>,
  schema: VcfSchema,
  levels: readonly number[],
): Coverage {
  const fields = fieldInfos(schema);
  const systems: SystemInfo[] = names.map((name) => {
    const length = lengths.get(name);
    return length !== undefined ? { name, length } : { name };
  });
  return { systems, shapes: names.map((system) => ({ system, shape: 'ranges' as const, fields })), levels };
}
