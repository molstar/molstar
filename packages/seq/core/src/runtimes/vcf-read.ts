/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { TabixLimits } from '@molstar/io/common/tabix';
import { SourceIssues } from '../data/issues.js';
import type { Coverage, OpenContext, OpenedSource, ShapeInfo, SourceEntry, SystemInfo } from '../data/protocol.js';
import { type GenotypeSchema, genotypeSchema } from './vcf-samples.js';
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
  /** The coarsest level listed by region (4,096: 4 Mb tiles); a level is listed when each of its tiles fits. */
  readonly maxLevel: number;
  /** Bins a tile (the store's `tileBins`, 1,024): a level's tiles are `level × tileBins` positions. */
  readonly tileBins: number;
  /** A matrix answer's bytes, counted before it is allocated (region-reads' `maxTileBytes`, 16 MiB). */
  readonly maxMatrixBytes: number;
  /** io's other tabix bounds of a read by region (`DefaultTabixLimits` for those not given). */
  readonly tabix?: Partial<TabixLimits>;
}

export const DefaultVcfOptions: VcfOptions = {
  maxLineBytes: 16 * 1024 * 1024,
  maxHeaderBytes: 16 * 1024 * 1024,
  maxInfoFields: 256,
  maxSamples: 100_000,
  maxAlleles: 256,
  maxLevel: 4096,
  tileBins: 1024,
  maxMatrixBytes: 16 * 1024 * 1024,
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

/**
 * What a header gives a source: the schema of its INFO keys, and its systems' lengths (`##contig`, EC-188), a length
 * past 2^31 − 1 refused by name (the system not read, the others are).
 */
export function headerSchema(o: VcfOpened, header: VcfHeader) {
  const { source, options, causes } = o;
  // one cause, counted once a line (its message then says the count once)
  for (let k = 0; k < header.unread; k++)
    causes.add('seq.source.vcf.header', 'Header lines not read (no ID, or an unknown Type)');
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

/**
 * A coverage of `names` (their lengths when known), each serving `ranges` with the schema's fields, and a `matrix` of
 * the samples held (its states) when there are any.
 */
export function coverageOf(
  names: readonly string[],
  lengths: ReadonlyMap<string, number | undefined>,
  schema: VcfSchema,
  levels: readonly number[],
  genotypes?: GenotypeSchema,
): Coverage {
  const fields = fieldInfos(schema);
  const systems: SystemInfo[] = names.map((name) => {
    const length = lengths.get(name);
    return length !== undefined ? { name, length } : { name };
  });
  const shapes: ShapeInfo[] = [];
  for (const system of names) {
    shapes.push({ system, shape: 'ranges', fields });
    if (genotypes) shapes.push({ system, shape: 'matrix', states: genotypes.states, fields: genotypes.fields });
  }
  return { systems, shapes, levels };
}

/** The samples a source holds (none: no matrix), said at open: names the file lacks. */
export function genotypesOf(o: VcfOpened, header: VcfHeader): GenotypeSchema | undefined {
  const samples = (o.source as { readonly samples?: readonly string[] }).samples;
  const g = genotypeSchema(header, samples, o.causes);
  return g.states.length ? g : undefined;
}
