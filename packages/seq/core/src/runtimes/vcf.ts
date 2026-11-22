/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { abortError, type FormatRuntime, type OpenContext, type SourceEntry } from '../data/protocol.js';
import { wholeLevel } from '../data/seq-format.js';
import { failure, IndexHint, refuse } from './internal/answers.js';
import { textBytes } from './internal/text.js';
import { wholeLines } from './internal/whole-lines.js';
import { readVcfHeader } from './vcf-header.js';
import { VcfCauses } from './vcf-issues.js';
import { vcfAnswer } from './vcf-answer.js';
import { parseLine, type VcfRows, VcfRowsBuilder } from './vcf-parse.js';
import {
  coverageOf,
  DefaultVcfOptions,
  genotypesOf,
  headerSchema,
  refusedIssues,
  type VcfOpened,
  type VcfOptions,
  type VcfSource,
} from './vcf-read.js';

// VCF sites (_spec/README.md "Data", VCF), the opt-in runtime `@molstar/seqx-core/runtimes/vcf`: a row an ALT allele
// (EC-196), every allele of a line with the line's extent by the index's rule, the anchor base kept; the fixed fields
// and the header's INFO keys typed by it. Read whole (plain, gzip or BGZF by its magic, within `limits`), or by region
// at the levels its index's bytes allow through a tabix `index` (`vcf-indexed.ts`, loaded then; EC-186), by one line
// parser either way. Positions
// as the format writes them (base 1, converted at the edge): the source's `base` and `ends` are not read (EC-195).
// Genotypes are read only for a `matrix` request, for the samples held (`vcf-genotypes.ts`).

export type { VcfContig, VcfDeclaration, VcfHeader, VcfHeaderLimits, VcfHeaderResult, VcfType } from './vcf-header.js';
export { DefaultVcfHeaderLimits, readVcfHeader } from './vcf-header.js';
export { VcfIssues } from './vcf-issues.js';
export { MissingInt, VcfCategories } from './vcf-parse.js';
export { DefaultVcfOptions, type VcfOptions, type VcfSource } from './vcf-read.js';

export const VcfFormats = ['vcf'] as const;

/** The bytes a whole read parses: inline `value` or the file at `uri`, inflated by its magic within `maxBytes`. */
async function wholeBytes(o: VcfOpened): Promise<Uint8Array> {
  const { source, ctx } = o;
  const at = source.address;
  let bytes: Uint8Array;
  if (typeof source.value === 'string') bytes = new TextEncoder().encode(source.value);
  else if (source.uri !== undefined) {
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch {
      if (ctx.signal.aborted) throw abortError();
      refuse(at, 'seq.source.unreadable', `The file ${source.uri} could not be read`);
    }
  } else refuse(at, 'seq.source.no_data', 'A vcf source needs a uri or a text value');
  const text = await textBytes(bytes, ctx.limits.maxBytes, ctx.signal);
  if (!text.ok)
    refuse(at, text.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', text.message, IndexHint);
  if (ctx.signal.aborted) throw abortError();
  return text.bytes;
}

async function openWhole(o: VcfOpened): Promise<VcfSource> {
  const { source, ctx, options, causes } = o;
  const at = source.address;
  const bytes = await wholeBytes(o);
  const lines = wholeLines(bytes, { meta: '#', maxLines: ctx.limits.maxRows, maxLineBytes: options.maxLineBytes });
  if (!lines.ok)
    refuse(at, lines.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', lines.message, IndexHint);
  const { header: h, lines: l } = lines.value;
  const headerEnd = h.lineEnd.length ? h.lineEnd[h.lineEnd.length - 1] + 1 : 0;
  const read = readVcfHeader(bytes.subarray(0, headerEnd), options);
  if (!read.ok) refuse(at, 'seq.source.malformed', read.message);
  const { schema, lengths, refused } = headerSchema(o, read.header);
  const refusedSet = new Set(refused);

  const builders = new Map<string, VcfRowsBuilder>();
  const builder = (system: string) => {
    let b = builders.get(system);
    if (!b && !refusedSet.has(system)) builders.set(system, (b = new VcfRowsBuilder(schema, lengths.get(system))));
    return b;
  };
  let made = 0;
  const maxRows = ctx.limits.maxRows;
  const encoding = { replaced: false };
  const lineOptions = { maxAlleles: options.maxAlleles, builder, rowsLeft: () => maxRows - made };
  for (let i = 0; i < l.lineStart.length; i++) {
    if (l.bytes[l.lineStart[i]] === 35) continue;
    const b = l.lineStart[i];
    const got = parseLine(l.bytes, b, l.lineEnd[i], `line ${l.lineNumber[i]}`, lineOptions, causes, encoding);
    if (got < 0) refuse(at, 'seq.source.too_large', `More than ${maxRows} alleles`, IndexHint);
    made += got;
  }
  const rows = new Map<string, VcfRows>();
  for (const [name, b] of builders) rows.set(name, b.finish(bytes));
  const names = [...lengths.keys()];
  for (const name of rows.keys()) if (!lengths.has(name)) names.push(name);
  let extent = 1;
  for (const name of names) extent = Math.max(extent, lengths.get(name) ?? 0, rows.get(name)?.maxEnd ?? 0);
  // read whole, sites are served at level 1 and at the level whose one tile holds the longest system (the inline rule)
  const whole = wholeLevel(extent);
  const header = read.header;
  const genotypes = genotypesOf(o, header);
  const coverage = coverageOf(names, lengths, schema, whole > 1 ? [1, whole] : [1], genotypes);
  const answer = { maxRows, genotypes, matrix: { maxBytes: options.maxMatrixBytes }, density: false };
  const served = new Set(names);
  return {
    coverage,
    issues: [...causes.issues(at), ...refusedIssues(at, refused)],
    header,
    async read(r, signal) {
      if (signal.aborted) throw abortError();
      if (!served.has(r.system)) return failure('not-served', `No sequence ${r.system}`);
      if (!coverage.levels.includes(r.level)) return failure('not-served', `No level ${r.level}`);
      const a = await vcfAnswer(rows.get(r.system), schema, r, answer);
      if (signal.aborted) throw abortError();
      return a;
    },
    dispose() {},
  };
}

async function open(source: SourceEntry, ctx: OpenContext, options: VcfOptions): Promise<VcfSource> {
  if (source.format !== 'vcf') refuse(source.address, 'seq.source.no_data', `Not VCF: ${source.format}`);
  const causes = new VcfCauses();
  if (source.base !== undefined || source.ends !== undefined)
    causes.add('seq.source.vcf.transform', 'VCF fixes base 1 and the index’s extents');
  const o: VcfOpened = { source, ctx, options, causes };
  if (source.index !== undefined && source.uri !== undefined)
    return (await import('./vcf-indexed.js')).openIndexed(o, source.uri, source.index);
  return openWhole(o);
}

/** A VCF runtime with its own bounds (`DefaultVcfOptions` for those not given). */
export function vcfRuntime(over: Partial<VcfOptions> = {}): FormatRuntime {
  const options = { ...DefaultVcfOptions, ...over };
  return { name: 'vcf', formats: VcfFormats, open: (source, ctx) => open(source, ctx, options) };
}

/** The runtime of `vcf` (sites) with the default bounds: opt in by listing it. */
export const VcfRuntime: FormatRuntime = vcfRuntime();
