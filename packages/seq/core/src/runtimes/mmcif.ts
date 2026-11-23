/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { CifFrame } from '@molstar/io/reader/cif/data-model';
import { SourceIssues } from '../data/issues.js';
import {
  abortError,
  type Coverage,
  type FormatRuntime,
  isAbortError,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
} from '../data/protocol.js';
import { SourceOpenError } from '../data/seq-runtime.js';
import { DefaultStructureLimits, type StructureLimits, structureNodes } from '../io/structure.js';
import type { MmcifReadOptions } from './mmcif-description.js';
import { MmcifIssues } from './mmcif-issues.js';

// mmCIF (_spec/README.md "Data", mmCIF), the opt-in runtime `@molstar/seqx-core/runtimes/mmcif`:
// a source of format `mmcif` read whole (gzipped or not) by `@molstar/io`'s text CIF parser, its categories read into
// a `StructureDescription` (`mmcif-description.ts`; both loaded when a source first opens, a chunk a viewer that lists
// the runtime and opens none never pays) and that into nodes by `structureNodes`, the one converter every adapter of a
// structure shares (_spec/README.md "Seq io"). It serves the nodes up front and, by region, the units'
// per-residue values by `field`: `b_factor`, `occupancy` and, in a model with them, `plddt` (the default when present).

export { MmcifIssues } from './mmcif-issues.js';

/** The source params (_spec/README.md "Data", mmCIF) a runtime is made with; a source's own `assembly` wins. */
export interface MmcifOptions {
  /** The assembly read; absent: the deposited model. */
  readonly assembly?: string;
  /** Waters read: `false`. */
  readonly water?: boolean;
  /** A molecule's name: `description` (the entity's own words), or `composition` for a branched one. */
  readonly name?: 'description' | 'composition';
  /** The bounds the description is read within, over `structureNodes`'s defaults. */
  readonly limits?: Partial<StructureLimits>;
}

/** The value fields a source may serve, the first one it has its default. */
const Fields = ['plddt', 'b_factor', 'occupancy'] as const;

/** Whether bytes are gzipped (RFC 1952's magic). */
const isGzip = (b: Uint8Array) => b.length >= 2 && b[0] === 0x1f && b[1] === 0x8b;

/** Inflates gzipped bytes, at most `limit` of them (else too large), chunks gathered rather than sized up front. */
async function gunzip(bytes: Uint8Array, limit: number): Promise<Uint8Array> {
  const stream = new Blob([bytes as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new DecompressionStream('gzip'));
  const reader = stream.getReader();
  const chunks: Uint8Array[] = [];
  let size = 0;
  for (;;) {
    const next = await reader.read();
    if (next.done) break;
    size += next.value.length;
    if (size > limit) {
      await reader.cancel().catch(() => {});
      throw new TooLarge(`inflates past ${limit} bytes`);
    }
    chunks.push(next.value);
  }
  const out = new Uint8Array(size);
  let o = 0;
  for (const c of chunks) {
    out.set(c, o);
    o += c.length;
  }
  return out;
}

class TooLarge extends Error {}

/** The first data block of a text CIF, by `@molstar/io`'s parser, loaded on first use. */
async function parseFrame(bytes: Uint8Array): Promise<CifFrame> {
  const { parseCifText } = await import('@molstar/io/reader/cif/text/parser');
  const parsed = await parseCifText(bytes).run();
  if (parsed.isError) throw new Error(parsed.message);
  const block = parsed.result.blocks[0];
  if (!block) throw new Error('no data block');
  return block;
}

/**
 * Opens mmCIF bytes (gzipped or not) as a source: its nodes up front, its units' values by region. Rejects with a
 * `SourceOpenError` when the bytes are past `limits.maxBytes`, do not parse, or lack the assembly asked.
 */
export async function openMmcif(
  bytes: Uint8Array,
  ctx: Pick<OpenContext, 'signal' | 'limits'>,
  source: SourceEntry = { address: '', format: 'mmcif' },
  options: MmcifOptions = {},
): Promise<OpenedSource> {
  const address = source.address;
  const fail = (code: 'seq.source.malformed' | 'seq.source.unreadable', message: string): never => {
    throw new SourceOpenError([issueOf(SourceIssues, code, message, { address })]);
  };
  let frame: CifFrame;
  try {
    if (bytes.length > ctx.limits.maxBytes) throw new TooLarge(`larger than ${ctx.limits.maxBytes} bytes`);
    const text = isGzip(bytes) ? await gunzip(bytes, ctx.limits.maxBytes) : bytes;
    if (ctx.signal.aborted) throw abortError();
    frame = await parseFrame(text);
  } catch (e) {
    if (ctx.signal.aborted || isAbortError(e)) throw abortError();
    const why = String((e as Error)?.message ?? e);
    return fail(
      e instanceof TooLarge ? 'seq.source.unreadable' : 'seq.source.malformed',
      `The mmCIF source cannot be read: ${why}`,
    );
  }
  if (ctx.signal.aborted) throw abortError();
  // the source's `assembly` param rides on its entry (resolve's `SourceKeys`), not in `SourceEntry`'s type
  const own = (source as { readonly assembly?: unknown }).assembly;
  const assembly = typeof own === 'string' ? own : options.assembly;
  const readOptions: MmcifReadOptions = {
    ...(assembly !== undefined ? { assembly } : {}),
    water: !!options.water,
    name: options.name ?? 'description',
    maxUnits: options.limits?.maxUnits ?? DefaultStructureLimits.maxUnits,
  };
  // the reading of categories is loaded with the parser, on first open (Size)
  const { MmcifAssemblyError, mmcifDescription } = await import('./mmcif-description.js');
  let described: ReturnType<typeof mmcifDescription>;
  try {
    described = mmcifDescription(frame, readOptions, address);
  } catch (e) {
    if (e instanceof MmcifAssemblyError)
      throw new SourceOpenError([issueOf(MmcifIssues, 'seq.io.mmcif.assembly', e.message, { address })]);
    // a file is untrusted input: whatever else its reading met is said, never rethrown
    return fail('seq.source.malformed', `The mmCIF source cannot be read: ${String((e as Error)?.message ?? e)}`);
  }
  const nodes = structureNodes(described.description, {
    source: address,
    water: readOptions.water,
    ...(options.limits ? { limits: options.limits } : {}),
  });
  const fields = Fields.filter((f) => nodes.values(f));
  const lengths = new Map(nodes.systems.map((s) => [s.name, s.length ?? 0]));
  const coverage: Coverage = {
    systems: nodes.systems,
    shapes: fields.length
      ? nodes.systems.map((s) => ({
          system: s.name,
          shape: 'values' as const,
          fields: fields.map((name) => ({ name, type: 'float32' as const })),
        }))
      : [],
    levels: [1],
    nodes: nodes.nodes,
  };
  const issues: Issue[] = [...described.issues, ...nodes.issues];
  const cache = new Map<string, ReadonlyMap<string, Float32Array> | undefined>();
  const valuesOf = (field: string) => {
    if (!cache.has(field)) cache.set(field, nodes.values(field));
    return cache.get(field);
  };
  const notServed = (why: string): RegionAnswer => ({
    ok: false,
    kind: 'not-served',
    issue: issueOf(SourceIssues, 'seq.source.not_served', why, { address }),
  });

  async function read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer> {
    if (signal.aborted) throw abortError();
    const length = lengths.get(request.system);
    if (length === undefined) return notServed(`The mmCIF source has no unit ${request.system}`);
    if (request.level !== 1) return notServed(`The mmCIF source serves level 1 only, not ${request.level}`);
    if (!request.shapes.includes('values')) return { ok: true, tables: {} };
    const field = request.field ?? fields[0];
    const values = field === undefined ? undefined : valuesOf(field)?.get(request.system);
    if (!values) return notServed(`The mmCIF source has no value ${request.field ?? ''} on ${request.system}`);
    const start = Math.max(0, Math.min(length, request.start));
    const end = Math.max(start, Math.min(length, request.end));
    return { ok: true, tables: { values: { form: 'positions', start, value: values.slice(start, end) } } };
  }

  return { coverage, issues, read, dispose() {} };
}

async function open(source: SourceEntry, ctx: OpenContext, options: MmcifOptions): Promise<OpenedSource> {
  const uri = source.uri;
  if (!uri)
    throw new SourceOpenError([
      issueOf(SourceIssues, 'seq.source.no_data', 'An mmCIF source needs a uri', { address: source.address }),
    ]);
  let bytes: Uint8Array;
  try {
    bytes = await ctx.assets.readBytes(uri, ctx.signal);
  } catch (e) {
    if (isAbortError(e) || ctx.signal.aborted) throw abortError();
    throw new SourceOpenError([
      issueOf(SourceIssues, 'seq.source.unreadable', `The file ${uri} could not be read`, { address: source.address }),
    ]);
  }
  return openMmcif(bytes, ctx, source, options);
}

/** An mmCIF runtime with its own source params and bounds (_spec/README.md "Data", mmCIF). */
export function mmcifRuntime(options: MmcifOptions = {}): FormatRuntime {
  return { name: 'mmcif', formats: ['mmcif'], open: (source, ctx) => open(source, ctx, options) };
}

/** The runtime of mmCIF (`format: 'mmcif'`), the deposited model without waters: opt in by listing it in a registry. */
export const MmcifRuntime: FormatRuntime = mmcifRuntime();
