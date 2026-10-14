/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load stages: `indexDoc`, `loadNode`, `buildComplex`, `assembleSet`, `computeProvider`, and `compile`, which runs
// them in the load order: sources (an annotation's `uri` too); complexes from providers; molecules and alignments;
// units; complexes by value; tables, named rows and mappings (`align` ones once the set holds their molecules);
// providers. A missing provider or aligner makes its node unavailable (a warning), and a node whose load reported an
// error is unavailable; nodes that refer to it are left out quietly, so nothing cascades.

import { AlignerSizeError, TreeSizeError } from "@molstar/seq-core/algorithms";
import {
  type Alignment,
  type AlignmentInput,
  type AnnotationInput,
  ColumnRoles,
  type Geometry,
  type Issue,
  KeyKind,
  type KeyTarget,
  type MappingInput,
  type Molecule,
  type MoleculeInput,
  MoleculeKind,
  type NumberingInput,
  type SeqSet,
  SeqSets,
  type SetInput,
  type TreeInput,
  Trees,
  type UnitInput,
} from "@molstar/seq-core/model";
import type { SeqQLQuery } from "@molstar/seq-core/query";
import {
  type AlignerInput,
  emptyRequest,
  type FormatRequest,
  type FormatResult,
  type Loaders,
  type Manifest,
  type SourceResult,
  type StageContext,
  type StreamDataKind,
  type Streamed,
  StreamNotServed,
  type StreamRequest,
  type StreamSource,
} from "../contract.ts";
import type { DocNode, SeqViewSpecDoc } from "../doc.ts";
import { equal } from "../equal.ts";
import { normalizeNodeParams } from "../normalize.ts";
import { queryRefs } from "../query.ts";
import type { SvsRegistries } from "../registry.ts";
import { type CheckContext, isNamespaced } from "../schema/field-schema.ts";
import type { IndexedNode } from "../schema/tree-schema.ts";
import { CellKind, cellKind, cellOk } from "../tables/cells.ts";
import { SvsTreeSchema } from "../tables/svs-tree.ts";
import { renameFields, table, typeColumns } from "../tables/table.ts";
import { validate } from "../validate.ts";
import { tableColumns } from "./formats.ts";
import { at, report, svsIssue } from "./issues.ts";
import type { AlignmentSummary, RealignInput } from "./node-providers.ts";

/**
 * The runtime property (seqspec:model#custom-data; seqspec:sources#the-loader) of what a streamed source serves: a lazy
 * alignment, a streamed annotation. `systems`: the reader's coordinate systems by ref (a complex's member units, else
 * the reader's own), each its manifest system's index; `labels`: row labels as they arrive (the viewer's stream manager
 * fills it).
 */
export const StreamProp = "svs:stream";
export interface StreamedData extends Streamed {
  readonly systems: ReadonlyMap<string, number>;
  readonly labels: Map<number, string>;
}
/** A request with the defaults: data (level -1), a view's, no labels; `system` null asks for the manifest. */
export function streamRequest(source: string, system: string | null, over: Partial<StreamRequest> = {}): StreamRequest {
  return {
    source,
    system,
    start: 0,
    stop: 0,
    row_start: 0,
    row_stop: 0,
    level: -1,
    labels: false,
    priority: "view",
    ...over,
  };
}

export interface DocIndex {
  readonly nodes: ReadonlyMap<string, IndexedNode>;
  /** Every node in document order. */
  readonly order: readonly IndexedNode[];
}

/**
 * A node's normalized params and `custom`, kept by the identity of its raw `params` (`indexDoc`'s cache), and the
 * registries' `version` it was normalized at: after a registration it is normalized again, an equal result keeping
 * its objects.
 */
export interface NormalizedParams {
  readonly kind: string;
  readonly custom: unknown;
  readonly version: number;
  readonly out: { readonly params: Record<string, unknown>; readonly custom: Record<string, unknown> | undefined };
}
/** The cache key of a node given no params, per kind. */
const noParams = new Map<string, object>();

function normalized(node: DocNode, reg: SvsRegistries, cache?: WeakMap<object, NormalizedParams>) {
  if (!cache) return normalizeNodeParams(node, SvsTreeSchema, reg);
  const key = node.params ?? noParams.get(node.kind) ?? (noParams.set(node.kind, {}).get(node.kind) as object);
  const hit = cache.get(key);
  const same = hit?.kind === node.kind;
  if (same && hit.version === reg.version && hit.custom === node.custom) return hit.out;
  const out = normalizeNodeParams(node, SvsTreeSchema, reg);
  // Normalized again (a changed `custom`, or any registration): what came out equal keeps its object.
  const keep = <T>(old: T | undefined, now: T): T => (same && old !== undefined && equal(old, now) ? old : now);
  const kept = { params: keep(hit?.out.params, out.params), custom: keep(hit?.out.custom, out.custom) };
  const entry = { kind: node.kind, custom: node.custom, version: reg.version, out: kept };
  cache.set(key, entry);
  return entry.out;
}

/**
 * Ref → node, parent, kind, path (implied refs filled). An unknown kind (a newer minor's) is left out with its subtree.
 * With `reg`, each node's params and `custom` are normalized here (`node` then carries them), so a raw document needs
 * no `normalize`; with `cache` (the caller's), a node whose raw `params` object is unchanged keeps its normalized one,
 * also across a registration that leaves it equal.
 */
export function indexDoc(
  doc: SeqViewSpecDoc,
  reg?: SvsRegistries,
  cache?: WeakMap<object, NormalizedParams>,
): DocIndex {
  const nodes = new Map<string, IndexedNode>();
  const order: IndexedNode[] = [];
  const walk = (raw: DocNode, parent: string | null, path: string) => {
    if (!Object.hasOwn(SvsTreeSchema.nodes, raw.kind)) return;
    const norm = reg ? normalized(raw, reg, cache) : null;
    const node = norm && (norm.params !== raw.params || norm.custom !== raw.custom) ? { ...raw, ...norm } : raw;
    const n = { node, ref: raw.ref ?? raw.kind, kind: raw.kind, params: node.params ?? {}, parent, path };
    nodes.set(n.ref, n);
    order.push(n);
    raw.children?.forEach((c, i) => {
      walk(c, n.ref, `${path}/children/${i}`);
    });
  };
  walk(doc.root, null, "/root");
  return { nodes, order };
}

type P = Readonly<Record<string, unknown>>;
export type SetPart = Partial<SetInput>;
/** A load stage's result: its part of the set, a source's result, and whether the caller is to mark it unavailable. */
export interface LoadResult extends Partial<SetInput> {
  /** A `source` node's result, for its readers (`StageContext.sources`). */
  readonly source?: SourceResult;
  /** Unavailable (reported, or quietly when it reads an unavailable node): add it to `StageContext.unavailable`. */
  readonly unavailable?: boolean;
}
/** A table's columns: values, or typed after `fields` (seqspec:format#values). */
type Column = readonly unknown[] | Int32Array | Float64Array;
type Table = { columns: Record<string, Column>; length: number };

const MiB = 1 << 20;
const ROWS_AT_A_TIME = 8;
const messageOf = (e: unknown) => (e instanceof Error ? e.message : String(e));
const encoder = new TextEncoder();
const head = (ref: string) => ref.split("/")[0];
/** The object without null and undefined values: model inputs take absent, never null. */
const clean = <T extends object>(o: T): T => Object.fromEntries(Object.entries(o).filter(([, v]) => v != null)) as T;
const node = (doc: DocIndex, n: DocNode) => doc.nodes.get(n.ref ?? n.kind) as IndexedNode;
const pathOf = (n: IndexedNode, param?: string) => (param ? `${n.path}/params/${param}` : n.path);
const UNAVAILABLE: LoadResult = { unavailable: true };
function unavailable(ctx: StageContext, n: IndexedNode, why: string): LoadResult {
  report(ctx, "svs.provider.unavailable", n.path, `"${n.ref}" is unavailable: ${why}`);
  return UNAVAILABLE;
}

const REF_PARAMS =
  "source_ref from_source sequence ligand from to target alignment of reference pairs partner options query items order";
/** Whether a ref param names an unavailable node (then this node is unavailable too, quietly). */
function blocked(ctx: StageContext, n: IndexedNode): boolean {
  const gone = ctx.unavailable;
  const hit = (v: unknown): boolean =>
    typeof v === "string"
      ? gone.has(head(v))
      : Array.isArray(v)
        ? v.some(hit)
        : v !== null && typeof v === "object" && Object.values(v).some(hit);
  return gone.size > 0 && REF_PARAMS.split(" ").some((k) => hit(n.params[k]));
}

/**
 * `StageContext.source` over the loaded sources: a source's default parse, by its ref or a data node's read from it (an
 * annotation's `uri` its implied source).
 */
export function sourceParse(doc: DocIndex, sources: ReadonlyMap<string, SourceResult>, ref: string) {
  const n = doc.nodes.get(head(ref));
  const implied = n && impliedSource(n, doc);
  const r =
    sources.get(n?.kind === "source" ? n.ref : implied ? implied.ref : String(n?.params.source_ref)) ??
    (implied ? sources.get(n.ref) : undefined);
  return r?.kind === "loaded" ? r.parsed : null;
}

/** A source's bytes, loaded once; its parses by request, each once, reported at its path. */
async function open(src: IndexedNode, ctx: StageContext, loaders: Loaders): Promise<SourceResult> {
  const p = src.params;
  const c = at(ctx, src.path);
  const fail = (code: string, m: string): SourceResult => {
    report(c, code, "", m);
    return { kind: code === "svs.provider.unavailable" ? "unavailable" : "failed" };
  };
  let uri = (p.uri as string | null) ?? null;
  let format = p.format as string | null;
  const service = p.service == null ? null : ctx.registries.get("services", p.service as string);
  const [id, resource] = [String(p.id), String(p.resource ?? service?.resources[0])];
  if (p.service != null) {
    if (!service) return fail("svs.provider.unavailable", `service "${p.service}" is not registered`);
    let r: ReturnType<typeof service.resolve>;
    try {
      r = service.resolve(id, resource);
    } catch (e) {
      return fail("svs.source.load", `${service.name}: ${messageOf(e)}`);
    }
    if (!r) return fail("svs.source.load", `${service.name} has no resource "${resource}" (${service.resources})`);
    uri = r.uri;
    format ??= r.format;
  }
  // A format no package registered: unavailable, so a document opens without the extension
  // (seqspec:sources#data-by-reference).
  const entry = format == null ? undefined : ctx.registries.get("formats", format);
  if (!entry) return fail("svs.provider.unavailable", `format "${format}" is not registered`);
  const req = {
    uri,
    hostFile: (p.host_file as string | null) ?? null,
    maxBytes: (p.max_bytes as number | null) ?? null,
  };
  let shared: FormatResult | null = null;
  let bytes: Uint8Array = new Uint8Array(0);
  try {
    if (p.value == null && loaders.parsed) shared = await loaders.parsed(req, entry.name, ctx.signal);
    if (!shared) {
      bytes = p.value != null ? encoder.encode(p.value as string) : await loaders.load(req, ctx.signal);
      const max = req.maxBytes ?? (p.value != null ? MiB : null);
      if (max != null && bytes.length > max) return fail("svs.source.size_bound", `${bytes.length} bytes, over ${max}`);
    }
  } catch (e) {
    ctx.signal.throwIfAborted();
    return fail("svs.source.load", `${uri ?? req.hostFile ?? "the value"}: ${messageOf(e)}`);
  }
  // A read that reports an error is unreadable (null): its readers are unavailable, quietly. A source's issues are
  // said once however many reads of it meet them (GV-02: a FASTA alignment's, twice, on two nodes).
  const said = new Set<string>();
  const parse = async (request: FormatRequest, rc: StageContext) => {
    const c = at(rc, src.path);
    let failed = false;
    const counted: StageContext = Object.assign(Object.create(c), {
      report: (i: Issue) => {
        failed ||= i.severity === "error";
        const key = `${i.code}\0${i.path}\0${i.message}`;
        if (said.has(key)) return;
        said.add(key);
        c.report(i);
      },
    });
    try {
      const result = shared ?? (await entry.read(bytes, request, counted));
      return failed ? null : service?.after ? service.after(result, id, resource, counted) : result;
    } catch (e) {
      rc.signal.throwIfAborted();
      report(rc, "svs.source.load", src.path, `${uri ?? req.hostFile ?? "the value"}: ${messageOf(e)}`);
      return null;
    }
  };
  const reads = new Map<string, Promise<FormatResult | null>>();
  const read = (request: FormatRequest, rc: StageContext) => {
    // One parse per source and distinct values of the request fields the format reads; `want` null is its first. A
    // record and a geometry choose a tag's reading only: N alignments over N records of one file parse it once.
    const key = JSON.stringify(
      entry.reads.map((k) =>
        k === "want"
          ? (request.want ?? entry.gives[0])
          : (k === "record" || k === "geometry") && request.tag == null
            ? null
            : request[k],
      ),
    );
    const r = reads.get(key) ?? parse(request, rc);
    reads.set(key, r);
    return r;
  };
  return { kind: "loaded", format: entry.name, parsed: null, read };
}

/** The readers' params a provider with no manifest of its own is told (`StreamSource.declared`). */
const SHAPE = ["row_count", "column_count", "length", "target"];
/**
 * A streamed source (seqspec:sources#streamed-sources), by `stream` or by `protocol` with `uri`: its provider
 * (`loaders.stream`) asked for the manifest, and nothing else; its readers' shape params `declared`, and with
 * `protocol` its `data_kind` the first reader's. No provider: unavailable; a manifest the provider does not serve:
 * unavailable; a failed or malformed one: failed (reported).
 */
async function streamSource(
  n: IndexedNode,
  ctx: StageContext,
  loaders: Loaders,
  doc?: DocIndex,
  parent: string | null = null,
): Promise<SourceResult> {
  const p = n.params;
  const readers = (doc?.order ?? []).filter((r) => r.params.source_ref === n.ref);
  const declared: Record<string, unknown> = { ...(p.declared as object) };
  for (const r of readers) for (const k of SHAPE) if (r.params[k] != null) declared[k] = r.params[k];
  if (p.coordinate_system != null) declared.target = p.coordinate_system;
  const first = readers[0];
  const kind = (p.data_kind ??
    (first?.kind === "alignment"
      ? "alignment"
      : first?.params.geometry === "values"
        ? "values"
        : first?.params.geometry === "ranges"
          ? "items"
          : null)) as StreamDataKind | null;
  // A protocol source nothing reads is not asked.
  if (kind === null) return { kind: "unavailable" };
  const source: StreamSource = {
    ref: n.ref,
    stream: (p.stream as string | null) ?? null,
    protocol: (p.protocol as string | null) ?? null,
    uri: (p.uri as string | null) ?? null,
    data_kind: kind,
    parent,
    declared,
  };
  const name = source.stream != null ? `stream "${source.stream}"` : `${source.protocol} source ${source.uri}`;
  const provider = loaders.stream?.(source) ?? null;
  if (!provider) {
    unavailable(ctx, n, `the host provides no ${name}`);
    return { kind: "unavailable" };
  }
  let m: Manifest;
  try {
    const patch = await provider.fetch(streamRequest(n.ref, null), ctx.signal);
    if (patch.kind !== "manifest") throw new Error(`its first answer is ${patch.kind}, not a manifest`);
    m = patch.manifest;
  } catch (e) {
    ctx.signal.throwIfAborted();
    if (e instanceof StreamNotServed) {
      unavailable(ctx, n, `${name}: ${messageOf(e)}`);
      return { kind: "unavailable" };
    }
    report(ctx, "svs.source.load", n.path, `${name}: the manifest: ${messageOf(e)}`);
    return { kind: "failed" };
  }
  const fault =
    !Array.isArray(m?.systems) || !Array.isArray(m.levels)
      ? "a manifest has systems and levels"
      : kind === "values" && !Array.isArray(m.domain)
        ? "a values manifest has a domain"
        : kind === "alignment" && !((m.rows ?? 0) > 0 && (m.columns ?? 0) > 0)
          ? "an alignment manifest has rows and columns"
          : null;
  if (fault) {
    report(ctx, "svs.source.format", n.path, `${name}: ${fault}`);
    return { kind: "failed" };
  }
  return { kind: "streamed", source, manifest: m, provider };
}

/** A node's streamed source by its `source_ref`; null: none. A reader of another data kind is reported. */
function streamedOf(n: IndexedNode, doc: DocIndex, ctx: StageContext, want: StreamDataKind | null): Streamed | null {
  const src = sourceOf(doc, n.params.source_ref);
  const sr = src && ctx.sources.get(src.ref);
  if (sr?.kind !== "streamed") return null;
  if (want === sr.source.data_kind) return sr;
  const m = `a ${sr.source.data_kind} stream is not read by this ${n.kind}${want ? ` (it reads ${want})` : ""}`;
  report(ctx, "svs.source.format", pathOf(n, "source_ref"), m);
  return null;
}

/**
 * A lazy molecule of a `letters` stream's coordinate system: its length and numbering the manifest's
 * (seqspec:model#molecules), the stream as its `svs:stream` (its one system, by its key: the viewer draws its letters
 * from the regions, ST3a).
 */
function lazySequence(key: string, s: Streamed, sys: Manifest["systems"][number], alphabet: unknown): MoleculeInput {
  const data: StreamedData = { ...s, systems: new Map([[key, s.manifest.systems.indexOf(sys)]]), labels: new Map() };
  return clean({
    key,
    kind: "polymer" as const,
    alphabet: (alphabet ?? s.manifest.alphabet ?? undefined) as MoleculeInput["alphabet"],
    lazy: { length: sys.extent[1] },
    label: sys.label ?? undefined,
    numberings: sys.numbering ? [numberingOf(sys.numbering, undefined)] : undefined,
    props: { [StreamProp]: data },
  });
}

/** A source node, or the source an annotation's `uri` implies: its default parse, or its streamed provider. */
async function loadSource(n: IndexedNode, ctx: StageContext, loaders: Loaders, doc?: DocIndex): Promise<SourceResult> {
  // A nested source (seqspec:sources#dynamic-nodes) sits under the source that added it.
  const parent = doc?.nodes.get(n.parent ?? "")?.kind === "source" ? n.parent : null;
  if (n.params.stream != null || n.params.protocol != null) return streamSource(n, ctx, loaders, doc, parent);
  const r = await open(n, ctx, loaders);
  return r.kind === "loaded" ? { ...r, parsed: await r.read(emptyRequest(), ctx) } : r;
}

/** The source a node reads (`source_ref`), or null (reported) when it names none. */
function sourceOf(doc: DocIndex, ref: unknown): IndexedNode | null {
  const n = typeof ref === "string" ? doc.nodes.get(ref) : undefined;
  return n?.kind === "source" ? n : null;
}

/**
 * The source a node implies: an annotation's inline `uri`, `<ref>#uri`; and the statistic a `column_stats` computes
 * over a streamed alignment with the default params, `<ref>#stats`, a `values` source of its own asked of the host's
 * stream providers like the alignment's, told its `measure`, the columns and the alignment. Null: it has none.
 */
function impliedSource(n: IndexedNode, doc: DocIndex): IndexedNode | null {
  const p = n.params;
  if (n.kind === "annotation" && p.uri != null)
    return { ...n, ref: `${n.ref}#uri`, kind: "source", params: { uri: p.uri, format: p.format } };
  if (n.kind !== "computed" || p.provider !== "column_stats" || p.measure === "consensus") return null;
  if (p.row_subset != null || p.reference_row != null || p.gaps !== "count") return null;
  const aln = doc.nodes.get(String(p.alignment));
  const src = sourceOf(doc, aln?.params.source_ref);
  if (!aln || !src || (src.params.stream == null && src.params.protocol == null)) return null;
  const declared = clean({ measure: p.measure, length: aln.params.column_count, target: aln.ref });
  return { ...n, ref: `${n.ref}#stats`, kind: "source", params: { ...src.params, data_kind: "values", declared } };
}

/** One record of a source's list: by id (`key`), by index, or the only one. */
function pick<T extends { key: string }>(ctx: StageContext, n: IndexedNode, list: T[], record: unknown): T | null {
  const r =
    typeof record === "string"
      ? list.find((x) => x.key === record)
      : record == null
        ? list.length === 1
          ? list[0]
          : null
        : list[Number(record)];
  if (r) return r;
  const ids = list.map((x, k) => x.key || k).join(", ");
  const why = record == null ? `has ${list.length} records; set record` : `has no record ${JSON.stringify(record)}`;
  report(ctx, "svs.source.record", pathOf(n, "record"), `the source ${why} (records: ${ids || "none"})`);
  return null;
}

/** A node's source read with its request (null: none, or unreadable, reported). */
async function read(n: IndexedNode, doc: DocIndex, ctx: StageContext, request: Partial<FormatRequest>) {
  const src = sourceOf(doc, n.params.source_ref);
  const sr = src && ctx.sources.get(src.ref);
  if (sr?.kind === "streamed") streamedOf(n, doc, ctx, null);
  if (sr?.kind !== "loaded") return null;
  const r = await sr.read(emptyRequest(request), ctx);
  if (r && r.kind !== request.want) {
    report(ctx, "svs.source.format", pathOf(n, "source_ref"), `the source gives ${r.kind}, not ${request.want}`);
    return null;
  }
  return r;
}

/** A `numbering` node's params (seqspec:model#numberings) as the model's input: a node's, or a manifest system's. */
function numberingOf(p: P, custom: object | undefined): NumberingInput {
  const base = { name: p.name as string, isDefault: p.default === true, custom };
  if (p.method === "offset")
    return { ...base, method: "offset", first: p.first as number, skipZero: p.skip_zero === true };
  if (p.method === "segments")
    return { ...base, method: "segments", segments: p.segments as [number, number, number][] };
  if (p.method === "names") return { ...base, method: "names", names: p.labels as (string | null)[] };
  return { ...base, method: "labels", labels: p.labels as (number | string | null)[] };
}
function numberings(n: IndexedNode, doc: DocIndex): NumberingInput[] {
  return (n.node.children ?? [])
    .filter((c) => c.kind === "numbering")
    .map((c) => numberingOf(node(doc, c).params, c.custom));
}

const provenance = (p: P) => (p.provenance ?? undefined) as MoleculeInput["provenance"];

async function molecule(n: IndexedNode, doc: DocIndex, ctx: StageContext): Promise<LoadResult> {
  const p = n.params;
  const ligand = n.kind === "ligand";
  let m: MoleculeInput = { key: n.ref };
  const streamed = ligand ? null : streamedOf(n, doc, ctx, "letters");
  if (p.value != null)
    m = ligand ? { key: n.ref, elements: p.value as string[] } : { key: n.ref, residues: p.value as string };
  else if (streamed) {
    const sys = pick(ctx, n, streamed.manifest.systems, p.record);
    if (!sys) return UNAVAILABLE;
    m = lazySequence(n.ref, streamed, sys, p.alphabet);
  } else {
    const r = await read(n, doc, ctx, { want: "molecules", record: p.record as string | number | null });
    const got = r?.kind === "molecules" ? pick(ctx, n, r.molecules, p.record) : null;
    if (!got) return UNAVAILABLE;
    m = got;
  }
  const input: MoleculeInput = clean({
    ...m,
    key: n.ref,
    kind: ligand ? "ligand" : "polymer",
    alphabet: (p.alphabet ?? m.alphabet) as MoleculeInput["alphabet"],
    cyclic: p.cyclic === true || m.cyclic,
    chemCompIds: (p.chem_comp_ids ?? m.chemCompIds) as (string | null)[],
    label: (p.label ?? m.label) as string,
    identifiers: [...(m.identifiers ?? []), ...((p.identifiers as MoleculeInput["identifiers"]) ?? [])],
    numberings: [...(m.numberings ?? []), ...numberings(n, doc)],
    provenance: provenance(p) ?? m.provenance,
    ownUnit: true,
    custom: n.node.custom,
  });
  return { molecules: [input] };
}

function unit(n: IndexedNode, doc: DocIndex): SetPart {
  const p = n.params;
  const parent = doc.nodes.get(n.parent ?? "");
  const origin = (parent?.kind === "complex" ? parent.params.origin : null) ?? p.origin;
  const u: UnitInput = clean({
    key: n.ref,
    molecule: (p.sequence ?? p.ligand) as string,
    range: p.range as [number, number],
    observed: p.observed as [number, number][],
    label: p.label as string,
    origin: origin as string,
    score: p.score as number,
    numberings: numberings(n, doc),
    custom: n.node.custom,
  });
  return { units: [u] };
}

/**
 * A tree computed by its registered node provider over the alignment `a` (null: a streamed alignment, whose rows are not
 * held): the provider's `TreeInput`, or undefined with the reason reported. A newer request aborts it: that propagates.
 */
async function providerTree(tn: IndexedNode, ctx: StageContext, a: AlignmentInput | null) {
  const name = String(tn.params.provider);
  const entry = ctx.registries.get("nodeProviders", name);
  const why = (m: string) => unavailable(ctx, tn, m);
  if (!entry) return void why(`node provider "${name}" is not registered`);
  if (entry.kind !== "tree") return void why(`node provider "${name}" computes a ${entry.kind}, not a tree`);
  if (!a) {
    report(ctx, "svs.provider.lazy", tn.path, `"${tn.ref}": a tree is not computed over a streamed alignment's rows`);
    return undefined;
  }
  try {
    return (await entry.compute(at(ctx, tn.path), [a], tn.params.options, ctx.signal)) as TreeInput | null;
  } catch (e) {
    ctx.signal.throwIfAborted();
    if (e instanceof TreeSizeError) report(ctx, "svs.provider.size", tn.path, `"${tn.ref}": ${e.message}`);
    else why(`node provider "${name}" failed: ${messageOf(e)}`);
    return undefined;
  }
}

async function treeOf(
  n: IndexedNode,
  doc: DocIndex,
  ctx: StageContext,
  labels: (string | null)[],
  held: AlignmentInput | null = null,
) {
  const t = (n.node.children ?? []).find((c) => c.kind === "tree");
  if (!t) return undefined;
  const tn = node(doc, t);
  if (tn.params.provider != null) return (await providerTree(tn, ctx, held)) ?? undefined;
  let r: FormatResult | null | undefined;
  try {
    r =
      tn.params.value != null
        ? await ctx.registries
            .get("formats", "newick")
            ?.read(encoder.encode(tn.params.value as string), emptyRequest(), at(ctx, tn.path))
        : await read(tn, doc, ctx, { want: "tree" });
  } catch (e) {
    ctx.signal.throwIfAborted();
    report(ctx, "svs.source.format", tn.path, messageOf(e));
    return undefined;
  }
  if (r?.kind !== "tree") return undefined;
  return { ...r.tree, leafRow: Trees.matchLeaves(r.tree, labels, tn.params.leaves as "labels" | "order") };
}

const named = new WeakMap<DocIndex, Promise<Map<string, Set<number>>>>();
/**
 * The rows of every lazy alignment that mappings and positional annotations name by `<alignment>/<row>`, found in one
 * pass over the document's nodes and tables, however many alignments there are.
 */
function namedRows(doc: DocIndex, ctx: StageContext, loaders: Loaders): Promise<Map<string, Set<number>>> {
  const out = new Map<string, Set<number>>();
  const see = (v: unknown) => {
    const m = typeof v === "string" ? /^(.+)\/(\d+)$/.exec(v) : null;
    if (m) out.set(m[1], (out.get(m[1]) ?? new Set()).add(Number(m[2])));
  };
  const scan = async () => {
    for (const n of doc.order) {
      if (n.kind === "mapping") [n.params.from, n.params.to].forEach(see);
      if (n.kind !== "annotation" || n.params.geometry === "rows") continue;
      see(n.params.target);
      const t = await tableOf(n, doc, ctx, loaders, false);
      for (const k of TARGETS) t?.columns[k]?.forEach(see);
    }
    return out;
  };
  let done = named.get(doc);
  if (done === undefined) {
    const own = scan();
    named.set(doc, own);
    own.catch(() => named.get(doc) === own && named.delete(doc));
    done = own;
  }
  return done;
}

/**
 * A lazy alignment: a streamed one (its manifest's shape, the declared counts checked against it); its reference row
 * and the rows named by mappings and positional annotations fetched whole, 8 at a time, labels and tiles from its source.
 */
async function lazyAlignment(n: IndexedNode, src: IndexedNode, doc: DocIndex, ctx: StageContext, loaders: Loaders) {
  const p = n.params;
  const s = streamedOf(n, doc, ctx, "alignment");
  if (!s) return UNAVAILABLE;
  const labels = new Map<number, string>();
  const sys = s.manifest.systems[0]?.key ?? n.ref;
  const name = s.source.stream != null ? `stream "${s.source.stream}"` : (s.source.uri ?? n.ref);
  const row = async (k: number, signal: AbortSignal) => {
    const r = await s.provider.fetch(
      streamRequest(src.ref, sys, { stop: s.manifest.columns ?? 0, row_start: k, row_stop: k + 1 }),
      signal,
    );
    if (r.kind !== "rows" || !r.rows) throw new Error(`a ${r.kind} patch, not rows`);
    labels.set(k, r.labels[0]);
    return { label: r.labels[0], row: r.rows[0] };
  };
  const rows = new Map<number, { label: string; gapped: string }>();
  // The scan is shared: one another caller began and aborted (the viewer's graph gives each node its signal) is redone.
  const scanned = await namedRows(doc, ctx, loaders).catch(() => {
    ctx.signal.throwIfAborted();
    return namedRows(doc, ctx, loaders);
  });
  const wanted = [...new Set([(p.reference_row as number | null) ?? 0, ...(scanned.get(n.ref) ?? [])])];
  wanted.sort((x, y) => x - y);
  // Eight at a time, their results in row order.
  const failed: string[] = [];
  for (let i = 0; i < wanted.length; i += ROWS_AT_A_TIME) {
    const batch = wanted.slice(i, i + ROWS_AT_A_TIME);
    const got = await Promise.all(
      batch.map(async (k) => {
        try {
          return await row(k, ctx.signal);
        } catch (e) {
          return e instanceof Error ? e : new Error(String(e));
        }
      }),
    );
    ctx.signal.throwIfAborted();
    got.forEach((r, j) => {
      if (r instanceof Error) failed.push(`row ${batch[j]}: ${r.message}`);
      else rows.set(batch[j], { label: r.label, gapped: r.row });
    });
  }
  // One issue for every row that failed, with the first few (GV-07: an outage gave one per row scanned).
  if (failed.length > 0) {
    const more = failed.length > 3 ? ` and ${failed.length - 3} more` : "";
    const message = `${name}: ${failed.length} row(s) failed: ${failed.slice(0, 3).join("; ")}${more}`;
    report(ctx, "svs.source.load", pathOf(n, "source_ref"), message);
  }
  const [rowCount, columnCount] = [s.manifest.rows ?? 0, s.manifest.columns ?? 0];
  for (const [k, v] of [
    ["row_count", rowCount],
    ["column_count", columnCount],
  ] as const)
    if (p[k] != null && p[k] !== v)
      report(ctx, "svs.alignment.shape", pathOf(n, k), `${k} is ${p[k]}; the manifest's ${v}`);
  const a: AlignmentInput = clean({
    key: n.ref,
    gapped: [],
    flavor: "aligned",
    lazy: { rowCount, columnCount, rows, labelOf: (k: number) => labels.get(k) ?? null },
    tree: await treeOf(n, doc, ctx, [], null),
    alphabet: p.alphabet as AlignmentInput["alphabet"],
    referenceRow: p.reference_row as number,
    parentRows: p.parent_rows as (number | null)[],
    custom: n.node.custom,
    props: { [StreamProp]: { ...s, systems: new Map([[n.ref, 0]]), labels } satisfies StreamedData },
  });
  return { alignments: [a] };
}

/**
 * What a computed alignment's `inputs` read in the set: per input its residues for the provider, and the key its new row
 * is bound to; null with the reason reported. A list names sequences, units (their range) or rows; one alignment key
 * names its rows, bound to the old rows.
 */
function realignInputs(n: IndexedNode, ctx: StageContext) {
  const set = ctx.set;
  const p = n.params;
  const say = (code: string, message: string) => {
    report(ctx, code, pathOf(n, "inputs"), message);
    return null;
  };
  if (!set) return say("svs.ref.derived_unknown", "nothing is loaded to align");
  const keys = typeof p.inputs === "string" ? [p.inputs] : (p.inputs as string[]);
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  const inputs: RealignInput[] = [];
  const targets: string[] = [];
  const see = (
    key: string,
    label: string | null,
    alphabet: AlignerInput["alphabet"],
    code: Uint8Array,
    bind: string,
  ) => {
    inputs.push({ key: set.strings.find(key), kind: MoleculeKind.polymer, alphabet, length: code.length, code, label });
    targets.push(bind);
  };
  const fromRows = (aln: Alignment, key: string, rows: Iterable<number>) => {
    for (const k of rows) {
      const slot = aln.rows.slotOf(k as never);
      if (slot < 0) return say("svs.provider.lazy", `"${key}/${k}" is a row that is not held (a streamed alignment)`);
      const code = aln.rows.residues.subarray(aln.rows.residueStart[slot], aln.rows.residueStart[slot + 1]);
      see(`${key}/${k}`, aln.labelText(k as never), aln.alphabet, code, `${key}/${k}`);
    }
    return true;
  };
  for (const key of keys) {
    if (key === n.ref || key.startsWith(`${n.ref}/`))
      return say("svs.ref.cycle", `"${key}" is this alignment's own: it cannot be aligned from itself`);
    if (!set.byName(key, t)) return say("svs.ref.derived_unknown", `"${key}" names nothing loaded to align`);
    if (t.kind === KeyKind.alignment) {
      const aln = set.alignments[t.id] as Alignment;
      if (
        !fromRows(
          aln,
          key,
          Array.from({ length: aln.rowCount }, (_, k) => k),
        )
      )
        return null;
    } else if (t.kind === KeyKind.row) {
      if (!fromRows(set.alignments[t.id] as Alignment, key.slice(0, key.lastIndexOf("/")), [t.row])) return null;
    } else if (t.kind === KeyKind.molecule || t.kind === KeyKind.unit) {
      const u = t.kind === KeyKind.unit;
      const m = set.molecules[u ? set.units.molecule[t.id] : t.id] as Molecule;
      if (m.kind !== MoleculeKind.polymer)
        return say("svs.ligand.unsupported", `"${key}" is a ligand: it cannot be aligned`);
      if (m.lazy) return say("svs.provider.lazy", `"${key}" is a streamed sequence: its letters are not held`);
      const code = u ? m.code.subarray(set.units.rangeStart[t.id], set.units.rangeStop[t.id]) : m.code;
      const label = m.label >= 0 ? set.strings.get(m.label) : null;
      see(key, label, m.alphabet, code, key);
    } else return say("svs.ref.wrong_kind", `"${key}" is no sequence, unit, row or alignment`);
  }
  return {
    inputs,
    targets,
    reference: typeof p.inputs === "string" ? ((set.alignments[t.id] as Alignment).referenceRow as number) : 0,
  };
}

/**
 * An alignment computed by its registered node provider (`realign`) from the set's sequences, units or rows: the
 * provider's `AlignmentInput` under this node's key, and one `identity` mapping per row from the new row to its input,
 * so whatever is bound to the inputs reaches the rows (seqspec:format#extensions; seqspec:model#alignments-and-trees).
 * Runs once the set holds the inputs, as an `align` mapping does; null from the provider is an unavailable node (its
 * reason reported).
 */
async function computedAlignment(n: IndexedNode, doc: DocIndex, ctx: StageContext): Promise<LoadResult> {
  const p = n.params;
  const name = String(p.provider);
  const entry = ctx.registries.get("nodeProviders", name);
  if (!entry) return unavailable(ctx, n, `node provider "${name}" is not registered`);
  if (entry.kind !== "alignment")
    return unavailable(ctx, n, `node provider "${name}" computes a ${entry.kind}, not an alignment`);
  const read = realignInputs(n, ctx);
  if (!read) return UNAVAILABLE;
  let result: Awaited<ReturnType<typeof entry.compute>>;
  try {
    result = await entry.compute(at(ctx, n.path), read.inputs, p.options, ctx.signal);
  } catch (e) {
    ctx.signal.throwIfAborted();
    if (e instanceof AlignerSizeError) report(ctx, "svs.provider.size", n.path, `"${n.ref}": ${e.message}`);
    else unavailable(ctx, n, `node provider "${name}" failed: ${messageOf(e)}`);
    return UNAVAILABLE;
  }
  if (!result) return UNAVAILABLE;
  const a: AlignmentInput = clean({
    ...(result as AlignmentInput),
    key: n.ref,
    labels: (p.labels ?? (result as AlignmentInput).labels) as string[],
    alphabet: (p.alphabet ?? (result as AlignmentInput).alphabet) as AlignmentInput["alphabet"],
    referenceRow: (p.reference_row ?? read.reference) as number,
    parentRows: p.parent_rows as (number | null)[],
    custom: n.node.custom,
  });
  const tree = await treeOf(n, doc, ctx, a.labels ?? a.gapped.map(() => null), a);
  if (tree) a.tree = tree;
  const mappings: MappingInput[] = read.targets.map((to, k) => ({
    key: `${n.ref}.row-${k}`,
    from: `${n.ref}/${k}`,
    to,
    method: "identity",
  }));
  return { alignments: [a], mappings };
}

async function alignment(n: IndexedNode, doc: DocIndex, ctx: StageContext, loaders: Loaders): Promise<LoadResult> {
  const p = n.params;
  if (p.provider != null) return computedAlignment(n, doc, ctx);
  const src = sourceOf(doc, p.source_ref);
  if (src?.params.protocol != null || src?.params.stream != null) return lazyAlignment(n, src, doc, ctx, loaders);
  let a: AlignmentInput;
  if (p.value != null) a = { key: n.ref, gapped: p.value as string[], flavor: p.flavor as "aligned" | "a3m" };
  else {
    const r = await read(n, doc, ctx, {
      want: "alignments",
      flavor: p.flavor as "aligned",
      record: p.record as number,
    });
    const got = r?.kind === "alignments" ? pick(ctx, n, r.alignments, p.record) : null;
    if (!got) return UNAVAILABLE;
    a = got;
  }
  a = clean({
    ...a,
    key: n.ref,
    labels: (p.labels ?? a.labels) as string[],
    alphabet: (p.alphabet ?? a.alphabet) as AlignmentInput["alphabet"],
    referenceRow: p.reference_row as number,
    parentRows: p.parent_rows as (number | null)[],
    custom: n.node.custom,
  });
  const width = a.flavor === "a3m" ? (a.gapped[0] ?? "").replace(/[a-z.]/g, "").length : (a.gapped[0]?.length ?? 0);
  for (const [k, v] of [
    ["row_count", a.gapped.length],
    ["column_count", width],
  ] as const) {
    if (p[k] != null && p[k] !== v)
      report(ctx, "svs.alignment.shape", pathOf(n, k), `${k} is ${p[k]}; the data has ${v}`);
  }
  const tree = await treeOf(n, doc, ctx, a.labels ?? a.gapped.map(() => null), a);
  if (tree) a.tree = tree;
  return { alignments: [a] };
}

/** What an aligner reads of an endpoint: a molecule (a unit's, whose space uses its molecule's indices) or a row. */
function alignerInput(set: SeqSet, key: string) {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  if (!set.byName(key, t)) return null;
  if (t.kind === KeyKind.molecule || t.kind === KeyKind.unit) {
    return set.molecules[t.kind === KeyKind.unit ? set.units.molecule[t.id] : t.id] ?? null;
  }
  const a = set.alignments[t.id];
  const slot = t.kind === KeyKind.row && a ? a.rows.slotOf(t.row as never) : -1;
  if (!a || slot < 0) return null;
  const code = a.rows.residues.subarray(a.rows.residueStart[slot], a.rows.residueStart[slot + 1]);
  return { key: set.strings.find(key), kind: 0, alphabet: a.alphabet, length: code.length, code };
}

async function mapping(n: IndexedNode, ctx: StageContext): Promise<LoadResult> {
  const p = n.params;
  const base = clean({
    key: n.ref,
    from: p.from as string,
    to: p.to as string,
    provenance: provenance(p),
    custom: n.node.custom,
  });
  let m: MappingInput;
  if (p.method === "identity") m = clean({ ...base, method: "identity", fromRange: p.from_range as [number, number] });
  else if (p.method === "segments")
    m = { ...base, method: "segments", segments: p.segments as [number, number, number | null][] };
  else if (p.method === "pairwise") m = { ...base, method: "pairwise", cigar: p.cigar as string };
  else if (p.method === "by_numbering") {
    m = clean({
      ...base,
      method: "by_numbering",
      numbering: p.numbering as string,
      toNumbering: p.to_numbering as string,
    });
  } else if (p.method === "codon") {
    m = clean({
      ...base,
      method: "codon",
      cds: p.cds as [number, number][],
      start: p.start as number,
      fromRange: p.from_range as [number, number],
      strand: p.strand as "+" | "-",
      geneticCode: p.genetic_code as number,
    });
  } else if (p.method === "provider") {
    return mappingProvider(n, ctx);
  } else {
    const entry = ctx.registries.get("aligners", p.aligner as string);
    if (!entry) return unavailable(ctx, n, `aligner "${p.aligner}" is not registered`);
    const [a, b] = ctx.set ? [alignerInput(ctx.set, base.from), alignerInput(ctx.set, base.to)] : [null, null];
    if (!a || !b) {
      const which = a ? "to" : "from";
      report(ctx, "svs.ref.derived_unknown", pathOf(n, which), `"${p[which]}" names nothing loaded to align`);
      return UNAVAILABLE;
    }
    let cigar: string | null;
    try {
      cigar = await entry.align(at(ctx, n.path), a, b, p.options);
    } catch (e) {
      ctx.signal.throwIfAborted();
      report(ctx, "svs.mapping.align", n.path, `${entry.name}: ${messageOf(e)}`);
      return UNAVAILABLE;
    }
    // A null is an aligner that did not align, its reason reported by the entry: the mapping is unavailable.
    if (cigar == null) return UNAVAILABLE;
    // A given provenance keeps its fields and gains the aligner and its options under `parameters`.
    const given = base.provenance;
    const parameters = { ...(given?.parameters as object | undefined), aligner: entry.name, options: p.options };
    m = {
      ...base,
      method: "pairwise",
      cigar,
      provenance: given ? { ...given, parameters } : { source: entry.name, method: "alignment", parameters },
    };
  }
  return { mappings: [m] };
}

/** An alignment's key and row count, as a mapping provider reads its endpoints; null with the reason reported. */
function summaryOf(n: IndexedNode, ctx: StageContext, param: "from" | "to"): AlignmentSummary | null {
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  const key = String(n.params[param]);
  const set = ctx.set;
  if (!set?.byName(key, t) || t.kind !== KeyKind.alignment) {
    const why = set?.byName(key, t) ? "is not an alignment" : "names nothing loaded";
    report(
      ctx,
      "svs.ref.derived_unknown",
      pathOf(n, param),
      `"${key}" ${why}: the provider pairs the rows of two alignments`,
    );
    return null;
  }
  const a = set.alignments[t.id] as Alignment;
  if (a.lazy) {
    report(ctx, "svs.provider.lazy", pathOf(n, param), `"${key}" is a streamed alignment: its rows are not held`);
    return null;
  }
  return { key, rowCount: a.rowCount };
}

/**
 * A mapping computed by its registered node provider (`codon_rows`) from its two endpoints, alignments read from the set.
 * The provider's mappings have suffixes for keys; this node's ref is put before them.
 */
async function mappingProvider(n: IndexedNode, ctx: StageContext): Promise<LoadResult> {
  const name = String(n.params.provider);
  const entry = ctx.registries.get("nodeProviders", name);
  if (!entry) return unavailable(ctx, n, `node provider "${name}" is not registered`);
  if (entry.kind !== "mapping")
    return unavailable(ctx, n, `node provider "${name}" computes a ${entry.kind}, not a mapping`);
  const [from, to] = [summaryOf(n, ctx, "from"), summaryOf(n, ctx, "to")];
  if (!from || !to) return UNAVAILABLE;
  let result: Awaited<ReturnType<typeof entry.compute>>;
  try {
    result = await entry.compute(at(ctx, n.path), [from, to], n.params.options, ctx.signal);
  } catch (e) {
    ctx.signal.throwIfAborted();
    unavailable(ctx, n, `node provider "${name}" failed: ${messageOf(e)}`);
    return UNAVAILABLE;
  }
  if (!result) return UNAVAILABLE;
  const given = provenance(n.params);
  const list = (Array.isArray(result) ? result : [result]) as MappingInput[];
  return {
    mappings: list.map((m) =>
      clean({ ...m, key: `${n.ref}.${m.key}`, provenance: given ?? m.provenance, custom: n.node.custom }),
    ),
  };
}

const TARGETS = ["target", "b_target"];

/**
 * An annotation's table: from a file, `fields` applied and text cells typed, targets prefixed; null when unreadable.
 * `loud` false: a source that gives no table is not reported again (a second look at it).
 */
async function tableOf(
  n: IndexedNode,
  doc: DocIndex,
  ctx: StageContext,
  loaders: Loaders,
  loud = true,
): Promise<Table | null> {
  const p = n.params;
  let t: Table | null = tableColumns(p.value);
  let prefix: string | null = null;
  if (p.value == null) {
    const request = emptyRequest({
      want: "table",
      tag: p.tag as string | null,
      categoryName: p.category_name as string | null,
      blockHeader: p.block_header as string | null,
      blockIndex: p.block_index as number,
      geometry: p.geometry as FormatRequest["geometry"],
    });
    let src = sourceOf(doc, p.source_ref);
    const inline = impliedSource(n, doc);
    src = inline ?? src;
    if (p.from_source != null) {
      const from = doc.nodes.get(head(p.from_source as string));
      src = sourceOf(doc, from?.params.source_ref);
      if (!src) {
        report(ctx, "svs.source.load", pathOf(n, "from_source"), `"${p.from_source}" was not read from a file`);
        return null;
      }
      request.record = (from?.params.record as number | null) ?? null;
      prefix = from?.ref ?? null;
    }
    const sr =
      inline && src === inline
        ? (ctx.sources.get(inline.ref) ?? (await loadSource(inline, ctx, loaders)))
        : src && ctx.sources.get(src.ref);
    const r = sr?.kind === "loaded" ? await sr.read(request, ctx) : null;
    if (r?.kind !== "table") {
      if (r && loud) report(ctx, "svs.source.format", n.path, `the source gives ${r.kind}, not a table`);
      return null;
    }
    // `fields` read the file's columns as read; a text file's cells are typed after them, then checked as an inline
    // table's are (CF-03: one issue per column, its rows left out).
    const named = renameFields(r.columns, (p.fields ?? {}) as Record<string, string | null>);
    t = { columns: (r.text ? typeColumns(named as never) : named) as Table["columns"], length: r.length };
    if (r.text) t = checkedCells(t, p.geometry as Geometry, n, ctx);
    prefix ??= (p.target as string | null) ?? null;
  }
  if (!t) return null;
  const cols = t.columns;
  // A target read from a file is a suffix of the node it came from, or of the target param.
  if (p.value == null && prefix) {
    for (const k of TARGETS) if (cols[k]) cols[k] = Array.from(cols[k], (v) => (v == null ? v : `${prefix}/${v}`));
  }
  return t;
}

/**
 * A file table's typed columns checked as an inline table's are (`table`'s static check: one `svs.columns.type` per
 * column, its first bad cell and the count), and the rows with a bad index or value cell left out, so a file and the
 * same table inline say one thing (CF-03: the file gave a `[NaN, 3)` per row).
 */
function checkedCells(t: Table, geometry: Geometry, n: IndexedNode, ctx: StageContext): Table {
  // Each index or value cell by the inline check's own rule (LA-10: a number that is no index, -1 or 2.5, passed);
  // a missing cell (null, NaN in a typed column) is not bad. `keep` is made at the first bad cell only.
  let keep: Uint8Array | null = null;
  for (const k of Object.keys(t.columns)) {
    const kind = cellKind(k, geometry);
    if (ColumnRoles[geometry]?.[k] === undefined || (kind !== CellKind.index && kind !== CellKind.value)) continue;
    const c = t.columns[k] as ArrayLike<unknown>;
    for (let i = 0; i < t.length; i++) {
      const v = c[i];
      if (v == null || (typeof v === "number" && Number.isNaN(v)) || cellOk(kind, v)) continue;
      keep ??= new Uint8Array(t.length).fill(1);
      keep[i] = 0;
    }
  }
  if (!keep) return t;
  const kept = keep;
  const check: CheckContext = {
    issues: [],
    registries: ctx.registries,
    refs: [],
    lenient: false,
    from: n.ref,
    depth: 0,
  };
  table(geometry).check(t.columns, `${n.path}/source`, check);
  for (const i of check.issues) if (i.code === "svs.columns.type") ctx.report(i);
  const columns: Table["columns"] = Object.create(null);
  for (const [k, c] of Object.entries(t.columns))
    columns[k] = Array.from(c as ArrayLike<unknown>).filter((_, i) => kept[i] === 1);
  return { columns, length: kept.reduce((n, x) => n + x, 0) };
}

async function annotation(n: IndexedNode, doc: DocIndex, ctx: StageContext, loaders: Loaders): Promise<LoadResult> {
  const p = n.params;
  const base = {
    key: n.ref,
    geometry: p.geometry as AnnotationInput["geometry"],
    semanticType: p.semantic_type as string,
    label: p.label as string,
    target: p.target as string,
    provenance: provenance(p),
    custom: n.node.custom,
  };
  // A streamed annotation holds no rows: its viewer asks its source for regions.
  const src = sourceOf(doc, p.source_ref);
  if (src?.params.stream != null || src?.params.protocol != null) return streamedAnnotation(n, doc, ctx, base);
  // An inline `uri` no driver loaded (the viewer's graph loads node by node): loaded here and returned as this node's
  // source, so `ctx.source(ref)` finds it by the annotation's ref (`sourceParse`).
  const inline = p.value == null ? impliedSource(n, doc) : null;
  const own = inline && !ctx.sources.has(inline.ref) ? await loadSource(inline, ctx, loaders) : null;
  const loaded = own?.kind === "loaded" && inline ? own : null;
  const sources = loaded && inline ? new Map(ctx.sources).set(inline.ref, loaded) : null;
  const t = await tableOf(n, doc, sources ? Object.assign(Object.create(ctx), { sources }) : ctx, loaders);
  if (!t) return UNAVAILABLE;
  // Keyed by data: a `__proto__` column stays a column when the rows are filtered below (CF-02).
  const cols: Table["columns"] = Object.assign(Object.create(null), t.columns);
  // Rows over unavailable nodes are left out quietly.
  const gone = ctx.unavailable;
  const keep =
    gone.size > 0
      ? Array.from({ length: t.length }, (_, i) => !TARGETS.some((k) => gone.has(head(String(cols[k]?.[i] ?? "")))))
      : [];
  if (keep.includes(false)) for (const k of Object.keys(cols)) cols[k] = cols[k].filter((_, i) => keep[i]);
  const a: AnnotationInput = clean({
    ...base,
    numbering: p.numbering as string,
    bNumbering: p.b_numbering as string,
    columns: cols,
    states: p.states as string[],
    links: (p.links as AnnotationInput["links"])?.length ? (p.links as AnnotationInput["links"]) : undefined,
  });
  return loaded ? { annotations: [a], source: loaded } : { annotations: [a] };
}

/**
 * A streamed `items` or `values` annotation: no rows, its coordinate systems (`coordinate_system`, else its target): a
 * complex's member units by the manifest's keys (`hg38/chr17`), or the one system keyed by the ref; a key that names
 * nothing is left out (`svs.stream.ref`).
 */
function streamedAnnotation(
  n: IndexedNode,
  doc: DocIndex,
  ctx: StageContext,
  base: Omit<AnnotationInput, "columns">,
): LoadResult {
  const p = n.params;
  const s = streamedOf(n, doc, ctx, p.geometry === "values" ? "values" : p.geometry === "ranges" ? "items" : null);
  if (!s) return UNAVAILABLE;
  const cs = (sourceOf(doc, p.source_ref)?.params.coordinate_system ?? p.target) as string | null;
  const c = cs == null ? undefined : doc.nodes.get(cs);
  if (cs == null) report(ctx, "svs.params.missing", pathOf(n, "target"), "a streamed annotation needs a target");
  // A complex's members: its units by value, or a letters stream's systems.
  const lettered = c?.kind === "complex" && c.params.provider == null ? sourceOf(doc, c.params.source_ref) : null;
  const ls = lettered && ctx.sources.get(lettered.ref);
  const members = new Set(
    c?.kind !== "complex"
      ? [cs]
      : ls?.kind === "streamed"
        ? ls.manifest.systems.map((x) => `${cs}/${x.key}`)
        : (c.node.children ?? []).map((u) => u.ref),
  );
  const systems = new Map<string, number>();
  s.manifest.systems.forEach((sys, i) => {
    const t = c?.kind === "complex" ? `${cs}/${sys.key}` : sys.key;
    if (members.has(t)) systems.set(t, i);
    else
      report(
        ctx,
        "svs.stream.ref",
        pathOf(n, "source_ref"),
        `the manifest's "${sys.key}" names no ${cs} coordinate system`,
      );
  });
  if (systems.size === 0) return UNAVAILABLE;
  const props = { [StreamProp]: { ...s, systems, labels: new Map() } satisfies StreamedData };
  return { annotations: [clean({ ...base, target: undefined, columns: {}, streamed: [...systems.keys()], props })] };
}

/**
 * A data node's part of the set: a molecule with its numberings, a unit, an alignment with its tree, a mapping (an
 * `align` one through its aligner, with `ctx.set` holding its endpoints), an annotation. A source gives its
 * `SourceResult`, which the caller passes to its readers in `ctx.sources`. Readers of a node in `ctx.unavailable` are
 * unavailable, quietly. Nothing is kept between calls.
 */
export async function loadNode(dn: DocNode, doc: DocIndex, ctx: StageContext, loaders: Loaders): Promise<LoadResult> {
  const n = node(doc, dn);
  if (n.kind === "source") {
    const source = await loadSource(n, ctx, loaders, doc);
    // A source that failed is unavailable like one no package reads: its readers are left out quietly.
    return { source, unavailable: source.kind === "unavailable" || source.kind === "failed" };
  }
  if (blocked(ctx, n)) return UNAVAILABLE;
  if (n.kind === "sequence" || n.kind === "ligand") return molecule(n, doc, ctx);
  if (n.kind === "unit") return unit(n, doc);
  if (n.kind === "alignment") return alignment(n, doc, ctx, loaders);
  if (n.kind === "mapping") return mapping(n, ctx);
  if (n.kind === "annotation") return annotation(n, doc, ctx, loaders);
  return {};
}

/** A complex: its provider's molecules and units, keys prefixed `R/`, or its `unit` children as one complex. */
export async function buildComplex(dn: DocNode, doc: DocIndex, ctx: StageContext): Promise<LoadResult> {
  const n = node(doc, dn);
  const p = n.params;
  const custom = n.node.custom;
  if (p.provider == null && p.source_ref != null) {
    // From a letters stream: per system of its manifest a unit `R/<key>` over a lazy `R/<key>/sequence`.
    const s = streamedOf(n, doc, ctx, "letters");
    if (!s) return UNAVAILABLE;
    const R = (k: string) => `${n.ref}/${k}`;
    const sys = s.manifest.systems;
    const molecules = sys.map((x) => ({ ...lazySequence(R(`${x.key}/sequence`), s, x, null), ownUnit: false }));
    const units = sys.map((x) =>
      clean({ key: R(x.key), molecule: R(`${x.key}/sequence`), label: x.label ?? x.key, origin: p.origin as string }),
    );
    return {
      molecules,
      units,
      complexes: [clean({ key: n.ref, units: units.map((u) => u.key), label: p.label as string, custom })],
    };
  }
  if (p.provider == null) {
    // By value, whatever the document's order: the units that are unavailable are left out quietly, and with none left
    // the complex is unavailable.
    const all = (n.node.children ?? []).filter((c) => c.kind === "unit").map((c) => node(doc, c).ref);
    const units = all.filter((u) => !ctx.unavailable.has(u));
    if (all.length > 0 && units.length === 0) return UNAVAILABLE;
    return { complexes: [clean({ key: n.ref, units, label: p.label as string, custom })] };
  }
  const entry = ctx.registries.get("complexProviders", p.provider as string);
  if (!entry) return unavailable(ctx, n, `complex provider "${p.provider}" is not registered`);
  const src = sourceOf(doc, p.source_ref);
  const sr = src && ctx.sources.get(src.ref);
  const source = sr?.kind === "loaded" ? sr.parsed : null;
  if (!src || !source) return UNAVAILABLE;
  // The format that read the source, a service's included: a provider never runs on one it does not read.
  const format = sr?.kind === "loaded" ? sr.format : null;
  if (format == null || !entry.formats.includes(format)) {
    report(ctx, "svs.source.format", pathOf(n, "source_ref"), `${entry.name} reads ${entry.formats}, not ${format}`);
    return UNAVAILABLE;
  }
  let out: Awaited<ReturnType<typeof entry.build>>;
  try {
    out = await entry.build(at(ctx, n.path), source, p.options);
  } catch (e) {
    ctx.signal.throwIfAborted();
    return unavailable(ctx, n, `provider "${p.provider}" failed: ${messageOf(e)}`);
  }
  const R = (k: string) => `${n.ref}/${k}`;
  const molecules = out.molecules.map((m) => ({ ...m, key: R(m.key), ownUnit: false }));
  const units = out.units.map((u) =>
    clean({
      ...u,
      key: R(u.key),
      label: u.label ?? u.key,
      molecule: R(u.molecule),
      origin: (p.origin ?? u.origin) as string,
    }),
  );
  const complex = clean({ key: n.ref, units: units.map((u) => u.key), label: p.label as string, custom });
  return { molecules, units, complexes: [complex] };
}

/** `SeqSets.build` or `withChanges` (a set's lineage kept), its issues under `svs.*` codes; `evict`: dynamic nodes. */
export function assembleSet(
  prev: SeqSet | null,
  upsert: SetPart,
  remove: string[],
  evict?: string[],
): { set: SeqSet; issues: Issue[] } {
  const r = prev ? SeqSets.withChanges(prev, { upsert, remove, evict }) : SeqSets.build({ molecules: [], ...upsert });
  return { set: r.set, issues: r.issues.map(svsIssue) };
}

/** Query nodes, each after the query nodes its expression names; a cycle (`svs.ref.cycle`) left in document order. */
function queryOrder(nodes: DocNode[], index: DocIndex): DocNode[] {
  const refs = new Map(nodes.map((d) => [d.ref as string, queryRefs(node(index, d).params.query as SeqQLQuery)]));
  const out: DocNode[] = [];
  const done = new Set<string>();
  for (let left = nodes; left.length > 0; ) {
    const ready = left.filter((d) => (refs.get(d.ref as string) ?? []).every((r) => done.has(r) || !refs.has(r)));
    for (const d of ready.length ? ready : left.slice(0, 1)) {
      out.push(d);
      done.add(d.ref as string);
    }
    left = left.filter((d) => !done.has(d.ref as string));
  }
  return out;
}

/**
 * A `computed` node's annotation by its registered provider; null when unavailable (reported, or quietly when it reads
 * a node in `ctx.unavailable`, or by the provider when it returns null): the caller marks it unavailable.
 */
export async function computeProvider(
  dn: DocNode,
  doc: DocIndex,
  ctx: StageContext,
  loaders: Loaders,
): Promise<AnnotationInput | null> {
  const n = node(doc, dn);
  const { provider, ...params } = n.params;
  const entry = ctx.registries.get("providers", provider as string);
  if (!entry) {
    unavailable(ctx, n, `provider "${provider}" is not registered`);
    return null;
  }
  if (blocked(ctx, n)) return null;
  // A statistic a streamed alignment's source serves (`impliedSource`): a streamed annotation over its columns.
  const implied = impliedSource(n, doc);
  if (implied) return servedStatistic(n, String(params.alignment), await loadSource(implied, ctx, loaders, doc));
  try {
    const a = await entry.compute(at(ctx, n.path), isNamespaced(provider) ? params.options : params);
    return a && { ...a, key: n.ref, custom: n.node.custom ?? a.custom };
  } catch (e) {
    ctx.signal.throwIfAborted();
    unavailable(ctx, n, `provider "${provider}" failed: ${messageOf(e)}`);
    return null;
  }
}

/** A statistic served by a streamed source (`<ref>#stats`) as a streamed `values` annotation over the alignment's columns. */
function servedStatistic(n: IndexedNode, alignment: string, sr: SourceResult): AnnotationInput | null {
  if (sr.kind !== "streamed") return null;
  const props = {
    [StreamProp]: { ...sr, systems: new Map([[alignment, 0]]), labels: new Map() } satisfies StreamedData,
  };
  return clean({
    key: n.ref,
    geometry: "values" as const,
    columns: {},
    streamed: [alignment],
    custom: n.node.custom,
    props,
  });
}

const LISTS = ["molecules", "units", "complexes", "alignments", "mappings", "annotations"] as const;
function merge(into: SetPart, part: SetPart) {
  for (const k of LISTS) if (part[k]?.length) into[k] = [...(into[k] ?? []), ...(part[k] as never[])] as never;
}

/**
 * Model issue paths (`/upsert/annotations/2/columns/beg/5`) as paths into the document: `assembleSet`'s issues for
 * `input` as `compile` reports them, for a driver that assembles itself.
 */
export function remapIssues(issues: Issue[], input: Partial<SetInput>, doc: DocIndex): Issue[] {
  return issues.map((i) => {
    const m = /^(?:\/upsert)?\/(\w+)\/(\d+)(.*)$/.exec(i.path);
    const key = m ? (input[m[1] as (typeof LISTS)[number]]?.[Number(m[2])] as { key?: string })?.key : undefined;
    const n = key ? doc.nodes.get(head(key)) : undefined;
    if (!n || !m) return i;
    const value = m[3].startsWith("/columns") && n.params.value != null;
    return { ...i, path: value ? `${n.path}/params/value${m[3].slice(8)}` : n.path };
  });
}

/**
 * Before assembling: a unit whose `sequence` names a ligand or `ligand` a polymer (known for derived refs only now,
 * `svs.ref.wrong_kind`) is left out of `input.units`. `input` holds the units' molecules.
 */
export function checkUnitKinds(input: Partial<SetInput>, doc: DocIndex, ctx: StageContext): void {
  const kinds = new Map((input.molecules ?? []).map((m) => [m.key, m.kind ?? "polymer"]));
  input.units = input.units?.filter((u) => {
    const n = doc.nodes.get(u.key);
    const param = n?.params.ligand != null ? "ligand" : "sequence";
    const got = kinds.get(u.molecule);
    if (!n || !got || (got === "ligand") === (param === "ligand")) return true;
    const message = `"${u.molecule}" is a ${got === "ligand" ? "ligand" : "polymer"}, not a ${param}`;
    report(ctx, "svs.ref.wrong_kind", pathOf(n, param), message);
    return false;
  });
}

/**
 * Validates, normalizes and loads a document into a set, stage by stage: a driver over the stage functions that keeps
 * the sources' results and the unavailable refs and passes them in.
 */
export async function compile(
  doc: SeqViewSpecDoc,
  reg: SvsRegistries,
  loaders: Loaders,
  signal: AbortSignal = new AbortController().signal,
): Promise<{ set: SeqSet; issues: Issue[] }> {
  const issues = validate(doc, reg, { base: loaders.base });
  if (issues.some((i) => i.severity === "error")) return { set: SeqSets.build({ molecules: [] }).set, issues };
  const index = indexDoc(doc, reg);
  const [sources, gone] = [new Map<string, SourceResult>(), new Set<string>()];
  const ctx = {
    signal,
    registries: reg,
    set: null as SeqSet | null,
    sources,
    unavailable: gone,
    source: (ref: string) => sourceParse(index, sources, ref),
    report: (i: Issue) => issues.push(i),
  };
  const of = (...kinds: string[]) => index.order.filter((n) => kinds.includes(n.kind)).map((n) => n.node);
  const stage = async (
    nodes: DocNode[],
    run: (d: DocNode) => Promise<LoadResult>,
    check?: (input: SetPart) => void,
  ) => {
    const input: SetPart = {};
    for (const d of nodes) {
      signal.throwIfAborted();
      const r = await run(d);
      if (r.unavailable) gone.add(d.ref ?? d.kind);
      merge(input, r);
    }
    check?.(input);
    if (!LISTS.some((k) => input[k]?.length) && ctx.set) return;
    const r = assembleSet(ctx.set, input, []);
    for (const i of remapIssues(r.issues, input, index)) issues.push(i);
    ctx.set = r.set;
  };
  const load = (d: DocNode) => loadNode(d, index, ctx, loaders);
  // The sources, and the ones the annotations' `uri` imply (read once, by whoever reads them).
  for (const n of index.order) {
    const implied = n.kind === "annotation" ? impliedSource(n, index) : null;
    if (n.kind !== "source" && !implied) continue;
    signal.throwIfAborted();
    const r: LoadResult = implied ? { source: await loadSource(implied, ctx, loaders, index) } : await load(n.node);
    const ref = implied ? implied.ref : n.ref;
    if (r.source) sources.set(ref, r.source);
    if (r.unavailable) gone.add(ref);
  }
  // Availability is decided bottom-up: complexes from providers, molecules and alignments, then the units that name
  // them; complexes by value last, once their units' availability is known.
  // A complex from a provider or a letters stream (`source_ref` without a provider) gives its molecules and units.
  // An alignment computed by a provider reads the set, so it is loaded once its inputs are in it (as `align` mappings).
  const given = (d: DocNode) => node(index, d).kind !== "alignment" || node(index, d).params.provider == null;
  const byProvider = (d: DocNode) => node(index, d).params.provider != null || node(index, d).params.source_ref != null;
  await stage(
    [...of("complex").filter(byProvider), ...of("sequence", "ligand", "alignment").filter(given), ...of("unit")],
    (d) => (d.kind === "complex" ? buildComplex(d, index, ctx) : load(d)),
    (input) => {
      // A unit left out for naming the wrong kind is unavailable, so a complex by value leaves it out too.
      const before = input.units ?? [];
      checkUnitKinds(input, index, ctx);
      const kept = new Set(input.units);
      for (const u of before) if (!kept.has(u)) gone.add(u.key);
    },
  );
  await stage(
    of("complex").filter((d) => !byProvider(d)),
    (d) => buildComplex(d, index, ctx),
  );
  // One stage each, in document order: a later one may align the rows of an earlier one.
  for (const d of of("alignment").filter((d) => !given(d))) await stage([d], load);
  await stage(of("annotation", "mapping"), load);
  const compute = async (d: DocNode): Promise<LoadResult> => {
    const a = await computeProvider(d, index, ctx, loaders);
    return a ? { annotations: [a] } : UNAVAILABLE;
  };
  const isQuery = (d: DocNode) => node(index, d).params.provider === "query";
  // `redundancy` orders by a key table, a computed node's result, so it runs after the others, in document order.
  const isLate = (d: DocNode) => node(index, d).params.provider === "redundancy";
  await stage(
    of("computed").filter((d) => !isQuery(d) && !isLate(d)),
    compute,
  );
  for (const d of of("computed").filter(isLate)) await stage([d], compute);
  // Queries last, each in a stage of its own, after the queries it names (seqspec:seqql#in-documents).
  for (const d of queryOrder(of("computed").filter(isQuery), index)) await stage([d], compute);
  return { set: ctx.set as SeqSet, issues };
}
