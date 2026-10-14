/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1b-a: dynamic nodes (seqspec:sources#dynamic-nodes). A manifest's or a patch's `nodes` become set inputs: each node
// is checked as a document node of its kind is (seqspec:format#validation), its ref derived (`<source>/<key>`), its
// refs naming the document's nodes, the set's keys, dynamic nodes held or the same patch's (any other is
// `svs.stream.ref`, and the node is left out), then read by the same `loadNode` as a document node. A node provider's
// result passes the same checks (seqspec:model#providers).

import type { AnnotationInput, CustomData, KeyResolver, KeyTarget, SeqSet } from "@molstar/seq-core/model";
import type { Loaders, StageContext, StreamNode } from "../contract.js";
import type { DocNode } from "../doc.js";
import { normalizeNodeParams } from "../normalize.js";
import { isRefSyntax } from "../schema/field-schema.js";
import type { IndexedNode } from "../schema/tree-schema.js";
import { SvsTreeSchema } from "../tables/svs-tree.js";
import { checkStreamNode } from "../validate.js";
import { type DocIndex, type LoadResult, loadNode, type StreamedData, StreamProp } from "./compile.js";
import { issueOf } from "./issues.js";

/** The kinds a patch may add; a `source` is a nested streamed source (ST1b-b). */
const KINDS: ReadonlySet<string> = new Set(["sequence", "ligand", "unit", "mapping", "annotation", "source"]);

/**
 * An item coordinate system (ST7, seqspec:sources#dynamic-nodes): a streamed item's own positions, a viewer's dynamic
 * node of its source made from the item's `map` (never a provider's: `#` kinds are refused from patches), read as a
 * molecule with a length and no letters and its own unit, marked `ItemSystemProp` so lists of the document's sequences
 * leave it out. Params: `length` (positive), `alphabet`, `label`.
 */
export const ItemSystemKind = "#item_system";
export const ItemSystemProp = "svs:item_system";

/** A map with more entries over it (a document's nodes and held nested sources, their results): `get` and `has` only. */
export class MapOver<V> extends Map<string, V> {
  private readonly base: ReadonlyMap<string, V>;
  private readonly more: ReadonlyMap<string, V>;
  constructor(base: ReadonlyMap<string, V>, more: ReadonlyMap<string, V>) {
    super();
    this.base = base;
    this.more = more;
  }
  get(ref: string) {
    return this.more.get(ref) ?? this.base.get(ref);
  }
  has(ref: string) {
    return this.more.has(ref) || this.base.has(ref);
  }
}

/**
 * A streamed source's level (the document's is 1, a source it nests 2, ...) and the `max_depth` of the document's
 * source at the top of its chain; nested sources are in `doc` with their parent source as `parent`.
 */
export function sourceLevel(doc: DocIndex, ref: string): { level: number; max: number } {
  let n = doc.nodes.get(ref);
  let level = 1;
  for (let k = 0; n && k < 64; k++) {
    const up = doc.nodes.get(n.parent ?? "");
    if (up?.kind !== "source") break;
    level++;
    n = up;
  }
  return { level, max: typeof n?.params.max_depth === "number" ? n.params.max_depth : 4 };
}

/** A dynamic node's ref: derived from its source's. */
export const nodeRef = (source: string, key: string) => `${source}/${key}`;

/**
 * One dynamic node of `source` as set inputs; null when it is left out (its issues reported through `ctx`). `known`:
 * whether a ref names something a dynamic node may name (the document's nodes, the set's keys, held dynamic nodes and
 * the same patch's). Issues are at `/nodes/<key>` of the patch that added it, naming its source. `doc`: the document's
 * index with the held nested sources (their parent source as `parent`), so a node reads a streamed source by ref (an
 * annotation's `source_ref`, its result in `ctx.sources`); a `source` node (nested) is read as a streamed source, its
 * manifest asked for, within its chain's `max_depth` (deeper: `svs.stream.depth`, not requested). `viewers`: the viewer's own node (an item coordinate system, ST7); a `#` kind
 * from anyone else is `svs.node.kind`.
 */
export async function loadStreamNode(
  source: string,
  sn: StreamNode,
  known: (ref: string) => boolean,
  ctx: StageContext,
  loaders: Loaders,
  doc?: DocIndex,
  viewers = false,
): Promise<LoadResult | null> {
  const ref = nodeRef(source, String(sn?.key));
  const path = `/nodes/${String(sn?.key)}`;
  const say = (code: string, at: string, message: string) =>
    ctx.report(issueOf(ctx.registries, code, at, `${source}: ${message}`));
  if (typeof sn?.key !== "string" || !isRefSyntax(ref, true)) {
    say("svs.ref.syntax", path, `a dynamic node's key is a ref suffix ([A-Za-z0-9_.:@-]+, joined by /): ${sn?.key}`);
    return null;
  }
  if (!viewers && String(sn?.kind).startsWith("#")) {
    say("svs.node.kind", `${path}/kind`, `"${sn.kind}" is the viewer's own kind: a provider's nodes may not have it`);
    return null;
  }
  if ((sn.kind as string) === ItemSystemKind) {
    const p = sn.params ?? {};
    const length = typeof p.length === "number" && p.length > 0 ? Math.floor(p.length) : 0;
    if (length === 0) {
      say("svs.params.missing", `${path}/params/length`, "an item coordinate system needs a positive length");
      return null;
    }
    const alphabet = (typeof p.alphabet === "string" ? p.alphabet : "generic") as "generic";
    const label = typeof p.label === "string" && p.label !== "" ? p.label : undefined;
    const props = { [ItemSystemProp]: true };
    return { molecules: [{ key: ref, kind: "polymer", alphabet, lazy: { length }, label, ownUnit: true, props }] };
  }
  if (!KINDS.has(sn.kind)) {
    say("svs.node.kind", `${path}/kind`, `"${sn.kind}" is no dynamic node's kind (${[...KINDS].join(", ")}, source)`);
    return null;
  }
  const raw: DocNode = { kind: sn.kind, ref, params: sn.params ?? {} };
  const checks = checkStreamNode(
    { node: raw, ref, kind: sn.kind, params: raw.params ?? {}, parent: null, path },
    ctx.registries,
    loaders.base,
  );
  for (const i of checks.issues) ctx.report({ ...i, message: `${source}: ${i.message}` });
  if (checks.issues.some((i) => i.severity === "error")) return null;
  for (const u of checks.refs) {
    if (known(u.ref)) continue;
    say("svs.stream.ref", u.path, `"${u.ref}" names no node of the document, held or added by the same patch`);
    return null;
  }
  const norm = normalizeNodeParams(raw, SvsTreeSchema, ctx.registries);
  const nested = sn.kind === "source";
  const n: IndexedNode = {
    node: { ...raw, ...norm },
    ref,
    kind: sn.kind,
    params: norm.params,
    parent: nested ? source : null,
    path,
  };
  const own = new Map([[ref, n]]);
  const index: DocIndex = { nodes: doc ? new MapOver(doc.nodes, own) : own, order: [n] };
  if (nested) {
    // A nested source streams, its data kind given (nothing in the patch is read to infer it).
    const p = n.params;
    if ((p.stream == null && p.protocol == null) || p.data_kind == null) {
      say("svs.params.missing", `${path}/params`, "a nested source streams: give stream or protocol, and data_kind");
      return null;
    }
    const { level, max } = sourceLevel(index, ref);
    if (level > max) {
      say("svs.stream.depth", path, `"${ref}" is level ${level} of streamed sources, deeper than max_depth ${max}`);
      return null;
    }
  }
  const part = await loadNode(n.node, index, ctx, loaders);
  return part.unavailable ? null : part;
}

/** Every key a part puts in the set (a molecule's own unit too): what an eviction takes. */
export function keysOfPart(part: Partial<LoadResult>, out: string[]) {
  for (const m of part.molecules ?? []) out.push(m.key, ...(m.ownUnit ? [`${m.key}/unit`] : []));
  for (const k of ["units", "complexes", "alignments", "mappings", "annotations"] as const)
    for (const x of part[k] ?? []) out.push(x.key);
}

/** The refs of the streamed sources a set's data reads (its streamed annotations' and lazy alignments'). */
const streamed = new WeakMap<SeqSet, ReadonlySet<string>>();
export function streamedSources(set: SeqSet): ReadonlySet<string> {
  let out = streamed.get(set);
  if (out) return out;
  const refs = new Set<string>();
  const see = (x: { readonly custom: CustomData } | null) => {
    const d = x?.custom.props.get(StreamProp) as StreamedData | undefined;
    if (d) refs.add(d.source.ref);
  };
  for (const a of set.annotations) see(a);
  for (const a of set.alignments) see(a);
  out = refs;
  streamed.set(set, out);
  return out;
}

/** Whether `ref` lies under one of `sources` (`<source>/<key>`): a dynamic node's ref. */
export function underSource(sources: ReadonlySet<string>, ref: string): boolean {
  for (let k = ref.indexOf("/"); k > 0; k = ref.indexOf("/", k + 1)) if (sources.has(ref.slice(0, k))) return true;
  return false;
}

/**
 * A data-mode resolver that also accepts refs under the streamed sources its data reads: dynamic
 * nodes, valid before they arrive. `inputs`: the annotations given as build inputs.
 */
export function withDynamic(keys: KeyResolver, sources: ReadonlySet<string>): KeyResolver {
  if (sources.size === 0) return keys;
  const t: KeyTarget = { kind: -1, id: -1, row: -1 };
  return {
    byName: (key, out) => keys.byName(key, out),
    molecules: keys.molecules,
    dynamic: (key) => !keys.byName(key, t) && underSource(sources, key),
  };
}

/** The streamed sources build inputs read (their annotations' `svs:stream`). */
export function inputSources(annotations: readonly AnnotationInput[] | undefined): ReadonlySet<string> {
  const refs = new Set<string>();
  for (const a of annotations ?? []) {
    const d = a.props?.[StreamProp] as StreamedData | undefined;
    if (d) refs.add(d.source.ref);
  }
  return refs;
}
