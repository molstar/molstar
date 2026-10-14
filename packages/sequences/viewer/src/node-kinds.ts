/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in node kinds (seqspec:viewer#extension-points; seqspec:viewer#updates-and-status): each wraps one stage
// function of the format package. Edges are the format's references (seqspec:format#what-a-change-affects): a data node
// reads a source's result and other data nodes' availability; assemblies read their parts; per-view stages read the
// section's axis, never its whole result. Keys
// hold every input a stage reads (VC1): the content token of the blocks' molecules and alignments for prepare and theme
// nodes whose entries read residues, the mappings and `of` polymers a layer is placed through.

import {
  assembleSet,
  buildComplex,
  checkUnitKinds,
  childrenOf,
  composeView,
  computeProvider,
  type DocIndex,
  type ForEach,
  forEachOf,
  groupMembers,
  type IndexedNode,
  ItemSystemKind,
  instantiate,
  isStreamedSource,
  keysOfPart,
  type LoadResult,
  letterMapParam,
  loadNode,
  loadStreamNode,
  MapOver,
  nodeRef,
  nodeRefs,
  type ResolvedLayer,
  type ResolvedSection,
  remapIssues,
  resolveLayer,
  resolveSection,
  type SchemeSubject,
  type SourceResult,
  type StageContext,
  type StreamedData,
  type StreamNode,
  StreamProp,
  type SvsRegistries,
  type ThemeRef,
  templateNodes,
} from '@molstar/seq-core/format';
import {
  type Annotation,
  type AnnotationInput,
  Issue,
  type KeyTarget,
  type Molecule,
  type SeqSet,
  type SetInput,
} from '@molstar/seq-core/model';
import { pickedTheme, schemeEntries, subjectOf, withRange } from './color-schemes.js';
import type { InputCollector, NodeKind, PreparedLayer, RunContext, Stage, ViewModel } from './graph.js';
import type { LayerInput, LayerRenderer } from './registry.js';
import { createTheme, type Theme, themeContext } from './themes.js';

/** An alignment a node provider computes from the set (`realign`): assembled with the mappings, after the base set. */
const computedAlignment = (n: IndexedNode) => n.kind === 'alignment' && n.params.provider != null;
/** A node that reads the set before it is assembled: an `align` or provider mapping, a computed alignment. */
const readsSet = (n: IndexedNode) =>
  computedAlignment(n) || (n.kind === 'mapping' && (n.params.method === 'align' || n.params.method === 'provider'));
const BASE = ['complex', 'sequence', 'ligand', 'unit', 'alignment'];
const head = (ref: string) => ref.split('/')[0];

/**
 * The heads of a node's refs that name another document node: core's `nodeRefs`, every typed ref its params hold
 * (seqspec:format#nodes-and-refs), the one source of the edges and keys below.
 */
function refsOf(n: IndexedNode, doc: DocIndex, svs: SvsRegistries): string[] {
  const out = new Set<string>();
  for (const u of nodeRefs(n, svs)) out.add(head(u.ref));
  return [...out].filter((r) => doc.nodes.has(r) && r !== n.ref);
}
function edges(n: IndexedNode, doc: DocIndex, out: InputCollector, svs: SvsRegistries) {
  for (const r of refsOf(n, doc, svs)) (doc.nodes.get(r)?.kind === 'source' ? out.result : out.status)(r);
}
/**
 * The computed alignments a node reading the set names, and the ones those name: `compile` assembles each computed
 * alignment after the ones before it in the document, and the mappings after them all, so these are in the set it reads.
 */
function computedNamed(n: IndexedNode, doc: DocIndex, svs: SvsRegistries): string[] {
  const out: string[] = [];
  // Its place in the document, only where it is read (a computed alignment names only the ones before it).
  const at = computedAlignment(n) ? doc.order.indexOf(n) : -1;
  const visit = (m: IndexedNode) => {
    for (const r of refsOf(m, doc, svs)) {
      const d = doc.nodes.get(r);
      if (!d || !computedAlignment(d) || out.includes(r)) continue;
      if (computedAlignment(n) && doc.order.indexOf(d) > at) continue;
      out.push(r);
      visit(d);
    }
  };
  if (readsSet(n)) visit(n);
  return out;
}
/** The stage a node reading the set runs in: the base set with the computed alignments it names (`computedNamed`). */
function withComputed(rc: RunContext, n: IndexedNode): StageContext {
  const part: Part = {};
  for (const r of computedNamed(n, rc.doc, rc.registry.svs)) add(part, rc.result(r));
  if (empty(part) || !rc.stage.set) return rc.stage;
  return { ...rc.stage, set: assembleSet(rc.stage.set, part, []).set };
}

/** What a ref names in a set: the object (or table) a node reading it depends on. */
const t: KeyTarget = { kind: -1, id: -1, row: -1 };
/** The set's table by `KeyKind` (molecule, unit, alignment, row: its alignment, mapping, annotation, complex). */
const ENTITY = ['molecules', 'units', 'alignments', 'alignments', 'mappings', 'annotations', 'complexes'] as const;
function entity(set: SeqSet | null, ref: string): unknown {
  if (!set?.byName(ref, t)) return undefined;
  const list = set[ENTITY[t.kind]];
  return ENTITY[t.kind] === 'units' ? list : (list as readonly unknown[])[t.id];
}

function load(name: string, stage: Stage = 'load'): NodeKind<LoadResult> {
  return {
    name,
    stage,
    inputs(_ref, n, doc, out, svs) {
      if (!n) return;
      edges(n, doc, out, svs);
      for (const c of childrenOf(doc, n)) if (c.kind === 'numbering' || c.kind === 'tree') edges(c, doc, out, svs);
      // Targets read from a file may name any data node; an aligner reads its endpoints in the set.
      if (n.kind === 'annotation' && n.params.value == null)
        for (const d of doc.order) if (BASE.includes(d.kind)) out.status(d.ref);
      if (readsSet(n)) out.after('#assemble.base');
      for (const r of computedNamed(n, doc, svs)) out.result(r);
      // A complex by value leaves out its unavailable units: it runs once their availability is known.
      if (n.kind === 'complex' && n.params.provider == null && n.params.source_ref == null)
        for (const c of childrenOf(doc, n)) if (c.kind === 'unit') out.status(c.ref);
    },
    key(rc, out) {
      const n = rc.node as IndexedNode;
      for (const c of n.node.children ?? []) out.push(rc.token(c.ref ?? c.kind));
      const parent = rc.doc.nodes.get(n.parent ?? '');
      if (parent?.kind === 'complex') out.push(rc.token(parent.ref));
      if (typeof n.params.from_source === 'string') out.push(rc.token(head(n.params.from_source)));
      // A lazy alignment fetches the rows the document's annotations and mappings name (`namedRows`).
      const src = n.kind === 'alignment' ? rc.doc.nodes.get(String(n.params.source_ref)) : undefined;
      if (src && (src.params.protocol != null || src.params.stream != null))
        for (const d of rc.doc.order) if (d.kind === 'annotation' || d.kind === 'mapping') out.push(rc.token(d.ref));
      if (n.params.method === 'align' || n.params.method === 'provider')
        out.push(entity(rc.stage.set, String(n.params.from)), entity(rc.stage.set, String(n.params.to)));
      // A computed alignment reads its inputs' residues in the set.
      if (computedAlignment(n))
        for (const r of typeof n.params.inputs === 'string' ? [n.params.inputs] : (n.params.inputs as string[]))
          out.push(entity(rc.stage.set, r));
    },
    run: (rc) => {
      const n = rc.node as IndexedNode;
      if (stage === 'complex') return buildComplex(n.node, rc.doc, rc.stage);
      return loadNode(n.node, rc.doc, readsSet(n) ? withComputed(rc, n) : rc.stage, rc.loaders);
    },
  };
}

/**
 * The computed nodes a computed node names: its inputs, since the set it reads holds the data nodes, not the computed
 * ones. A query names any computed node, as the load runs every other provider before the queries and the queries in
 * order (`queryOrder`): queries (`(count cysteines)`) and other providers' answers. A `redundancy` names the key tables
 * its `order` reads (a `row_stats`).
 */
function answersNamed(n: IndexedNode | null, doc: DocIndex, svs: SvsRegistries): string[] {
  if (n?.params.provider === 'redundancy')
    return keyTables(n.params).filter((r) => r !== n.ref && doc.nodes.get(r)?.kind === 'computed');
  if (n?.params.provider !== 'query') return [];
  return refsOf(n, doc, svs).filter((r) => r !== n.ref && doc.nodes.get(r)?.kind === 'computed');
}
/** The key tables an `order` and a `filter` of a node's params read (keyed sets, seqspec:format#layer-types). */
function keyTables(p: Readonly<Record<string, unknown>>): string[] {
  const out = new Set<string>();
  for (const k of [p.order ?? [], p.filter ?? []].flat() as unknown[])
    if (k && typeof k === 'object' && typeof (k as { data?: unknown }).data === 'string')
      out.add(head((k as { data: string }).data));
  return [...out];
}

const computed: NodeKind<AnnotationInput | null> = {
  name: 'computed',
  stage: 'compute',
  inputs(_ref, n, doc, out, svs) {
    out.after('#assemble.mapped');
    if (n) edges(n, doc, out, svs);
    for (const r of answersNamed(n, doc, svs)) out.result(r);
  },
  // The provider reads the set: what its refs name there (a query's names too), not the set's identity; a query also
  // the answers of the computed nodes it names.
  key(rc, out) {
    const svs = rc.registry.svs;
    for (const r of refsOf(rc.node as IndexedNode, rc.doc, svs)) out.push(entity(rc.stage.set, r));
    for (const r of answersNamed(rc.node, rc.doc, svs)) out.push(rc.result(r));
  },
  run: (rc) => {
    // A query that names a computed node reads the set with that node's answer in it, as `compile` does.
    const named = answersNamed(rc.node, rc.doc, rc.registry.svs).flatMap((r) =>
      rc.result(r) ? [rc.result(r) as AnnotationInput] : [],
    );
    const stage =
      named.length && rc.stage.set
        ? { ...rc.stage, set: assembleSet(rc.stage.set, { annotations: named }, []).set }
        : rc.stage;
    return computeProvider((rc.node as IndexedNode).node, rc.doc, stage, rc.loaders);
  },
};

// --- assemblies -------------------------------------------------------------------------------------------------

/** An assembly's result: its set and the parts in it, by the node that gave each. */
interface Assembly {
  readonly set: SeqSet;
  readonly parts: ReadonlyMap<string, unknown>;
  /** The document's mappings over a streamed source's dynamic nodes (`overDynamic`), held back for `#assemble.stream`. */
  readonly deferred?: ReadonlyMap<string, unknown>;
}
type Part = Partial<SetInput>;
const LISTS = ['molecules', 'units', 'complexes', 'alignments', 'mappings', 'annotations'] as const;
const asPart = (v: unknown): Part =>
  v && typeof v === 'object' ? ('geometry' in v ? { annotations: [v as AnnotationInput] } : (v as Part)) : {};
function add(into: Part, v: unknown, keys?: Set<string>) {
  const p = asPart(v);
  for (const k of LISTS) {
    const list = p[k] as { key: string }[] | undefined;
    if (!list?.length) continue;
    if (keys) for (const x of list) keys.add(x.key);
    else into[k] = [...(into[k] ?? []), ...list] as never;
  }
}
const empty = (p: Part) => LISTS.every((k) => !p[k]?.length);
/** Parts without the units `checkUnits` left out, by part and the keys left out (kept while those stay the same). */
const withoutUnits = new WeakMap<object, { keys: string; part: Part }>();
/**
 * `compile`'s check before the base set (core's `checkUnitKinds`): a unit naming a molecule of the wrong kind is reported
 * and left out of the parts, and out of its complex by value (a complex left with none is left out).
 */
function checkUnits(parts: Map<string, unknown>, rc: RunContext) {
  const all: Part = {};
  for (const v of parts.values()) add(all, v);
  const units = all.units ?? [];
  checkUnitKinds(all, rc.doc, rc.stage);
  if ((all.units ?? []).length === units.length) return;
  const kept = new Set(all.units);
  const bad = new Set(units.flatMap((u) => (kept.has(u) ? [] : [u.key])));
  const keys = [...bad].join(' ');
  for (const [r, v] of parts) {
    const p = asPart(v);
    if (!p.units?.some((u) => bad.has(u.key)) && !p.complexes?.some((c) => c.units.some((u) => bad.has(u)))) continue;
    let w = withoutUnits.get(v as object);
    if (w?.keys !== keys) {
      const complexes = p.complexes?.flatMap((c) => {
        const left = c.units.filter((u) => !bad.has(u));
        return left.length ? [{ ...c, units: left }] : [];
      });
      w = { keys, part: { ...p, units: p.units?.filter((u) => !bad.has(u.key)), complexes } };
      withoutUnits.set(v as object, w);
    }
    parts.set(r, w.part);
  }
}
/**
 * A document mapping naming a streamed source's dynamic node (an item coordinate system, ST7): not in the document's
 * set (its end is not there), upserted by `#assemble.stream` while both ends are held.
 */
function overDynamic(doc: DocIndex, v: unknown): boolean {
  const p = asPart(v);
  for (const m of p.mappings ?? []) if (underStream(doc, m.from) || underStream(doc, m.to)) return true;
  for (const a of p.annotations ?? []) if (typeof a.target === 'string' && underStream(doc, a.target)) return true;
  return false;
}
/** A ref naming a streamed source's dynamic node (`<source>/<key>`). */
function underStream(doc: DocIndex, r: unknown): boolean {
  if (typeof r !== 'string') return false;
  const h = r.slice(0, Math.max(0, r.indexOf('/')));
  return h !== '' && isStreamedSource(doc.nodes.get(h));
}
/** An annotation whose per-row `target` or `b_target` names a dynamic node: not supported yet (its rows say so). */
function rowsOverDynamic(doc: DocIndex, v: unknown): string | null {
  for (const a of asPart(v).annotations ?? [])
    for (const c of ['target', 'b_target']) {
      const col = a.columns?.[c] as ArrayLike<unknown> | undefined;
      for (let i = 0; col && i < col.length; i++) if (underStream(doc, col[i])) return String(col[i]);
    }
  return null;
}
/**
 * What names a streamed item (an item coordinate system is read only then, ST7): the document's refs (section blocks,
 * mapping ends, annotation targets, layers' data), a `/unit` ref naming its item; and the sources whose every item a
 * template asks for (`for_each: { kind: "item_system" }`).
 */
interface ItemNames {
  readonly refs: ReadonlySet<string>;
  readonly sources: ReadonlySet<string>;
}
const itemNames = new WeakMap<DocIndex, ItemNames>();
function itemNamesOf(doc: DocIndex, svs: SvsRegistries): ItemNames {
  let n = itemNames.get(doc);
  if (n) return n;
  const sources = new Set<string>();
  for (const d of doc.order) {
    const fe = forEachOf(d);
    if (fe?.kind === 'item_system') sources.add(fe.source);
  }
  n = { refs: namedUnder(doc, svs), sources };
  itemNames.set(doc, n);
  return n;
}
/**
 * The refs under a source the document's nodes name (a `/` in them; `/unit` dropped), read from core's `nodeRefs`
 * as the graph's edges are: what item systems are made for, and what nested sources are asked for. Once per index.
 */
const under = new WeakMap<DocIndex, ReadonlySet<string>>();
function namedUnder(doc: DocIndex, svs: SvsRegistries): ReadonlySet<string> {
  let refs = under.get(doc);
  if (refs) return refs;
  const out = new Set<string>();
  for (const d of doc.order)
    for (const u of nodeRefs(d, svs))
      if (u.ref.includes('/')) out.add(u.ref.endsWith('/unit') ? u.ref.slice(0, -5) : u.ref);
  refs = out;
  under.set(doc, refs);
  return refs;
}

/**
 * `assembleSet` over the parts whose identity changed. Its first input, when an assembly, seeds it: the first run
 * continues that set's lineage with its own parts (as `compile` does); later runs apply their own changes, the seed's
 * included, to their previous set, so the three sets share one lineage and untouched tables.
 */
function assembly(name: string, seed: string | null, groups: string[][], computed = false): NodeKind<Assembly> {
  const units = groups.some((g) => g.includes('unit'));
  return {
    name,
    stage: 'assemble',
    inputs(_ref, _n, doc, out) {
      if (seed) out.result(seed);
      // `computed`: the alignments a node provider computes belong to this assembly, and the given ones to the base's.
      for (const g of groups)
        for (const d of doc.order)
          if (g.includes(d.kind) && (d.kind !== 'alignment' || computedAlignment(d) === computed)) out.result(d.ref);
    },
    run(rc) {
      const prev = rc.prev as Assembly | undefined;
      const first = rc.result(rc.inputs[0]) as Assembly | undefined;
      const base = first && typeof first === 'object' && 'parts' in first ? first : null;
      const own = base ? rc.inputs.slice(1) : rc.inputs;
      const parts = new Map(base?.parts);
      const deferred = new Map(base?.deferred);
      for (const r of own) {
        const v = rc.result(r);
        const row = rowsOverDynamic(rc.doc, v);
        if (row !== null)
          rc.stage.report(
            Issue.warning(
              'svs.provider.unavailable',
              `${rc.doc.nodes.get(r)?.path ?? ''}/params`,
              `"${r}": a row's target "${row}" names a streamed item: per-row targets cannot wait for one yet; give the annotation one target`,
            ),
          );
        if (overDynamic(rc.doc, v)) deferred.set(r, v);
        else parts.set(r, v);
      }
      if (units) checkUnits(parts, rc);
      // The model's issues name the document node that gave the part, as `compile` reports them; a change of the
      // seed's that this one applies again is the seed's to report, once.
      const said = (i: Issue) => `${i.code} ${i.path} ${i.message}`;
      const seeds = new Set(base ? (rc.status(rc.inputs[0])?.issues ?? []).map(said) : []);
      const done = (r: { set: SeqSet; issues: Issue[] }, input: Part) => {
        for (const i of remapIssues(r.issues, input, rc.doc)) if (!seeds.has(said(i))) rc.stage.report(i);
        return { set: r.set, parts, deferred };
      };
      if (!prev || (base && base.set.lineage !== prev.set.lineage)) {
        const input: Part = {};
        for (const r of own) if (parts.has(r)) add(input, parts.get(r));
        if (base && empty(input)) return { set: base.set, parts, deferred };
        return done(assembleSet(base?.set ?? null, input, []), input);
      }
      const upsert: Part = {};
      const remove = diffParts(prev.parts, parts, prev.set, upsert);
      const same =
        prev.deferred?.size === deferred.size && [...deferred].every(([k, v]) => prev.deferred?.get(k) === v);
      if (empty(upsert) && remove.length === 0) return same ? prev : { ...prev, deferred };
      return done(assembleSet(prev.set, upsert, remove), upsert);
    },
  };
}

/**
 * The parts whose identity changed from `prev` to `next`, into `upsert`; the keys to remove (left a changed part and in
 * none of `next`'s, and held by `set`: a key that never built is not removed again).
 */
function diffParts(prev: ReadonlyMap<string, unknown>, next: ReadonlyMap<string, unknown>, set: SeqSet, upsert: Part) {
  const keep = new Set<string>();
  for (const [ref, v] of next) {
    if (prev.get(ref) !== v) add(upsert, v);
    add(upsert, v, keep);
  }
  const remove = new Set<string>();
  for (const [ref, v] of prev) if (next.get(ref) !== v) add({}, v, remove);
  for (const k of keep) remove.delete(k);
  return [...remove].filter((k) => set.byName(k, t));
}

/**
 * A held nested source (seqspec:sources#dynamic-nodes): its node and index node, and once requested its result (null:
 * left out).
 */
interface Nested {
  readonly source: string;
  readonly node: StreamNode;
  readonly index: IndexedNode;
  result: SourceResult | null | undefined;
}
type Held = { readonly source: string; readonly node: StreamNode };
/** `#assemble.stream`'s entry for a document mapping over dynamic nodes (`overDynamic`): its part is its own. */
const DEFERRED: StreamNode = { kind: 'mapping', key: '', params: {} };
/**
 * `#assemble.stream`'s result: the final set over `base` (`#assemble`'s), per dynamic node its node and part, the held
 * nested sources, and per source a version bumped when its nodes change (`own`), or those of any source under it (`tree`).
 */
export interface StreamAssembly extends Assembly {
  readonly base: Assembly;
  readonly nodes: ReadonlyMap<string, Held & { readonly part: Part | null; readonly issues?: readonly Issue[] }>;
  readonly nested: ReadonlyMap<string, Nested>;
  readonly own: ReadonlyMap<string, number>;
  readonly tree: ReadonlyMap<string, number>;
  /** The held-back document nodes' issues, by ref, while they are in. */
  readonly deferredIssues?: ReadonlyMap<string, readonly Issue[]>;
}
/** The manifests' nodes (known up front): the document's streamed sources', the set's, the nested. */
function manifestNodes(rc: RunContext, set: SeqSet, nested: ReadonlyMap<string, Nested>, out: Map<string, Held>) {
  const seen = new Set<string>();
  const from = (ref: string, nodes: readonly StreamNode[] | undefined) => {
    if (seen.has(ref)) return;
    seen.add(ref);
    for (const node of nodes ?? []) out.set(nodeRef(ref, String(node?.key)), { source: ref, node });
  };
  for (const [ref, sr] of rc.stage.sources) if (sr.kind === 'streamed') from(ref, sr.manifest.nodes);
  for (const a of set.annotations) {
    const d = a?.streamed ? (a.custom.props.get(StreamProp) as StreamedData | undefined) : undefined;
    if (d) from(d.source.ref, d.manifest.nodes);
  }
  for (const [ref, x] of nested) if (x.result?.kind === 'streamed') from(ref, x.result.manifest.nodes);
}

/** What a document asks of nested sources: the refs its layers and focus name, and its templates (once per index). */
const demands = new WeakMap<DocIndex, { refs: ReadonlySet<string>; templates: ForEach[] }>();
function demandsOf(doc: DocIndex, svs: SvsRegistries) {
  let d = demands.get(doc);
  if (d) return d;
  d = { refs: namedUnder(doc, svs), templates: [] };
  for (const n of doc.order) {
    const fe = forEachOf(n);
    if (fe) d.templates.push(fe);
  }
  demands.set(doc, d);
  return d;
}
/**
 * A nested source is requested only when something needs its data: a held node reads it, a document ref lies under it,
 * or a template selects its nodes (its source, or one above it within `depth`).
 */
function needed(
  ref: string,
  reads: ReadonlySet<string>,
  doc: DocIndex,
  nested: ReadonlyMap<string, Nested>,
  svs: SvsRegistries,
): boolean {
  if (reads.has(ref)) return true;
  const d = demandsOf(doc, svs);
  const under = `${ref}/`;
  for (const r of d.refs) if (r.startsWith(under)) return true;
  for (const fe of d.templates) {
    if (fe.source === ref) return true;
    if (!ref.startsWith(`${fe.source}/`)) continue;
    let level = 1;
    for (let s = nested.get(ref)?.source; s && s !== fe.source; s = nested.get(s)?.source) level++;
    if (fe.depth === null || level <= fe.depth) return true;
  }
  return false;
}

const NO_VERSIONS: ReadonlyMap<string, number> = new Map();
/**
 * `#assemble.stream`: the document's set with the dynamic nodes held, the manifests' and the stream manager's
 * (`rc.dynamic`), upserted and evicted as they come and go; each node checked and read once (its issues reported then).
 * Nested sources first: each requested (its manifest) when needed, its manifest's nodes then wanted too, to
 * `max_depth`. Without dynamic nodes it is `#assemble`'s result itself. A change of the document's set is applied to
 * its own previous set, so the final set keeps one lineage and every id, the dynamic nodes' included.
 */
const streamAssembly: NodeKind<Assembly> = {
  name: '#assemble.stream',
  stage: 'assemble',
  inputs(_ref, _n, doc, out) {
    out.result('#assemble');
    for (const d of doc.order) if (isStreamedSource(d)) out.result(d.ref);
  },
  key: (rc, out) => out.push(rc.dynamic?.version ?? 0),
  async run(rc) {
    const base = rc.result('#assemble') as Assembly;
    const p = rc.prev as Partial<StreamAssembly> | undefined;
    const was = p?.nodes && p.base ? (p as StreamAssembly) : null;
    const want = new Map<string, Held>();
    const nested = new Map<string, Nested>();
    // The document's index and sources with the nested ones over them: a node reads a nested source by ref.
    const index = new Map<string, IndexedNode>();
    const results = new Map<string, SourceResult>();
    const doc: DocIndex = { nodes: new MapOver(rc.doc.nodes, index), order: rc.doc.order };
    const stage: StageContext = Object.assign(Object.create(rc.stage), {
      sources: new MapOver(rc.stage.sources, results),
    });
    // A node's refs name the document's nodes, the set's keys, or dynamic nodes held or of the same patch (and theirs).
    const known = (r: string) =>
      want.has(r) ||
      want.has(r.slice(0, Math.max(0, r.lastIndexOf('/')))) ||
      rc.doc.nodes.has(r) ||
      base.set.byName(r, t);
    for (let round = 0; round < 64; round++) {
      want.clear();
      manifestNodes(rc, base.set, nested, want);
      const seen = nested.size;
      for (const [ref, e] of rc.dynamic?.entries ?? []) if (!want.has(ref)) want.set(ref, e);
      const reads = new Set<string>();
      for (const [ref, w] of want) {
        const r = w.node?.params?.source_ref;
        if (typeof r === 'string') reads.add(r);
        if (w.node?.kind !== 'source' || nested.get(ref)?.node === w.node) continue;
        const old = was?.nested.get(ref);
        const params = w.node.params ?? {};
        const node = { kind: 'source', ref, params };
        const index1: IndexedNode = { node, ref, kind: 'source', params, parent: w.source, path: '' };
        const x: Nested =
          old?.node === w.node ? old : { source: w.source, node: w.node, index: index1, result: undefined };
        nested.set(ref, x);
        index.set(ref, x.index);
        if (x.result) results.set(ref, x.result);
      }
      // The sources of a depth are requested together (their answers can take a while), not one after another.
      let more = false;
      const asking: Promise<void>[] = [];
      for (const [ref, x] of nested) {
        if (x.result !== undefined || !needed(ref, reads, rc.doc, nested, rc.registry.svs)) continue;
        asking.push(
          loadStreamNode(x.source, x.node, known, stage, rc.loaders, doc).then((r) => {
            x.result = r?.source ?? null;
            if (x.result) results.set(ref, x.result);
            more ||= x.result?.kind === 'streamed' && (x.result.manifest.nodes?.length ?? 0) > 0;
          }),
        );
      }
      await Promise.all(asking);
      // Again while a source with a manifest joined after its manifest's nodes were collected.
      if (!more && nested.size === seen) break;
    }
    if (want.size === 0 && !was) return base;
    const nodes = new Map<string, Held & { part: Part | null; issues?: readonly Issue[] }>();
    const upsert: Part = {};
    const evict: string[] = [];
    const changed = new Set<string>();
    const names = itemNamesOf(rc.doc, rc.registry.svs);
    for (const [ref, w] of want) {
      // An item's coordinate system and its `.map` only when something names the item (ST7): else its row's `map`.
      const own = (w as { own?: boolean }).own === true;
      const item = (w.node?.kind as string) === ItemSystemKind ? ref : ref.endsWith('.map') ? ref.slice(0, -4) : ref;
      if (own && !names.refs.has(item) && !names.sources.has(w.source)) continue;
      const old = was?.nodes.get(ref);
      if (old?.node === w.node) {
        // Its issues again: a run reports what it has (a node read once, its faults kept while it is held).
        for (const i of old.issues ?? []) rc.stage.report(i);
        nodes.set(ref, old);
        continue;
      }
      const caught: Issue[] = [];
      const ctx = Object.assign(Object.create(stage) as StageContext, {
        report: (i: Issue) => {
          caught.push(i);
          stage.report(i);
        },
      });
      const part =
        w.node?.kind === 'source' ? null : await loadStreamNode(w.source, w.node, known, ctx, rc.loaders, doc, own);
      nodes.set(ref, { source: w.source, node: w.node, part: part && asPart(part), issues: caught });
      if (part) add(upsert, part);
      if (part) changed.add(w.source);
    }
    // The document's mappings over dynamic nodes, while both their ends are held (upserted again when they return).
    const upserted = new Set<string>();
    // Held: in the document's set, or a dynamic node read into this one (its `/unit` by its molecule's part).
    const live = (k: string) => nodes.get(k)?.part != null;
    const held = (k: string) => base.set.byName(k, t) || live(k) || (k.endsWith('/unit') && live(k.slice(0, -5)));
    for (const [ref, v] of base.deferred ?? []) {
      const part = asPart(v);
      const ends = (part.mappings ?? []).every((m) => held(m.from) && held(m.to));
      const targets = (part.annotations ?? []).every((a) => typeof a.target !== 'string' || held(a.target));
      if (!ends || !targets) continue;
      const old = was?.nodes.get(ref);
      nodes.set(ref, { source: '', node: DEFERRED, part });
      if (old?.part !== part) {
        add(upsert, part);
        upserted.add(ref);
      }
    }
    for (const [ref, old] of was?.nodes ?? [])
      if (!nodes.has(ref) && old.part) {
        keysOfPart(old.part, evict);
        changed.add(old.source);
      }
    // Versions per source, and per source over its subtree, for the templates that read them.
    let [own, tree] = [was?.own ?? NO_VERSIONS, was?.tree ?? NO_VERSIONS];
    if (changed.size > 0) {
      const [o, r] = [new Map(own), new Map(tree)];
      for (const s of changed) {
        o.set(s, (o.get(s) ?? 0) + 1);
        for (let a: string | undefined = s, k = 0; a && k < 64; a = (nested.get(a) ?? was?.nested.get(a))?.source, k++)
          r.set(a, (r.get(a) ?? 0) + 1);
      }
      [own, tree] = [o, r];
    }
    // The document's held-back nodes' issues, kept while they are in (a later run upserts nothing of them).
    const deferredIssues = new Map<string, readonly Issue[]>();
    const done = (r: { set: SeqSet; issues: readonly Issue[] }): StreamAssembly => {
      for (const i of r.issues) rc.stage.report(i);
      for (const [ref] of base.deferred ?? []) {
        if (!nodes.has(ref)) continue;
        const fresh = r.issues.filter((i) => i.message.includes(`"${ref}"`));
        const kept = fresh.length > 0 || upserted.has(ref) ? fresh : (was?.deferredIssues?.get(ref) ?? []);
        if (kept.length > 0) deferredIssues.set(ref, kept);
        if (fresh.length === 0) for (const i of kept) rc.stage.report(i);
      }
      return { set: r.set, parts: base.parts, base, nodes, nested, own, tree, deferredIssues };
    };
    if (!was || was.set.lineage !== base.set.lineage) {
      const all: Part = {};
      for (const v of nodes.values()) if (v.part) add(all, v.part);
      return done(assembleSet(base.set, all, []));
    }
    // A key upserted in this run is neither removed nor evicted (a mapping retargeted between held back and not).
    const removed = was.base === base ? [] : diffParts(was.base.parts, base.parts, was.set, upsert);
    const keep = new Set<string>();
    add({}, upsert, keep);
    const remove = removed.filter((k) => !keep.has(k));
    const gone = evict.filter((k) => !keep.has(k));
    if (empty(upsert) && remove.length === 0 && gone.length === 0) {
      for (const [ref, issues] of was.deferredIssues ?? [])
        if (nodes.has(ref)) for (const i of issues) rc.stage.report(i);
      return was.base === base && nodes.size === was.nodes.size && nested.size === was.nested.size
        ? was
        : { ...was, base, nodes, nested, own, tree };
    }
    return done(assembleSet(was.set, upsert, remove, gone));
  },
};

const entry: NodeKind = { name: '#entry', stage: 'load', inputs() {}, run: (rc) => rc.payload };
const wholeSet: NodeKind<Assembly> = {
  name: '#set',
  stage: 'assemble',
  inputs() {},
  run: (rc) => ({ set: rc.payload as SeqSet, parts: new Map() }),
};

// --- per-view stages --------------------------------------------------------------------------------------------

const sameBytes = (a: Uint8Array | null, b: Uint8Array | null) =>
  a === b || (!!a && !!b && a.length === b.length && a.every((x, k) => x === b[k]));
/** A molecule rebuilt with the same residues and atoms (a numbering edit) is the same content; an alignment by identity. */
const sameContent = (a: unknown, b: unknown) => {
  const [m, n] = [a as Molecule | null, b as Molecule | null];
  return a === b || (!!m?.code && !!n?.code && sameBytes(m.code, n.code) && sameBytes(m.element, n.element));
};
/**
 * A section's result: the resolved section; its axis, the previous one while the `axisKey` is equal; and its content,
 * a token kept while every block's molecule has the same residues and its alignment is the same object (`objects`), so
 * a same-length residue edit reruns what reads residues and nothing else.
 */
export interface SectionResult {
  readonly section: ResolvedSection;
  readonly axis: ResolvedSection;
  readonly lineage: number;
  readonly content: object;
  readonly objects: readonly unknown[];
}
const setOf = (rc: RunContext) => rc.stage.set as SeqSet;
const sectionOf = (doc: DocIndex, n: IndexedNode) => {
  let p = doc.nodes.get(n.parent ?? '');
  while (p && p.kind !== 'section') p = doc.nodes.get(p.parent ?? '');
  return p?.ref ?? '';
};

const section: NodeKind<SectionResult> = {
  name: 'section',
  stage: 'section',
  inputs: (_ref, _n, _doc, out) => out.after('#assemble.stream'),
  // Its view's params (layout), and what it reads of its subtree: blocks, focus, groups, tracks, and of layers only
  // the type, the data (a track's label) and a `rows` layer's params (lanes), so a set on a layer leaves it alone. Of
  // the set (ST1b-b): its data nodes' tables, and the annotations it names (labels, `positions`, a focus,
  // `highlights`), so a dynamic node arriving reruns no section that does not read it.
  key(rc, out) {
    const n = rc.node as IndexedNode;
    const set = setOf(rc);
    out.push(rc.token(n.parent ?? ''), set?.lineage, set?.molecules, set?.units, set?.complexes);
    out.push(set?.alignments, set?.mappings);
    for (const r of [n.params.positions ?? []].flat() as (string | { data: string })[])
      out.push(entity(set, typeof r === 'string' ? r : r.data));
    // Its `highlights`' data (LS1): a bookmark added gives the view a set that has it.
    for (const h of (n.params.highlights as { data: string }[] | undefined) ?? []) out.push(entity(set, h.data));
    const walk = (p: IndexedNode) => {
      for (const c of childrenOf(rc.doc, p)) {
        if (c.kind === 'annotation') continue;
        if (c.kind !== 'layer') {
          out.push(c.ref, rc.token(c.ref));
          if (c.kind === 'focus') out.push(entity(set, String(c.params.annotation)));
          walk(c);
        } else {
          out.push(c.ref, c.params.type, String(c.params.data), c.params.type === 'rows' ? rc.token(c.ref) : 0);
          for (const r of [c.params.data ?? []].flat()) out.push(entity(set, String(r)));
          // A `rows` layer's order and filter read key tables: a changed table reorders the lanes (KS1).
          if (c.params.type === 'rows') for (const r of keyTables(c.params)) out.push(entity(set, r));
        }
      }
    };
    walk(n);
  },
  run(rc) {
    const set = setOf(rc);
    const r = resolveSection(rc.doc, set, rc.ref, rc.stage.registries);
    for (const i of r.issues) rc.stage.report(i);
    const prev = rc.prev as SectionResult | undefined;
    const keep = prev && prev.lineage === set.lineage && prev.axis.axisKey === r.section.axisKey;
    const objects = r.section.blocks.flatMap((b) => [
      b.unit >= 0 ? set.molecules[set.units.molecule[b.unit]] : null,
      b.alignment >= 0 ? set.alignments[b.alignment] : null,
    ]);
    const same =
      keep && objects.length === prev.objects.length && objects.every((o, i) => sameContent(o, prev.objects[i]));
    const content = same ? prev.content : {};
    return { section: r.section, axis: keep ? prev.axis : r.section, lineage: set.lineage, content, objects };
  },
};

const dataOf = (n: IndexedNode, doc: DocIndex) => {
  const own = childrenOf(doc, n).filter((c) => c.kind === 'annotation');
  const refs = n.params.data != null ? [n.params.data].flat().map((r) => String(r)) : own.map((c) => c.ref);
  return refs;
};
const resultOf = (rc: RunContext, n: IndexedNode) => rc.result(sectionOf(rc.doc, n)) as SectionResult | undefined;
const axisOf = (rc: RunContext, n: IndexedNode) => resultOf(rc, n)?.axis ?? null;
/** A theme's name and those nested in a `rules` theme's rules and default. */
function themeNames(r: ThemeRef | null | undefined, out: string[]) {
  if (!r || typeof r.theme !== 'string') return;
  out.push(r.theme);
  if (r.theme !== 'rules') return;
  const q = r.params as { rules?: { theme?: ThemeRef | null }[]; default?: ThemeRef | null } | undefined;
  for (const x of q?.rules ?? []) themeNames(x.theme, out);
  themeNames(q?.default, out);
}

const layer: NodeKind<ResolvedLayer | null> = {
  name: 'layer',
  stage: 'layer',
  inputs(_ref, n, doc, out) {
    if (!n) return;
    out.after(sectionOf(doc, n));
    out.after('#assemble.stream');
    // A `domain_group`'s members' data too: the layer's themes read their union (seqspec:format#themes).
    for (const m of [n, ...groupMembers(doc, n).filter((x) => x !== n)])
      for (const r of dataOf(m, doc)) if (doc.nodes.has(head(r))) out.status(head(r));
  },
  // The axis (not the section's result), what its data and its `of` name in the set, the mappings, alignments and units
  // its paths and bindings go through, the first block's alphabet and laziness, and the registry entries along its
  // type's fallback chain and its channels' themes.
  key(rc, out) {
    const n = rc.node as IndexedNode;
    const set = setOf(rc);
    const axis = axisOf(rc, n);
    out.push(axis, set?.lineage, set?.mappings, set?.alignments, set?.units.count);
    for (const r of [...dataOf(n, rc.doc), ...[n.params.of ?? []].flat()]) out.push(entity(set, String(r)));
    // A `rows` layer's data columns and key tables (computed ones among them) too: a recomputed table redraws it.
    for (const c of (n.params.data_columns ?? []) as { data: string }[]) out.push(entity(set, c.data));
    for (const r of keyTables(n.params)) out.push(entity(set, r));
    for (const m of groupMembers(rc.doc, n)) {
      if (m !== n) for (const r of dataOf(m, rc.doc)) out.push(entity(set, String(r)));
      out.push(m.ref, m.params.domain);
    }
    const b = axis?.blocks[0];
    const aln = b && set.alignments[b.alignment];
    out.push(b && (aln?.alphabet ?? set.molecules[set.units.molecule[b.unit]]?.alphabet), aln?.lazy);
    const svs = rc.stage.registries;
    // `letter_map` (an `emphasis`'s: its track's letters layer's), its entry, and the blocks' alphabets it is checked against.
    const map = letterMapParam(rc.doc, n);
    if (map != null) out.push(map, typeof map === 'string' ? svs.get('letterMaps', map) : null);
    for (const x of map != null ? (axis?.blocks ?? []) : [])
      out.push(set.alignments[x.alignment]?.alphabet ?? set.molecules[set.units.molecule[x.unit]]?.alphabet);
    const names: string[] = [];
    let type: unknown = n.params.type;
    for (let k = 0; typeof type === 'string' && k < 9; k++) {
      const e = svs.get('layerTypes', type);
      out.push(e, rc.registry.get('renderers', type));
      for (const c of e?.channels ?? []) themeNames((n.params[c.name] ?? c.default) as ThemeRef | null, names);
      type = k === 0 && typeof n.params.fallback === 'string' ? n.params.fallback : e?.fallback;
    }
    for (const name of names) out.push(svs.get('themes', name), rc.registry.get('themes', name));
  },
  run(rc) {
    const n = rc.node as IndexedNode;
    const axis = axisOf(rc, n);
    const failed = dataOf(n, rc.doc).find((r) => rc.failed(head(r)));
    if (!axis || failed) {
      const message = `"${n.ref}" is omitted: its data "${failed}" is unavailable`;
      if (failed) rc.stage.report({ code: 'svs.render.omitted', severity: 'warning', path: n.path, message });
      return null;
    }
    const r = resolveLayer(axis, n.ref, rc.doc, setOf(rc), rc.stage.registries, rc.supported(axis.layout));
    for (const i of r.issues) rc.stage.report(i);
    return r.layer;
  },
};

/** One layer's input to its renderer and themes: the layer against the axis it was resolved on. */
function layerInput(rc: RunContext, l = rc.result((rc.node as IndexedNode).ref) as ResolvedLayer | null | undefined) {
  const axis = axisOf(rc, rc.node as IndexedNode);
  const set = setOf(rc);
  if (!l || !axis || l.type < 0) return null;
  const annotations = l.data.map((id) => set.annotations[id] as Annotation);
  return { set, section: axis, layer: l, annotations, params: l.params } as LayerInput;
}
/** A layer's prepared state by its renderer (null: none drawn). */
function prepared(rc: RunContext, input: LayerInput | null): Prepared | null {
  const renderer = input && rc.registry.get('renderers', input.set.strings.get(input.layer.type));
  if (!input || !renderer) return null;
  // What a prepare leaves out is reported on the layer's node.
  const path = (rc.node as IndexedNode).path;
  const report = (i: Issue) => rc.stage.report({ ...i, path: i.path || path });
  return { renderer, state: renderer.prepare ? renderer.prepare(input, { registry: rc.registry, report }) : null };
}
/** A theme per channel of the layer's type, null where its reference is (themes.ts). */
function themed(rc: RunContext, input: LayerInput | null, report: (i: Issue) => void): Themes {
  if (!input) return [];
  const entry = rc.stage.registries.get('layerTypes', input.set.strings.get(input.layer.type));
  const tc = themeContext(input, rc.skin.dark, rc.registry, report);
  // A color scheme picked at runtime takes the color channel's place (seqspec:interaction#color-schemes-at-runtime).
  const scheme = rc.colorScheme(input.layer.ref);
  const range = rc.colorRange(input.layer.ref);
  const drawn = rc.registry.get('renderers', input.set.strings.get(input.layer.type))?.sized;
  const themes = input.layer.themes.map((given, i) => {
    const channel = entry?.channels[i];
    const color = channel?.kind === 'color';
    const picked = scheme !== undefined && color ? pickedTheme(given, scheme) : given;
    const ref = range && color && picked ? withRange(picked, range) : picked;
    // A size channel its renderer does not draw: said once, drawn at the viewer's size.
    if (ref && channel?.kind === 'size' && !drawn?.includes(channel.name)) {
      report(
        Issue.warning('svs.render.theme', input.layer.ref, `"${channel.name}" is drawn at the viewer's size here`),
      );
      return null;
    }
    return ref && channel ? createTheme(tc, ref, channel, rc.registry, report) : null;
  });
  subjects.set(themes, subjectOf(tc));
  return themes;
}
/** The section's content token, unless every entry the node uses declares `content: false`. */
const contentOf = (rc: RunContext, entries: readonly ({ readonly content?: boolean } | undefined)[]) =>
  entries.some((e) => e && e.content !== false) ? resultOf(rc, rc.node as IndexedNode)?.content : null;
const layerOf = (rc: RunContext) => rc.result((rc.node as IndexedNode).ref) as ResolvedLayer | null;
interface Prepared {
  readonly renderer: LayerRenderer;
  readonly state: unknown;
}

// The prepare and theme nodes read their section's result through the layer, which waits for it.
const prepare: NodeKind<Prepared | null> = {
  name: '#prepare',
  stage: 'prepare',
  inputs: (_ref, n, _doc, out) => n && out.result(n.ref),
  key(rc, out) {
    const l = layerOf(rc);
    const renderer = l && l.type >= 0 ? rc.registry.get('renderers', setOf(rc).strings.get(l.type)) : null;
    out.push(renderer, contentOf(rc, [renderer ?? undefined]));
  },
  run: (rc) => prepared(rc, layerInput(rc)),
  dispose: (r) => r?.renderer.dispose?.(r.state),
};

type Themes = readonly (Theme | null)[];
/**
 * What each theme list's layer color is about (the scheme picker's), kept beside the list: a property on the array
 * would give the draws' theme lists another shape.
 */
const subjects = new WeakMap<Themes, SchemeSubject | null>();
const providersOf = (rc: RunContext, l: ResolvedLayer | null) => {
  const names: string[] = [];
  for (const r of l?.themes ?? []) themeNames(r, names);
  return names.map((name) => rc.registry.get('themes', name));
};
/** A theme per channel of the layer's type (null where its reference is), down its fallbacks (themes.ts). */
const theme: NodeKind<Themes> = {
  name: '#theme',
  stage: 'theme',
  inputs: (_ref, n, _doc, out) => n && out.result(n.ref),
  key(rc, out) {
    const l = layerOf(rc);
    const providers = providersOf(rc, l);
    const scheme = rc.colorScheme(rc.node?.ref ?? '');
    // A scheme picked at runtime is drawn by the `scheme` theme: its content counts as the document's themes' do.
    if (scheme !== undefined) providers.push(rc.registry.get('themes', 'scheme'));
    // What its themes read in the set besides the layer's data: a `rules` column or profile (a computed node's answer).
    const set = setOf(rc);
    for (const r of refsOf(rc.node as IndexedNode, rc.doc, rc.registry.svs)) out.push(entity(set, r));
    out.push(
      rc.colorRange(rc.node?.ref ?? ''),
      rc.skin,
      contentOf(rc, providers),
      ...providers,
      ...schemeEntries(rc.registry.svs, l?.themes ?? [], scheme),
    );
  },
  run: (rc) => themed(rc, layerInput(rc), (i) => rc.stage.report({ ...i, path: (rc.node as IndexedNode).path })),
  equals: (a, b) => a.length === b.length && a.every((t, i) => t === b[i]),
};

const layersOf = (doc: DocIndex, s: IndexedNode): IndexedNode[] =>
  childrenOf(doc, s).flatMap((c) => (c.kind === 'layer' ? [c] : c.kind === 'annotation' ? [] : layersOf(doc, c)));

/** A template instance (seqspec:format#placement): its layers resolved with its node, prepared and themed; shared by
 * successive results while its node and what it reads are unchanged, disposed when the last lets go. */
interface Instance {
  readonly node: string;
  readonly entity: unknown;
  readonly layers: readonly ResolvedLayer[];
  readonly prepared: readonly (Prepared | null)[];
  readonly themes: readonly Themes[];
  /** Shown: it places something, or the template keeps every node (`reaching: false`). */
  readonly kept: boolean;
  /** Per layer, the color scheme picked for it at runtime (by its instance ref, as the legend's picker sets it). */
  readonly schemes: readonly (string | undefined)[];
  /** Likewise its color range set at runtime. */
  readonly ranges: readonly (readonly [number, number] | undefined)[];
  holders: number;
}
interface TemplateResult {
  readonly env: readonly unknown[];
  readonly all: ReadonlyMap<string, Instance>;
  /** The instances shown, in order. */
  readonly shown: readonly Instance[];
}
const streamOf = (rc: RunContext) => rc.result('#assemble.stream') as Partial<StreamAssembly> | undefined;
const sameList = (a: readonly unknown[], b: readonly unknown[]) =>
  a.length === b.length && a.every((x, i) => x === b[i]);
/** What every instance reads besides its node: the axis and content, the set's tables paths go through, registries,
 * the skin, the template's params. */
function templateEnv(rc: RunContext, n: IndexedNode): unknown[] {
  const set = setOf(rc);
  const axis = axisOf(rc, n);
  const env: unknown[] = [axis, resultOf(rc, n)?.content, set?.lineage, set?.mappings, set?.alignments, set?.units];
  env.push(rc.skin, rc.registry.version.value, rc.stage.registries.version, axis && rc.supported(axis.layout));
  for (const c of [n, ...childrenOf(rc.doc, n)]) env.push(rc.token(c.ref));
  return env;
}

/**
 * `#template:<track>` (ST1b-b): a track template's instances, one per node `templateNodes` selects, each resolved,
 * prepared and themed once and kept while its node is the same object; an instance placing nothing is left out under
 * `reaching`. It reruns when its source's nodes change (`own`, or `tree` over nested sources), never
 * for another source's, so an arrival reruns the views its templates cover and nothing else.
 */
const template: NodeKind<TemplateResult> = {
  name: '#template',
  stage: 'layer',
  inputs(_ref, n, doc, out) {
    if (!n) return;
    out.after(sectionOf(doc, n));
    out.after('#assemble.stream');
  },
  key(rc, out) {
    const n = rc.node as IndexedNode;
    const fe = forEachOf(n);
    const sa = streamOf(rc);
    for (const x of templateEnv(rc, n)) out.push(x);
    out.push(fe && (fe.depth === 0 ? sa?.own : sa?.tree)?.get(fe.source));
    // The key tables its order and filter read (KS1): a changed table reorders the instances.
    if (fe) for (const r of keyTables(fe as never)) out.push(entity(setOf(rc), r));
    // A scheme picked on an instance reruns its template; only that instance is themed again.
    out.push(rc.colorSchemes, rc.colorRanges);
  },
  run(rc) {
    const n = rc.node as IndexedNode;
    const fe = forEachOf(n);
    const axis = axisOf(rc, n);
    const set = setOf(rc);
    const env = templateEnv(rc, n);
    const prev = rc.prev as TemplateResult | undefined;
    const reuse = prev && sameList(prev.env, env) ? prev.all : null;
    const all = new Map<string, Instance>();
    const shown: Instance[] = [];
    if (!fe || !axis) return { env, all, shown };
    const ls = childrenOf(rc.doc, n).filter((c) => c.kind === 'layer');
    // An instance's issues are its template layer's, reported once a run (when first read).
    const seen = new Set<string>();
    const report = (i: Issue) => {
      const k = `${i.code} ${i.path}`;
      if (!seen.has(k)) rc.stage.report(i);
      seen.add(k);
    };
    const supported = rc.supported(axis.layout);
    const say = (where: string, message: string) =>
      report({ code: 'svs.layer.data', severity: 'error', path: `${n.path}/params/for_each/${where}`, message });
    for (const node of templateNodes(n.ref, rc.doc, set, streamOf(rc)?.nested, rc.stage.registries, say)) {
      const e = entity(set, node);
      let x = reuse?.get(node);
      const picked = (r: string) => rc.colorScheme(r);
      if (
        !x ||
        x.entity !== e ||
        x.layers.some((l, k) => picked(l.ref) !== x?.schemes[k] || rc.colorRange(l.ref) !== x?.ranges[k])
      ) {
        const layers = ls.map((l) => {
          const r = resolveLayer(axis, l.ref, rc.doc, set, rc.stage.registries, supported, node);
          for (const i of r.issues) report(i);
          return r.layer;
        });
        const inputs = layers.map((l) => layerInput(rc, l));
        const themes = inputs.map((i, k) => themed(rc, i, (x: Issue) => report({ ...x, path: ls[k].path })));
        const kept = !fe.reaching || layers.some((l) => l.placements.length > 0);
        const schemes = layers.map((l) => picked(l.ref));
        const ranges = layers.map((l) => rc.colorRange(l.ref));
        x = {
          node,
          entity: e,
          layers,
          prepared: inputs.map((i) => prepared(rc, i)),
          themes,
          kept,
          schemes,
          ranges,
          holders: 0,
        };
      }
      x.holders++;
      all.set(node, x);
      if (x.kept) shown.push(x);
    }
    return { env, all, shown };
  },
  equals: (a, b) => a.all.size === b.all.size && sameList(a.shown, b.shown),
  dispose(r) {
    for (const x of r.all.values()) if (--x.holders === 0) for (const p of x.prepared) p?.renderer.dispose?.(p.state);
  },
};

const view: NodeKind<ViewModel> = {
  name: '#view',
  stage: 'view',
  inputs(_ref, n, doc, out) {
    if (!n) return;
    out.after('#assemble.stream');
    for (const s of childrenOf(doc, n).filter((c) => c.kind === 'section')) {
      out.result(s.ref);
      const templates = new Set<string>();
      for (const l of layersOf(doc, s)) {
        const track = doc.nodes.get(l.parent ?? '');
        if (track && forEachOf(track)) templates.add(`#template:${track.ref}`);
        else for (const r of [l.ref, `#prepare:${l.ref}`, `#theme:${l.ref}`]) out.result(r);
      }
      for (const r of templates) out.result(r);
    }
  },
  run(rc) {
    const n = rc.node as IndexedNode;
    const set = setOf(rc);
    const sections: ResolvedSection[] = [];
    const layers: ResolvedLayer[][] = [];
    const prepared: PreparedLayer[] = [];
    for (const c of childrenOf(rc.doc, n)) {
      const r = c.kind === 'section' ? (rc.result(c.ref) as SectionResult | undefined) : undefined;
      if (!r) continue;
      // Each template's instances at its place.
      const instances = new Map<number, string[]>();
      const own = new Map<string, { l: ResolvedLayer; p: Prepared | null; themes: Themes }>();
      r.section.tracks.forEach((tr, k) => {
        const res = rc.result(`#template:${set.strings.get(tr.ref)}`) as TemplateResult | undefined;
        if (!res) return;
        instances.set(
          k,
          res.shown.map((x) => x.node),
        );
        for (const x of res.shown)
          for (let i = 0; i < x.layers.length; i++)
            own.set(x.layers[i].ref, { l: x.layers[i], p: x.prepared[i], themes: x.themes[i] });
      });
      const s = instantiate(r.section, instances, rc.doc, set);
      const ls: ResolvedLayer[] = [];
      for (const { ref } of s.layers) {
        const x = own.get(ref);
        const l = x ? x.l : (rc.result(ref) as ResolvedLayer | null | undefined);
        if (!l) continue;
        ls.push(l);
        const p = x ? x.p : (rc.result(`#prepare:${ref}`) as Prepared | null | undefined);
        const themes = x ? x.themes : ((rc.result(`#theme:${ref}`) as Themes | undefined) ?? []);
        const subject = subjects.get(themes) ?? null;
        prepared.push({
          renderer: l.type >= 0 ? set.strings.get(l.type) : '',
          state: p?.state ?? null,
          themes,
          subject,
        });
      }
      sections.push(s);
      layers.push(ls);
    }
    const issues: Issue[] = [];
    const resolved = composeView(n.ref, sections, layers, set, issues);
    for (const i of issues) rc.stage.report(i);
    return { ref: n.ref, resolved, prepared, set };
  },
};

/** Every built-in node kind: one per document kind that loads or resolves, and per derived stage. */
export const DefaultNodeKinds: readonly NodeKind[] = [
  ...['source', 'sequence', 'ligand', 'unit', 'alignment', 'mapping', 'annotation'].map((k) => load(k)),
  load('complex', 'complex'),
  computed,
  assembly('#assemble.base', null, [['complex'], ['sequence', 'ligand', 'unit', 'alignment']]),
  assembly('#assemble.mapped', '#assemble.base', [['alignment', 'annotation', 'mapping']], true),
  assembly('#assemble', '#assemble.mapped', [['computed']]),
  streamAssembly,
  entry,
  wholeSet,
  section,
  layer,
  prepare,
  theme,
  template,
  view,
] as NodeKind[];
