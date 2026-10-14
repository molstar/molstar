/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static phase of validation (seqspec:format#validation): every issue in one pass, stable codes, JSON pointers,
// messages that say what is wrong and what would be valid. Nothing is repaired.

import { Issue, KeyKind, type KeyResolver, type KeyTarget, MoleculeKind } from "@molstar/seq-core/model";
import { type DocNode, type SeqViewSpecDoc, SvsFormat } from "./doc.js";
import { snapshotIssues } from "./patch/patch.js";
import { registeredNamespaces, type SvsRegistries } from "./registry.js";
import {
  float,
  isNamespace,
  isPlainObject,
  isRefSyntax,
  json,
  list,
  literal,
  nullable,
  ptr,
  type RefUse,
  report,
  reportUnknown,
  scratch,
  str,
  TEMPLATE_NODE,
  typeError,
} from "./schema/field-schema.js";
import { checkParams, OptionalField, object, RequiredField, SimpleParamsSchema } from "./schema/params-schema.js";
import { type IndexedNode, nodeList, type TreeContext, type TreeSchema } from "./schema/tree-schema.js";
import { SvsTreeSchema } from "./tables/svs-tree.js";

const optStr = (description: string) => OptionalField(nullable(str), null, description);
const descriptionFormat = OptionalField(literal("markdown", "plaintext"), "markdown", "");

/** The envelope (seqspec:format#the-document); `format`, `version` and `root` are checked apart. */
export const EnvelopeSchema = SimpleParamsSchema({
  format: RequiredField(json, '"seq_view_spec"'),
  version: RequiredField(json, "the format version the document was written for (semver)"),
  revision: optStr("the owner's name for this version of the document; a patch's base"),
  metadata: OptionalField(
    object({
      title: optStr(""),
      description: optStr("markdown, sanitized"),
      description_format: descriptionFormat,
      timestamp: optStr(""),
      authors: OptionalField(nullable(list(str)), null, ""),
    }),
    {},
    "never read for meaning",
  ),
  requires: OptionalField(list(str), [], "extension namespaces needed to interpret the document"),
  root: RequiredField(json, "the document's tree"),
  snapshots: OptionalField(
    list(
      object({
        ref: RequiredField(str, "unique among the snapshots"),
        title: optStr(""),
        description: optStr(""),
        description_format: descriptionFormat,
        ops: OptionalField(list(object()), [], "a patch over the base tree, without base"),
        duration_ms: OptionalField(float, 1000, "how long a player stays on it"),
        transition_ms: OptionalField(float, 0, "the transition into it"),
      }),
    ),
    [],
    "an ordered story over the tree",
  ),
});

/** How deep a document's JSON may nest (seqspec:format#the-document): a deeper one is refused before it is read. */
export const MAX_JSON_DEPTH = 256;

/** Each object's or array's height, once: documents are immutable once passed, so an edit rereads only its spine. */
const HEIGHT = new WeakMap<object, number>();

interface Frame {
  readonly o: Record<string, unknown> | unknown[];
  readonly keys: readonly string[] | null;
  i: number;
  h: number;
  readonly path: string;
}

const isNode = (v: unknown): v is Record<string, unknown> | unknown[] =>
  v !== null && typeof v === "object" && !ArrayBuffer.isView(v);

/**
 * The path of a value nested deeper than `max` objects and arrays (typed arrays are leaves), or null; found without
 * recursing, so a hostile document of any depth (JSON.parse takes 20,000) yields an issue, not a stack overflow (CM-03).
 */
export function tooDeep(value: unknown, max = MAX_JSON_DEPTH): string | null {
  if (!isNode(value)) return null;
  const stack: Frame[] = [];
  const push = (o: Record<string, unknown> | unknown[], path: string): boolean => {
    const known = HEIGHT.get(o);
    const parent = stack[stack.length - 1];
    if (known !== undefined && stack.length + known <= max) {
      if (parent) parent.h = Math.max(parent.h, known);
      return true;
    }
    if (stack.length >= max) return false;
    stack.push({ o, keys: Array.isArray(o) ? null : Object.keys(o), i: 0, h: 0, path });
    return true;
  };
  if (!push(value, "")) return "";
  while (stack.length > 0) {
    const f = stack[stack.length - 1];
    const n = f.keys ? f.keys.length : (f.o as unknown[]).length;
    if (f.i < n) {
      const at = f.i++;
      const v = f.keys ? (f.o as Record<string, unknown>)[f.keys[at]] : (f.o as unknown[])[at];
      if (!isNode(v)) continue;
      // The key's text only for a nested value (LA-18: a string per element of a long array of numbers).
      const path = `${f.path}/${f.keys ? f.keys[at] : at}`;
      if (!push(v, path)) return path;
      continue;
    }
    stack.pop();
    HEIGHT.set(f.o, f.h + 1);
    const parent = stack[stack.length - 1];
    if (parent) parent.h = Math.max(parent.h, f.h + 1);
  }
  return null;
}

/**
 * `format` and `version` (seqspec:format#versioning), and the depth bound. Null: the document cannot be read further;
 * else whether it is read leniently.
 */
export function checkEnvelope(doc: unknown, issues: Issue[], strict = false): boolean | null {
  if (!isPlainObject(doc) || doc.format !== SvsFormat) {
    issues.push(Issue.error("svs.format.invalid", "/format", `not a SeqViewSpec: format must be "${SvsFormat}"`));
    return null;
  }
  const deep = tooDeep(doc);
  if (deep !== null) {
    const message = `the document nests deeper than ${MAX_JSON_DEPTH} objects and lists here: it cannot be read`;
    issues.push(Issue.error("svs.format.invalid", deep, message));
    return null;
  }
  const m = typeof doc.version === "string" ? /^(\d+)\.(\d+)\.(\d+)$/.exec(doc.version) : null;
  if (!m) {
    issues.push(Issue.error("svs.version.missing", "/version", 'version must be the format\'s semver, "1.0.0"'));
    return null;
  }
  if (m[1] !== "1") {
    issues.push(
      Issue.error("svs.version.unsupported_major", "/version", `version ${doc.version}: this reader reads 1.x`),
    );
    return null;
  }
  if (m[2] === "0") return false;
  const message = `written for ${doc.version}, read as 1.0.0: unknown params, kinds, types and themes are warnings`;
  issues.push(Issue.warning("svs.version.newer_minor", "/version", message));
  return !strict;
}

export interface ValidateOptions {
  /** A newer minor read strictly: unknowns are errors. */
  strict?: boolean;
  /** Refs may name these keys (a set's, or a resolver's over build inputs): data without a document. */
  set?: KeyResolver;
  /**
   * The document is read from an archive, its main file in this directory of it (`""`: the root): a source's `uri` may
   * be a path inside it (seqspec:sources#archives). `readSvsx` gives it.
   */
  base?: string;
}

/** Every static issue of a document. */
export function validate(doc: unknown, reg: SvsRegistries, opts: ValidateOptions = {}): Issue[] {
  const issues: Issue[] = [];
  const lenient = checkEnvelope(doc, issues, opts.strict);
  if (lenient === null) return issues;
  const d = doc as Record<string, unknown>;
  const nodes = new Map<string, IndexedNode>();
  const ctx: TreeContext = { issues, registries: reg, refs: [], lenient, from: "", depth: 0, nodes, base: opts.base };
  checkParams(EnvelopeSchema, d, "", ctx);
  const schema: TreeSchema = SvsTreeSchema;
  const order: IndexedNode[] = [];
  if (d.root !== undefined) indexNode(d.root, null, "/root", schema, ctx, order);
  for (const n of order) checkNode(n, schema, ctx);
  resolveRefs(ctx, opts.set);
  const known = registeredNamespaces(reg);
  if (Array.isArray(d.requires)) {
    d.requires.forEach((ns, i) => {
      if (typeof ns === "string" && !known.has(ns)) {
        report(ctx, "svs.requires.unknown", `/requires/${i}`, `namespace "${ns}" is required but not registered`);
      }
    });
  }
  if (Array.isArray(d.snapshots)) checkSnapshots(d.snapshots, ctx);
  if (Array.isArray(d.snapshots) && !issues.some((i) => i.severity === "error")) {
    for (const i of snapshotIssues(doc as SeqViewSpecDoc, reg, { set: opts.set, base: opts.base })) issues.push(i);
  }
  return issues;
}

const NODE_KEYS = new Set(["kind", "ref", "params", "custom", "children"]);

function indexNode(
  v: unknown,
  parent: IndexedNode | null,
  path: string,
  schema: TreeSchema,
  ctx: TreeContext,
  order: IndexedNode[],
) {
  if (!isPlainObject(v) || typeof v.kind !== "string")
    return typeError(ctx, path, "a node {kind, ref, params, ...}", v);
  const kind = v.kind;
  const ns = Object.hasOwn(schema.nodes, kind) ? schema.nodes[kind] : undefined;
  if (!ns) {
    const kinds = Object.keys(schema.nodes).join(", ");
    return reportUnknown(ctx, "svs.node.kind", ptr(path, "kind"), `unknown kind "${kind}"; kinds: ${kinds}`);
  }
  if (parent === null ? kind !== schema.rootKind : !ns.parent.includes(parent.kind)) {
    const allowed = ns.parent.length > 0 ? ns.parent.join(", ") : "none: it is the root";
    report(
      ctx,
      "svs.node.parent",
      path,
      `${kind} cannot sit ${parent ? `under ${parent.kind}` : "at the root"}; its parents: ${allowed}`,
    );
  }
  let ref = v.ref;
  if (ref === undefined && ns.ref === "implied") ref = kind;
  else if (ref === undefined) report(ctx, "svs.ref.missing", path, `this ${kind} needs a ref`);
  else if (typeof ref !== "string" || !isRefSyntax(ref, false) || (ns.ref === "implied" && ref !== kind)) {
    const why = ns.ref === "implied" ? `the ${kind}'s ref is "${kind}"` : "a ref matches [A-Za-z0-9][A-Za-z0-9_.-]*";
    report(ctx, "svs.ref.syntax", ptr(path, "ref"), `${JSON.stringify(ref)}: ${why}`);
    ref = undefined;
  } else if (ctx.nodes.has(ref)) {
    report(ctx, "svs.ref.duplicate", ptr(path, "ref"), `"${ref}" is already the ref of ${ctx.nodes.get(ref)?.path}`);
    ref = undefined;
  }
  for (const k of Object.keys(v)) {
    if (!NODE_KEYS.has(k))
      report(ctx, "svs.params.unknown", ptr(path, k), `a node has kind, ref, params, custom, children`);
  }
  const params = isPlainObject(v.params) ? v.params : {};
  if (v.params !== undefined && !isPlainObject(v.params)) typeError(ctx, ptr(path, "params"), "object", v.params);
  const n: IndexedNode = {
    node: v as unknown as DocNode,
    ref: (ref as string) ?? "",
    kind,
    params,
    parent: parent?.ref ?? null,
    path,
  };
  if (n.ref) (ctx.nodes as Map<string, IndexedNode>).set(n.ref, n);
  order.push(n);
  if (v.children !== undefined && !Array.isArray(v.children)) {
    typeError(ctx, ptr(path, "children"), "a list of nodes", v.children);
  }
  const kids = nodeList(v.children);
  kids.forEach((c, i) => {
    indexNode(c, n, ptr(ptr(path, "children"), i), schema, ctx, order);
  });
  for (const [k, child] of Object.entries(schema.nodes)) {
    if (!child.count || !child.parent.includes(kind)) continue;
    const c = kids.filter((x) => isPlainObject(x) && x.kind === k).length;
    const ok = child.count === "one" ? c === 1 : child.count === "some" ? c >= 1 : c <= 1;
    const need = child.count === "one" ? "exactly one" : child.count === "some" ? "one or more" : "at most one";
    if (!ok) report(ctx, "svs.node.count", path, `this ${kind} has ${c} ${k} children; it needs ${need}`);
  }
}

/**
 * A node's own checks (its params, custom and inline tables), memoized per `params` object. The key is what they read
 * besides the node's own objects: its kind, its `custom` object, the registries and their version, and whether the
 * document is read leniently. Documents are immutable once passed, so an unchanged `params` object gives the same
 * issues and refs; they are kept with the node's path and ref and moved to its current ones on reuse (a sibling
 * inserted before it shifts its path).
 */
interface Own {
  readonly kind: string;
  readonly custom: unknown;
  readonly reg: SvsRegistries;
  readonly version: number;
  readonly lenient: boolean;
  readonly base: string | undefined;
  path: string;
  from: string;
  issues: Issue[];
  refs: RefUse[];
}
const OWN = new WeakMap<object, Own[]>();
/** Memo entries kept per `params` object: one per kind, registries, leniency and archive base that read it lately. */
const OWN_KEPT = 4;

function ownChecks(n: IndexedNode, schema: TreeSchema, ctx: TreeContext): Own {
  const params = n.node.params;
  const custom = n.node.custom;
  const reg = ctx.registries;
  const memo = isPlainObject(params) ? (OWN.get(params) ?? []) : null;
  const hit = memo?.find(
    (o) =>
      o.kind === n.kind &&
      o.custom === custom &&
      o.reg === reg &&
      o.version === reg.version &&
      o.lenient === ctx.lenient &&
      o.base === ctx.base,
  );
  if (hit) {
    if (hit.path !== n.path || hit.from !== n.ref) {
      const was = hit.path.length;
      const move = (p: string) => (p.startsWith(hit.path) ? n.path + p.slice(was) : p);
      hit.issues = hit.issues.map((i) => ({ ...i, path: move(i.path) }));
      hit.refs = hit.refs.map((r) => ({ ...r, path: move(r.path), from: n.ref }));
      hit.path = n.path;
      hit.from = n.ref;
    }
    return hit;
  }
  const sub: TreeContext = { ...ctx, issues: [], refs: [], from: n.ref, depth: 0 };
  checkOwn(n, schema, sub);
  const own: Own = {
    kind: n.kind,
    custom,
    reg,
    version: reg.version,
    lenient: ctx.lenient,
    base: ctx.base,
    path: n.path,
    from: n.ref,
    issues: sub.issues,
    refs: sub.refs,
  };
  if (memo && params) {
    if (memo.length >= OWN_KEPT) memo.shift();
    memo.push(own);
    OWN.set(params, memo);
  }
  return own;
}

function checkNode(n: IndexedNode, schema: TreeSchema, ctx: TreeContext) {
  const own = ownChecks(n, schema, ctx);
  for (const i of own.issues) ctx.issues.push(i);
  for (const r of own.refs) ctx.refs.push(r);
  ctx.from = n.ref;
  // The kind's own rules read other nodes (`svs.layer.data`, `.split`): run on every call.
  schema.nodes[n.kind].check?.(n, ctx);
}

function checkOwn(n: IndexedNode, schema: TreeSchema, ctx: TreeContext) {
  const ns = schema.nodes[n.kind];
  const at = ptr(n.path, "params");
  checkParams(ns.params, n.params, at, ctx);
  const owned = ns.owned?.schema(n.params, ctx.registries);
  if (ns.owned && owned) checkParams(owned, n.params[ns.owned.param] ?? {}, ptr(at, ns.owned.param), ctx);
  for (const set of ns.exclusive ?? []) {
    const given = set.filter((k) => n.params[k] != null);
    if (given.length !== 1) {
      report(
        ctx,
        "svs.params.exclusive",
        at,
        `give exactly one of ${set.join(", ")} (given: ${given.join(", ") || "none"})`,
      );
    }
  }
  const custom = n.node.custom;
  if (custom !== undefined) checkCustom(n, custom, ptr(n.path, "custom"), ctx);
}

/** `custom` by namespace (seqspec:format#extensions): registered tables checked, other namespaces passed through. */
function checkCustom(n: IndexedNode, custom: unknown, path: string, ctx: TreeContext) {
  if (!isPlainObject(custom)) return typeError(ctx, path, "an object keyed by namespace", custom, "svs.custom.invalid");
  for (const [ns, value] of Object.entries(custom)) {
    const at = ptr(path, ns);
    if (!isNamespace(ns) || ns === "svs" || ns === "seqv") {
      report(
        ctx,
        "svs.custom.namespace",
        at,
        `"${ns}" is not a namespace (a short name or reverse-DNS; not svs, seqv)`,
      );
    } else if (!isPlainObject(value)) {
      typeError(ctx, at, `an object of "${ns}" fields`, value, "svs.custom.invalid");
    } else {
      const entry = ctx.registries.get("custom", `${ns}/${n.kind}`);
      if (!entry) {
        if (ctx.registries.list("custom").some((e) => e.namespace === ns)) {
          report(ctx, "svs.custom.invalid", at, `namespace "${ns}" has no fields on ${n.kind} nodes`);
        }
        continue;
      }
      const sub = scratch(ctx);
      checkParams(entry.schema, value, at, sub);
      for (const r of sub.refs) ctx.refs.push(r);
      for (const i of sub.issues)
        ctx.issues.push(i.code.startsWith("svs.params.") ? { ...i, code: "svs.custom.invalid" } : i);
    }
  }
}

/**
 * A dynamic node's kinds (seqspec:sources#dynamic-nodes): what a streamed source's derived refs may name, known when
 * held.
 */
const DYNAMIC = ["sequence", "ligand", "unit", "mapping", "annotation", "source"];
/** A `source` node that streams (seqspec:sources#streamed-sources): its derived refs are its dynamic nodes. */
export const isStreamedSource = (n: { kind: string; params: Record<string, unknown> } | undefined) =>
  n?.kind === "source" && (n.params.stream != null || n.params.protocol != null);

/** What a derived ref's suffix can name, by the node it is under (seqspec:format#nodes-and-refs). */
function yields(n: IndexedNode, suffix: string): string[] {
  if ((n.kind === "sequence" || n.kind === "ligand") && suffix === "unit") return ["unit"];
  if (isStreamedSource(n)) return DYNAMIC;
  // A provider's complex, or a letters stream's (`source_ref`, its units by the manifest's keys; ST4b).
  if (n.kind === "complex" && (n.params.provider != null || n.params.source_ref != null))
    return ["sequence", "ligand", "unit"];
  return n.kind === "alignment" && /^[0-9]+$/.test(suffix) ? ["row"] : [];
}

/** By KeyKind; a molecule is a sequence unless a ligand. Shared: a table's target column resolves a ref per row. */
const KEY_KINDS: readonly (readonly string[])[] = [
  ["sequence"],
  ["unit"],
  ["alignment"],
  ["row"],
  ["mapping"],
  ["annotation", "computed"],
  ["complex"],
];
const LIGAND: readonly string[] = ["ligand"];

/** The kinds of node a set's key can stand for: a molecule is a sequence or a ligand, an annotation key either kind. */
function kindsOfKey(set: KeyResolver, t: KeyTarget): readonly string[] {
  if (t.kind === KeyKind.molecule && set.molecules[t.id]?.kind === MoleculeKind.ligand) return LIGAND;
  return KEY_KINDS[t.kind] ?? KEY_KINDS[KeyKind.unit];
}
const overlaps = (a: readonly string[], b: readonly string[]) => {
  for (let i = 0; i < a.length; i++) if (b.includes(a[i])) return true;
  return false;
};

function resolveRefs(ctx: TreeContext, set: KeyResolver | undefined) {
  const edges = new Map<string, Set<string>>();
  const target: KeyTarget = { kind: -1, id: -1, row: -1 };
  // The hint for an unresolved ref, once per kinds list (not a scan of every node per ref).
  const hints = new Map<readonly string[], string>();
  const hintOf = (kinds: readonly string[]) => {
    let hint = hints.get(kinds);
    if (hint === undefined) {
      const like = [...ctx.nodes.values()].filter((x) => kinds.includes(x.kind)).map((x) => x.ref);
      hint = like.length > 0 ? `; ${kinds.join(" or ")} refs here: ${like.slice(0, 6).join(", ")}` : "";
      hints.set(kinds, hint);
    }
    return hint;
  };
  for (const u of ctx.refs) {
    if (u.ref === TEMPLATE_NODE) {
      // A template's node: only in a layer of a track with `for_each` (seqspec:format#placement).
      const track = ctx.nodes.get(ctx.nodes.get(u.from)?.parent ?? "");
      if (track?.params.for_each == null)
        report(ctx, "svs.ref.syntax", u.path, `"${u.ref}" stands for a template's node: only in a track with for_each`);
      continue;
    }
    const rows = u.count !== undefined && u.count > 1 ? ` (${u.count} rows)` : "";
    const slash = u.ref.indexOf("/");
    const head = slash < 0 ? u.ref : u.ref.slice(0, slash);
    const n = ctx.nodes.get(head);
    if (!n) {
      if (set?.byName(u.ref, target)) {
        const kinds = kindsOfKey(set, target);
        if (!overlaps(kinds, u.kinds)) {
          const message = `"${u.ref}" is a ${kinds[0]}; expected ${u.kinds.join(" or ")}${rows}`;
          report(ctx, "svs.ref.wrong_kind", u.path, message);
        }
        continue;
      }
      // A streamed source's dynamic node, not held yet: its kind is known when it arrives.
      if (set?.dynamic?.(u.ref) && overlaps(DYNAMIC, u.kinds)) continue;
      report(ctx, "svs.ref.unresolved", u.path, `"${u.ref}" names no node${hintOf(u.kinds)}${rows}`);
      continue;
    }
    const kinds = slash < 0 ? [n.kind] : yields(n, u.ref.slice(slash + 1));
    if (!kinds.some((k) => u.kinds.includes(k))) {
      const is =
        slash < 0
          ? `of kind ${n.kind}`
          : kinds.length > 0
            ? `a derived ${kinds.join(" or ")}`
            : `nothing ${n.kind} nodes yield`;
      report(ctx, "svs.ref.wrong_kind", u.path, `"${u.ref}" is ${is}; expected ${u.kinds.join(" or ")}${rows}`);
      continue;
    }
    // A section's `follow` links windows, not data: its loops are resolve's (`svs.section.follow`), never a cycle here.
    if (u.kinds.length === 1 && u.kinds[0] === "section") continue;
    // A set per node: a table's target column is a ref per row (a genome's GC: 3M), so no copy per ref (was quadratic).
    const to = edges.get(u.from);
    if (to) to.add(head);
    else edges.set(u.from, new Set([head]));
  }
  const state = new Map<string, 1 | 2>();
  const visit = (r: string, trail: string[]) => {
    state.set(r, 1);
    for (const next of edges.get(r) ?? []) {
      if (state.get(next) === 1) {
        const cycle = [...trail.slice(trail.indexOf(next)), next].join(" -> ");
        report(ctx, "svs.ref.cycle", ctx.nodes.get(next)?.path ?? "", `refs form a cycle: ${cycle}`);
      } else if (!state.has(next)) visit(next, [...trail, next]);
    }
    state.set(r, 2);
  };
  for (const r of edges.keys()) if (!state.has(r)) visit(r, [r]);
}

/**
 * A dynamic node's own checks (as a document node of its kind is), on its raw params; its refs listed for the caller to
 * resolve (against the document, the set, held dynamic nodes and the same patch's). `base`: the document's archive base
 * (`ValidateOptions.base`).
 */
export function checkStreamNode(
  n: IndexedNode,
  reg: SvsRegistries,
  base?: string,
): { issues: Issue[]; refs: readonly RefUse[] } {
  const nodes = new Map([[n.ref, n]]);
  const ctx: TreeContext = {
    issues: [],
    registries: reg,
    refs: [],
    lenient: false,
    from: n.ref,
    depth: 0,
    nodes,
    base,
  };
  const schema: TreeSchema = SvsTreeSchema;
  const deep = tooDeep(n.params);
  if (deep !== null) {
    const message = `params nest deeper than ${MAX_JSON_DEPTH} objects and lists: the node cannot be read`;
    return { issues: [Issue.error("svs.params.type", `${n.path}/params${deep}`, message)], refs: [] };
  }
  checkOwn(n, schema, ctx);
  schema.nodes[n.kind].check?.(n, ctx);
  return { issues: ctx.issues, refs: ctx.refs };
}

/** Options of `nodeRefs`: read as `validate` read the document, so its memo answers. */
export interface NodeRefsOptions {
  /** A newer minor version (`checkEnvelope`'s answer). */
  readonly lenient?: boolean;
  /** The document's archive base (`ValidateOptions.base`). */
  readonly base?: string;
}

/**
 * The refs a node's params and `custom` name (seqspec:format#nodes-and-refs), each with the kinds it may name and its
 * path: every typed ref its schema records (data and targets, table target columns, a query's names with `sel`
 * expanded, a theme's `column` and `profile`, `inputs`, `items`, `follow`, `default_view`, …), as validation reads
 * them. The one source of a node's edges for a driver that assembles itself; memoized per params object with
 * validation's own checks. `n` as `indexDoc` gives it (implied refs filled); `reg`: the document's registries.
 */
export function nodeRefs(n: IndexedNode, reg: SvsRegistries, opts: NodeRefsOptions = {}): readonly RefUse[] {
  const ctx: TreeContext = {
    issues: [],
    registries: reg,
    refs: [],
    lenient: opts.lenient ?? false,
    from: n.ref,
    depth: 0,
    nodes: new Map([[n.ref, n]]),
    base: opts.base,
  };
  return ownChecks(n, SvsTreeSchema, ctx).refs;
}

/** Snapshot refs are unique; their ops are checked over a valid base (`snapshotIssues`). */
function checkSnapshots(snapshots: unknown[], ctx: TreeContext) {
  const refs = new Set<string>();
  snapshots.forEach((s, k) => {
    if (!isPlainObject(s) || typeof s.ref !== "string") return;
    if (!isRefSyntax(s.ref, false) || s.ref === "next" || s.ref === "prev") {
      const why = isRefSyntax(s.ref, false)
        ? "is reserved (a link's steps)"
        : "is not a ref: [A-Za-z0-9][A-Za-z0-9_.-]*";
      report(ctx, "svs.ref.syntax", `/snapshots/${k}/ref`, `snapshot "${s.ref}" ${why}`);
    }
    if (refs.has(s.ref))
      report(ctx, "svs.ref.duplicate", `/snapshots/${k}/ref`, `snapshot "${s.ref}" is already defined`);
    refs.add(s.ref);
  });
}
