/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Patches over a tree whose nodes all have refs (seqspec:format#patches), generic over a TreeSchema and its registries
// (implied refs, discriminators, normalized comparison), so MVS could take it as it is: ops by ref with anchors, all
// or nothing, the inverse, every untouched subtree and `params`/`custom` object shared; and `diffTree`. SVS's
// validation, revisions and snapshots are patch.ts's.

import { columnNumbers, Issue } from "@molstar/seq-core/model";
import type { DocNode, Op } from "../doc.js";
import { equal } from "../equal.js";
import { normalizeNodeParams } from "../normalize.js";
import type { SvsRegistries } from "../registry.js";
import { isPlainObject } from "../schema/field-schema.js";
import { caseOf } from "../schema/params-schema.js";
import type { TreeSchema } from "../schema/tree-schema.js";

export interface TreeOptions {
  readonly schema: TreeSchema;
  readonly reg: SvsRegistries;
}
type Obj = Record<string, unknown>;

/** A node while ops apply: its original object, reused when nothing in its subtree changed. */
interface Draft {
  readonly node: DocNode;
  readonly ref: string;
  params: Obj | undefined;
  custom: Obj | undefined;
  children: string[];
  parent: string | null;
  dirty: boolean;
}
interface Tree {
  readonly map: Map<string, Draft>;
  readonly root: string;
  readonly opts: TreeOptions;
}

const implied = (opts: TreeOptions, kind: string) =>
  Object.hasOwn(opts.schema.nodes, kind) && opts.schema.nodes[kind].ref === "implied";
/** A node's ref (`refless` and `checkOp` have seen to it that it has one, or that its kind implies it). */
const refOf = (n: DocNode) => n.ref ?? n.kind;

function add(t: Tree, n: DocNode, parent: string | null): string {
  const ref = refOf(n);
  const d: Draft = { node: n, ref, params: n.params, custom: n.custom, children: [], parent, dirty: false };
  t.map.set(ref, d);
  d.children = (n.children ?? []).map((c) => add(t, c, ref));
  return ref;
}
function draft(root: DocNode, opts: TreeOptions): Tree {
  const t = { map: new Map<string, Draft>(), root: "", opts };
  return Object.assign(t, { root: add(t, root, null) });
}

/** The node as it stands: the original object when nothing under it changed. */
function build(t: Tree, ref: string): DocNode {
  const d = t.map.get(ref) as Draft;
  const kids = d.children.map((c) => build(t, c));
  const old = d.node.children;
  if (!d.dirty && kids.length === (old?.length ?? 0) && kids.every((k, i) => k === old?.[i])) return d.node;
  const out: DocNode = { ...d.node, params: d.params, custom: d.custom, children: kids };
  if (d.params === undefined) delete out.params;
  if (d.custom === undefined) delete out.custom;
  if (!old && kids.length === 0) delete out.children;
  return out;
}

function drop(t: Tree, ref: string) {
  for (const c of t.map.get(ref)?.children ?? []) drop(t, c);
  t.map.delete(ref);
}

/** Where `after` puts a node among `parent`'s children; -1: not a child. */
function anchor(parent: Draft, after: string | undefined): number {
  if (after === undefined || after === "last") return parent.children.length;
  if (after === "first") return 0;
  const j = parent.children.indexOf(after);
  return j < 0 ? -1 : j + 1;
}
const before = (parent: Draft, ref: string) => {
  const k = parent.children.indexOf(ref);
  return k === 0 ? "first" : parent.children[k - 1];
};
function subtreeRefs(n: DocNode, opts: TreeOptions, out: string[] = []): string[] {
  if (n.ref !== undefined || implied(opts, n.kind)) out.push(refOf(n));
  for (const c of n.children ?? []) subtreeRefs(c, opts, out);
  return out;
}

/** Shallow merge: each key replaces, `null` removes; an empty result is none (`undefined`). */
function merge(cur: Obj | undefined, change: Obj, old: Obj): Obj | undefined {
  const next: Obj = { ...cur };
  for (const [k, v] of Object.entries(change)) {
    old[k] = cur?.[k] ?? null;
    if (v === null) delete next[k];
    else next[k] = v;
  }
  return Object.keys(next).length === 0 ? undefined : next;
}

/** A discriminator switch nulls the old case's params the new case lacks. */
function switched(t: Tree, kind: string, cur: Obj | undefined, change: Obj): Obj {
  const ps = Object.hasOwn(t.opts.schema.nodes, kind) ? t.opts.schema.nodes[kind].params : undefined;
  if (ps?.type !== "union" || !(ps.discriminator in change) || change[ps.discriminator] === cur?.[ps.discriminator])
    return change;
  const [from, to] = [caseOf(ps, cur?.[ps.discriminator]), caseOf(ps, change[ps.discriminator])];
  const out = { ...change };
  for (const k of Object.keys(cur ?? {})) {
    if (from && k in from.fields && !(to && k in to.fields) && !(k in change)) out[k] = null;
  }
  return out;
}

type Fail = (code: string, message: string) => null;
type Why = { code: string; message: string };

/** The path of the first node (preorder) that is not an object with a kind and a ref (implied for some kinds); null: none. */
function refless(n: unknown, opts: TreeOptions, path: string, seen = new Set<string>()): string | null {
  if (!isPlainObject(n) || typeof n.kind !== "string") return path;
  if (typeof n.ref !== "string" && !implied(opts, n.kind)) return path;
  // A ref met twice (a second `data`, or a `root` under the root) would make a node its own ancestor in the draft.
  const ref = typeof n.ref === "string" ? n.ref : n.kind;
  if (seen.has(ref)) return path;
  seen.add(ref);
  if (n.children === undefined) return null;
  if (!Array.isArray(n.children)) return `${path}/children`;
  for (let i = 0; i < n.children.length; i++) {
    const p = refless(n.children[i], opts, `${path}/children/${i}`, seen);
    if (p !== null) return p;
  }
  return null;
}

/** Why a node an op inserts is malformed (its subtree too), or null: `where` names it in the message. */
function checkNode(n: unknown, opts: TreeOptions, where: string): Why | null {
  const bad = (message: string, code = "svs.params.type") => ({ code, message: `${where}: ${message}` });
  if (!isPlainObject(n) || typeof n.kind !== "string") return bad("a node is {kind, ref, params, custom, children}");
  if (typeof n.ref !== "string" && !implied(opts, n.kind)) return bad(`this ${n.kind} needs a ref`, "svs.ref.missing");
  for (const k of ["params", "custom"]) if (n[k] !== undefined && !isPlainObject(n[k])) return bad(`${k} is an object`);
  if (n.children === undefined) return null;
  if (!Array.isArray(n.children)) return bad("children is a list of nodes");
  for (let i = 0; i < n.children.length; i++) {
    const why = checkNode(n.children[i], opts, `${where}, child ${i}`);
    if (why) return why;
  }
  return null;
}

/** Why an op is malformed (every shape is checked before any op applies), or null. */
function checkOp(op: unknown, opts: TreeOptions): Why | null {
  if (!isPlainObject(op) || !["set", "remove", "insert", "move"].includes(op.op as string)) {
    return { code: "svs.params.type", message: 'an op is {op: "set" | "remove" | "insert" | "move", ...}' };
  }
  for (const k of ["params", "custom"]) {
    if (op[k] !== undefined && !isPlainObject(op[k])) {
      return { code: "svs.params.type", message: `${op.op}: ${k} is an object (null inside it removes a key)` };
    }
  }
  for (const k of op.op === "insert" ? ["parent"] : op.op === "move" ? ["ref", "parent"] : ["ref"]) {
    if (typeof op[k] !== "string") return { code: "svs.params.type", message: `${op.op}: ${k} is a node's ref` };
  }
  if (op.after !== undefined && typeof op.after !== "string") {
    return { code: "svs.params.type", message: `${op.op}: after is a sibling's ref, "first" or "last"` };
  }
  return op.op === "insert" ? checkNode(op.node, opts, "insert node") : null;
}

/** One op on the draft; its inverse, or null (reported). `touch` is told every ref the op changes. */
function applyOp(t: Tree, op: Op, fail: Fail, touch: (ref: string) => void): Op | null {
  const get = (ref: unknown, what: string) => {
    const d = typeof ref === "string" ? t.map.get(ref) : undefined;
    if (!d) fail("svs.patch.ref", `${what} ${JSON.stringify(ref)} names no node (ops never address derived refs)`);
    return d;
  };
  const fixed = (d: Draft) =>
    d.parent === null || implied(t.opts, d.node.kind)
      ? fail("svs.patch.ref", `"${d.ref}" (${d.node.kind}) may be set, never removed or moved`)
      : d;
  if (op.op === "set") {
    const d = get(op.ref, "ref");
    if (!d) return null;
    const inv: Op & { op: "set" } = { op: "set", ref: d.ref };
    if (isPlainObject(op.params)) {
      inv.params = {};
      d.params = merge(d.params, switched(t, d.node.kind, d.params, op.params), inv.params);
    }
    if (isPlainObject(op.custom)) {
      inv.custom = {};
      d.custom = merge(d.custom, op.custom, inv.custom);
    }
    d.dirty = true;
    touch(d.ref);
    return inv;
  }
  if (op.op === "remove") {
    const d = get(op.ref, "ref");
    if (!d || !fixed(d)) return null;
    const p = t.map.get(d.parent as string) as Draft;
    const inv: Op = { op: "insert", parent: p.ref, after: before(p, d.ref), node: build(t, d.ref) };
    p.children = p.children.filter((c) => c !== d.ref);
    p.dirty = true;
    touch(p.ref);
    drop(t, d.ref);
    return inv;
  }
  const p = get(op.parent, "parent");
  if (!p) return null;
  if (op.op === "insert") {
    if (!isPlainObject(op.node) || typeof op.node.kind !== "string") return fail("svs.params.type", "insert a node");
    const refs = subtreeRefs(op.node, t.opts);
    const taken = refs.find((r, i) => t.map.has(r) || refs.indexOf(r) !== i);
    if (taken !== undefined) return fail("svs.patch.ref", `"${taken}" is already a ref of the document`);
    const at = anchor(p, op.after);
    if (at < 0) return fail("svs.patch.anchor", `"${op.after}" is not a child of "${p.ref}"`);
    p.children.splice(at, 0, add(t, op.node, p.ref));
    p.dirty = true;
    for (const r of [p.ref, ...refs]) touch(r);
    return { op: "remove", ref: refOf(op.node) };
  }
  const d = get(op.ref, "ref");
  if (!d || !fixed(d)) return null;
  for (let a: Draft | undefined = p; a; a = a.parent === null ? undefined : t.map.get(a.parent)) {
    if (a === d) return fail("svs.patch.cycle", `"${d.ref}" cannot move into its own subtree ("${p.ref}")`);
  }
  const old = t.map.get(d.parent as string) as Draft;
  const inv: Op = { op: "move", ref: d.ref, parent: old.ref, after: before(old, d.ref) };
  old.children = old.children.filter((c) => c !== d.ref);
  const at = anchor(p, op.after);
  if (at < 0) return fail("svs.patch.anchor", `"${op.after}" is not a child of "${p.ref}"`);
  p.children.splice(at, 0, d.ref);
  d.parent = p.ref;
  old.dirty = p.dirty = true;
  for (const r of [d.ref, old.ref, p.ref]) touch(r);
  return inv;
}

/**
 * Applies ops in order (intermediate states need not be valid). All or nothing: every op's shape is checked first
 * (`svs.params.type`, pathed `/ops/<k>`); then, on the first error the original root comes back with the issue. Else
 * the new root, the inverse ops and `touched`: for each ref some op changed (its parent too), the last op's index.
 */
export function patchTree(root: DocNode, ops: readonly Op[], opts: TreeOptions) {
  const touched = new Map<string, number>();
  const inverse: Op[] = [];
  const issues: Issue[] = [];
  if (!Array.isArray(ops)) issues.push(Issue.error("svs.params.type", "/ops", "ops is a list of ops"));
  else {
    ops.forEach((op, k) => {
      const why = checkOp(op, opts);
      if (why) issues.push(Issue.error(why.code, `/ops/${k}`, why.message));
    });
  }
  const bare = issues.length === 0 ? refless(root, opts, "/root") : null;
  if (bare !== null) {
    issues.push(
      Issue.error(
        "svs.ref.missing",
        "/ops",
        `a tree is patched by ref, and the node at ${bare} has none or repeats one`,
      ),
    );
  }
  if (issues.length > 0) return { root, inverse: null, issues, touched };
  const t = draft(root, opts);
  for (let k = 0; k < ops.length && issues.length === 0; k++) {
    const fail: Fail = (code, message) => {
      issues.push(Issue.error(code, `/ops/${k}`, message));
      return null;
    };
    const inv = applyOp(t, ops[k], fail, (ref) => touched.set(ref, k));
    if (inv) inverse.push(inv);
  }
  if (issues.length > 0) return { root, inverse: null, issues, touched };
  return { root: build(t, t.root), inverse: inverse.reverse(), issues, touched };
}

/** JSON with keys sorted: equal values give equal text. */
export function canonical(x: unknown): string {
  if (x === null || typeof x !== "object") return JSON.stringify(x) ?? "null";
  if (ArrayBuffer.isView(x)) return canonical(columnNumbers(x as Float64Array, null));
  if (Array.isArray(x)) return `[${x.map(canonical).join(",")}]`;
  const keys = Object.keys(x).sort();
  return `{${keys.map((k) => `${JSON.stringify(k)}:${canonical((x as Obj)[k])}`).join(",")}}`;
}

/**
 * The keys whose values differ, as `b`'s value or `null`; `na`, `nb` are compared, `b`'s raw values written. Compared
 * by `equal`, in place (CF-04: as canonical text, a 3,000² grid took 5 s).
 */
function changed(na: Obj, nb: Obj, b: Obj | undefined): Obj | null {
  const out: Obj = {};
  for (const k of new Set([...Object.keys(na), ...Object.keys(nb)])) {
    if (!equal(na[k], nb[k])) out[k] = b?.[k] ?? null;
  }
  return Object.keys(out).length > 0 ? out : null;
}

/**
 * The ops that turn `a` into `b`: nodes matched by ref and kind; a ref `a` holds under another kind removed first; `b`
 * walked in preorder, each node inserted (with its descendants not found) or moved when its parent or preceding sibling
 * differs, then set when its normalized params or `custom` differ; then every node `b` lacks removed, top-most first.
 * Ops are applied to a draft as they are made, so each is valid where it stands. Null: a tree holds a node without a
 * ref (other than the implied ones), which no op could name.
 */
export function diffTree(a: DocNode, b: DocNode, opts: TreeOptions): Op[] | null {
  if (refless(a, opts, "") !== null || refless(b, opts, "") !== null) return null;
  // The same tree: nothing to do (CF-04: the viewer rebases on every `setDocument`).
  if (a === b) return [];
  const t = draft(a, opts);
  const want = new Map<string, DocNode>();
  const index = (n: DocNode) => {
    want.set(refOf(n), n);
    n.children?.forEach(index);
  };
  index(b);
  const ops: Op[] = [];
  const run = (op: Op) => {
    applyOp(
      t,
      op,
      () => null,
      () => {},
    );
    ops.push(op);
  };
  const prune = (ref: string, keep: (d: Draft) => boolean) => {
    const d = t.map.get(ref) as Draft;
    if (!keep(d)) run({ op: "remove", ref });
    else for (const c of d.children.slice()) prune(c, keep);
  };
  prune(t.root, (d) => !want.has(d.ref) || want.get(d.ref)?.kind === d.node.kind);
  const fresh = (n: DocNode): DocNode => {
    const kids = n.children?.filter((c) => !t.map.has(refOf(c))).map(fresh);
    return !kids || (kids.length === n.children?.length && kids.every((k, i) => k === n.children?.[i]))
      ? n
      : { ...n, children: kids };
  };
  const walk = (n: DocNode, ref: string) => {
    const d = t.map.get(ref) as Draft;
    // Shared params and custom (the same objects) are equal: not normalized or compared at all (CF-04).
    if (n.params !== d.params || n.custom !== d.custom) {
      const na = normalizeNodeParams({ ...d.node, params: d.params, custom: d.custom }, opts.schema, opts.reg);
      const nb = normalizeNodeParams(n, opts.schema, opts.reg);
      const params = changed(na.params, nb.params, n.params);
      const custom = changed(na.custom ?? {}, nb.custom ?? {}, n.custom);
      if (params || custom) run({ op: "set", ref, ...(params && { params }), ...(custom && { custom }) });
    }
    (n.children ?? []).forEach((c, k) => {
      const r = refOf(c);
      const prev = k === 0 ? "first" : refOf((n.children as DocNode[])[k - 1]);
      const cur = t.map.get(r);
      if (!cur) run({ op: "insert", parent: ref, after: prev, node: fresh(c) });
      else if (cur.parent !== ref || before(d, r) !== prev) run({ op: "move", ref: r, parent: ref, after: prev });
      walk(c, r);
    });
  };
  walk(b, refOf(b));
  prune(t.root, (d) => want.has(d.ref));
  return ops;
}
