/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// SVS patches and snapshots (seqspec:format#patches, seqspec:format#snapshots): `applyPatch` (revisions, the result
// validated, snapshots re-checked), `diff`, `snapshotDoc`, and the static checks of snapshot ops `validate` runs.

import { Issue, type KeyResolver } from "@molstar/seq-core/model";
import type { DocNode, Op, Patch, SeqViewSpecDoc } from "../doc.js";
import { KernelRegistry } from "../kernel-registry.js";
import { issueOf } from "../load/issues.js";
import { createRegistries, register, type SvsRegistries } from "../registry.js";
import { isPlainObject } from "../schema/field-schema.js";
import { nodeList } from "../schema/tree-schema.js";
import { MvsRegistry } from "../tables/kernel.js";
import { SvsTreeSchema } from "../tables/svs-tree.js";
import { validate } from "../validate.js";
import { diffTree, patchTree } from "./tree-patch.js";

const PatchFormat = "seq_view_spec_patch";
const opts = (reg: SvsRegistries) => ({ schema: SvsTreeSchema, reg });
const at = (prefix: string, issues: Issue[]) => issues.map((i) => ({ ...i, path: prefix + i.path }));
const errors = (issues: Issue[]) => issues.filter((i) => i.severity === "error");

/** What a patch's result is validated with: `set`, a set or a resolver whose keys refs may name. */
export interface PatchOptions {
  set?: KeyResolver;
  /** The document's base inside an archive (`ValidateOptions.base`). */
  base?: string;
}

/** The ref of the node a path into the result names: the deepest node on its way down `children`; null: none. */
function nodeAt(root: DocNode, path: string): string | null {
  const seg = path.split("/").slice(1);
  if (seg[0] !== "root") return null;
  let n = root;
  for (let i = 1; seg[i] === "children" && n.children?.[Number(seg[i + 1])]; i += 2) n = n.children[Number(seg[i + 1])];
  return n.ref ?? n.kind;
}

/**
 * An issue of the result's validation pathed at the last op that touched its node (else `/ops`), its message naming
 * the node's ref and its path in the result, which the caller never gets.
 */
function repath(issues: Issue[], root: DocNode, touched: ReadonlyMap<string, number>): Issue[] {
  return issues.map((i) => {
    const ref = nodeAt(root, i.path);
    const k = ref === null ? undefined : touched.get(ref);
    const where = ref === null ? i.path : `${ref}: ${i.path}`;
    return { ...i, path: k === undefined ? "/ops" : `/ops/${k}`, message: `${i.message} (at ${where})` };
  });
}

/** Ops over the base tree, the result validated (its snapshots aside); the new root, or the issues. */
function over(doc: SeqViewSpecDoc, ops: readonly Op[], reg: SvsRegistries, o: PatchOptions = {}) {
  const r = patchTree(doc.root, ops, opts(reg));
  if (r.issues.length > 0) return { ...r, invalid: [] as Issue[] };
  const invalid = errors(validate({ ...doc, root: r.root, snapshots: [] }, reg, { set: o.set, base: o.base }));
  return { ...r, invalid: repath(invalid, r.root, r.touched) };
}

/**
 * A patch: all or nothing; on success the new document (sharing every untouched node, `params` and `custom`)
 * and the inverse, with a warning per snapshot that no longer applies; on failure the original and the issues.
 */
export function applyPatch(
  doc: SeqViewSpecDoc,
  patch: Patch,
  reg: SvsRegistries,
  o: PatchOptions = {},
): { doc: SeqViewSpecDoc; inverse: Patch | null; issues: Issue[] } {
  const fail = (issues: Issue[]) => ({ doc, inverse: null, issues });
  if (!isPlainObject(patch))
    return fail([Issue.error("svs.params.type", "", "a patch is {format, base, revision, ops}")]);
  if (patch.format != null && patch.format !== PatchFormat) {
    return fail([
      Issue.error("svs.format.invalid", "/format", `a patch's format is "${PatchFormat}", not "${patch.format}"`),
    ]);
  }
  if (patch.base != null && patch.base !== (doc.revision ?? undefined)) {
    return fail([
      Issue.error(
        "svs.patch.base",
        "/base",
        `for revision "${patch.base}"; the document is at "${doc.revision ?? "none"}"`,
      ),
    ]);
  }
  const r = over(doc, patch.ops, reg, o);
  if (!r.inverse || r.invalid.length > 0) return fail([...r.issues, ...r.invalid]);
  const next: SeqViewSpecDoc = { ...doc, root: r.root, revision: patch.revision };
  if (patch.revision === undefined) delete next.revision;
  const issues = (next.snapshots ?? []).flatMap((s, k) => {
    const bad = isPlainObject(s) && Array.isArray(s.ops) ? over(next, s.ops, reg, o) : null;
    if (!bad || (bad.issues.length === 0 && bad.invalid.length === 0)) return [];
    const why = [...bad.issues, ...bad.invalid][0];
    return [
      issueOf(reg, "svs.snapshot.invalid", `/snapshots/${k}`, `snapshot "${s.ref}" no longer applies: ${why.message}`),
    ];
  });
  const inverse: Patch = { format: PatchFormat, ops: r.inverse };
  if (patch.revision !== undefined) inverse.base = patch.revision;
  if (doc.revision != null) inverse.revision = doc.revision;
  return { doc: next, inverse, issues };
}

let kernel: SvsRegistries | undefined;
function kernelRegistries(): SvsRegistries {
  if (!kernel) {
    const reg = createRegistries();
    register(reg, [KernelRegistry, MvsRegistry]);
    kernel = reg;
  }
  return kernel;
}

/** The patch from `a` to `b`: `applyPatch(a, diff(a, b))` equals `b` up to normalization; null: a node has no ref. */
export function diff(a: SeqViewSpecDoc, b: SeqViewSpecDoc, reg?: SvsRegistries): Patch | null {
  const ops = diffTree(a.root, b.root, opts(reg ?? kernelRegistries()));
  if (!ops) return null;
  const out: Patch = { format: PatchFormat, ops };
  if (a.revision != null) out.base = a.revision;
  if (b.revision != null) out.revision = b.revision;
  return out;
}

/**
 * `next` sharing every node, `params` and `custom` object it has equal to `base`'s (`applyPatch(base, diff(base,
 * next))`), so a host that replaces a whole document keeps identity (seqspec:format#what-a-change-affects); `next`
 * itself when `diff` is null.
 */
export function rebase(base: SeqViewSpecDoc, next: SeqViewSpecDoc, reg: SvsRegistries): SeqViewSpecDoc {
  const ops = diffTree(base.root, next.root, opts(reg));
  const r = ops && patchTree(base.root, ops, opts(reg));
  return r && r.issues.length === 0 ? { ...next, root: r.root } : next;
}

/** The document a snapshot shows (without snapshots: theirs are over the base): its ops over the base; on
 * failure the base and the issues. */
export function snapshotDoc(
  doc: SeqViewSpecDoc,
  snapshot: string,
  reg: SvsRegistries,
  o: PatchOptions = {},
): { doc: SeqViewSpecDoc; issues: Issue[] } {
  const k = (doc.snapshots ?? []).findIndex((s) => s != null && s.ref === snapshot);
  if (k < 0) return { doc, issues: [Issue.error("svs.ref.unresolved", "/snapshots", `no snapshot "${snapshot}"`)] };
  const r = over(doc, doc.snapshots?.[k].ops ?? [], reg, o);
  const issues = at(`/snapshots/${k}`, [...r.issues, ...r.invalid]);
  if (issues.length === 0) return { doc: { ...doc, root: r.root, snapshots: [] }, issues };
  const message = `snapshot "${snapshot}" does not apply: ${issues[0].message}`;
  return { doc, issues: [...issues, issueOf(reg, "svs.snapshot.invalid", `/snapshots/${k}`, message)] };
}

/**
 * The static checks of snapshots (seqspec:format#validation): each one's ops apply to the base (`svs.patch.*`, under
 * the snapshot) and give a valid document (else `svs.snapshot.invalid`); a ref inserted by several has one kind and one
 * parent.
 */
export function snapshotIssues(doc: SeqViewSpecDoc, reg: SvsRegistries, o: PatchOptions = {}): Issue[] {
  const issues: Issue[] = [];
  const inserted = new Map<string, string>();
  const note = (n: DocNode, parent: string, path: string) => {
    const where = `${n.kind} under ${parent}`;
    if (typeof n.ref === "string" && (inserted.get(n.ref) ?? where) !== where) {
      issues.push(Issue.error("svs.snapshot.ref", path, `"${n.ref}" is inserted as ${inserted.get(n.ref)} elsewhere`));
    }
    if (typeof n.ref === "string") inserted.set(n.ref, where);
    for (const c of nodeList(n.children)) if (isPlainObject(c)) note(c, String(n.ref), path);
  };
  (doc.snapshots ?? []).forEach((s, k) => {
    if (!isPlainObject(s) || !Array.isArray(s.ops)) return;
    for (const t of ["duration_ms", "transition_ms"] as const) {
      if (typeof s[t] === "number" && s[t] < 0) {
        issues.push(Issue.error("svs.params.type", `/snapshots/${k}/${t}`, `${t} is a duration, >= 0`));
      }
    }
    const r = over(doc, s.ops, reg, o);
    for (const i of at(`/snapshots/${k}`, r.issues)) issues.push(i);
    if (r.invalid.length > 0) {
      const message = `snapshot "${s.ref}" gives an invalid document: ${r.invalid[0].message}`;
      issues.push(issueOf(reg, "svs.snapshot.invalid", `/snapshots/${k}`, message));
    }
    s.ops.forEach((op, j) => {
      if (isPlainObject(op) && op.op === "insert" && isPlainObject(op.node)) {
        note(op.node as DocNode, String(op.parent), `/snapshots/${k}/ops/${j}`);
      }
    });
  });
  return issues;
}
