/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The document (seqspec:viewer#the-document): base → snapshot → local patches → preview, each through `applyPatch`, so
// untouched subtrees keep their identity through every layer; patches proposed by the UI and routed by the document
// port (with the snapshot shown and whether the patch is masked); local patches tagged with the snapshot they were made
// over; the history of patches the context commits itself, kept across a host's echo of its own commits; the views of
// the shown document. In data mode every validation resolves refs against the host's `set()`.

import {
  applyPatch,
  equal,
  indexDoc,
  type Op,
  type Patch,
  rebase,
  type SeqViewSpecDoc,
  type SvsRegistries,
  snapshotDoc,
} from "@molstar/seq-core/format";
import { Issue, type KeyResolver } from "@molstar/seq-core/model";
import type { GraphResult } from "./context.js";
import type { Origin, SeqState, WritableSeqState } from "./state.js";
import { state } from "./state.js";

type Doc = SeqViewSpecDoc;
export type Routed = "commit" | "local" | "refuse";
/** What the port routes by. */
export interface RouteContext {
  /** The snapshot shown; null: the base. */
  readonly snapshot: string | null;
  /** The patch targets a node only a snapshot or a local patch adds, or a field a local patch also sets. */
  readonly masked: boolean;
}
export interface DocumentPort {
  /** Where a UI patch goes: committed to the document, kept local to this viewer, or refused. */
  route(patch: Patch, label: string, at: RouteContext): Routed;
  commit(patch: Patch, label: string): Promise<{ ok: boolean; issues: readonly Issue[] }>;
  /** True: the context keeps undo for committed patches. */
  readonly history: boolean;
}
/** A viewer-local patch and the snapshot it was made over: it applies only while that one is shown (null: the base). */
export interface LocalPatch {
  readonly patch: Patch;
  readonly snapshot: string | null;
}
export interface ProposeResult {
  readonly routed: Routed;
  readonly ok: boolean;
  readonly issues: readonly Issue[];
}
export interface DocumentManager {
  readonly base: SeqState<SeqViewSpecDoc | null>;
  readonly local: SeqState<readonly LocalPatch[]>;
  readonly preview: SeqState<Patch | null>;
  /** The snapshot layer, which `snapshots.apply` sets; null: the base. */
  readonly snapshot: SeqState<string | null>;
  readonly shown: SeqState<SeqViewSpecDoc | null>;
  readonly revision: SeqState<string | null>;
  propose(patch: Patch, label: string, origin: Origin): Promise<ProposeResult>;
  /** Applied over the current top first: an invalid preview is `ok: false` and not drawn. */
  setPreview(patch: Patch | null): Promise<GraphResult>;
  clearLocal(): void;
}
export interface HistoryManager {
  readonly canUndo: SeqState<boolean>;
  readonly canRedo: SeqState<boolean>;
  /** False when there was nothing to undo, commits are pending, or its patch is gone or no longer applies (reported). */
  undo(): Promise<boolean>;
  redo(): Promise<boolean>;
}
export interface ViewInfo {
  readonly ref: string;
  readonly label: string;
}
export interface ViewManager {
  readonly list: SeqState<readonly ViewInfo[]>;
  readonly current: SeqState<string | null>;
  /** Runtime, not a patch. */
  setCurrent(ref: string, origin?: Origin): void;
}

export interface DocumentHost {
  svs(): SvsRegistries;
  /** Data mode: what refs resolve against (the incoming set, or a resolver over the build inputs); else undefined. */
  set(): KeyResolver | undefined;
  /** An archive's document: the base its relative URIs resolve against (`ValidateOptions.base`); else undefined. */
  base(): string | undefined;
  readonly port: DocumentPort | null;
  historyLimit(): number;
  report(issue: Issue): void;
  patched(e: { patch: Patch; label: string; routed: Routed; origin: Origin }): void;
  /** Runs the graph; the result once it settles. */
  settle(issues: readonly Issue[], ok: boolean): Promise<GraphResult>;
}
interface Entry {
  readonly patch: Patch;
  readonly inverse: Patch;
  readonly label: string;
  readonly layer: "base" | "local";
  /** A local entry's patch as it is in `local`, removed by identity on undo. */
  readonly local?: LocalPatch;
}
type Applied =
  | { ok: true; issues: readonly Issue[]; inverse: Patch | null; local?: LocalPatch }
  | { ok: false; issues: readonly Issue[] };
type Check = ReturnType<typeof applyPatch>;

const FORMAT = "seq_view_spec_patch";
/** A patch over a layer above the base: its revisions are the base's business. */
const opsOf = (p: Patch): Patch => ({ format: FORMAT, ops: p.ops });
const noDocument = () => Issue.error("svs.patch.ref", "", "there is no document to patch");
/** The snapshot an `svs.snapshot.invalid` issue names (its path `/snapshots/<k>`); null: another issue. */
const invalidRef = (doc: Doc | null, i: Issue) =>
  i.code === "svs.snapshot.invalid" ? (doc?.snapshots?.[Number(i.path.split("/")[2])]?.ref ?? "") : null;
const invalidOf = (doc: Doc | null, issues: readonly Issue[]) =>
  new Set(issues.flatMap((i) => invalidRef(doc, i) ?? []));
/** The document's node refs, the implied `root` and `data` included (core's index). */
const refsOf = (doc: Doc | null): ReadonlySet<string> => new Set(doc ? indexDoc(doc).nodes.keys() : []);

export class Documents implements DocumentManager {
  readonly host: DocumentHost;
  readonly port: DocumentPort;
  readonly base = state<Doc | null>(null);
  readonly local = state<readonly LocalPatch[]>([]);
  readonly preview = state<Patch | null>(null);
  readonly snapshot = state<string | null>(null);
  readonly shown = state<Doc | null>(null);
  readonly revision = state<string | null>(null);
  readonly history: History;
  /** The layers' results, kept by their inputs' identity, so a layer change re-applies only what is above it. */
  private snap: { base: Doc | null; ref: string | null; doc: Doc | null } = { base: null, ref: null, doc: null };
  private localDocs: { lp: LocalPatch; doc: Doc }[] = [];
  /** The document under the preview. */
  private top: Doc | null = null;
  private previewDoc: { over: Doc | null; patch: Patch | null; doc: Doc | null } = {
    over: null,
    patch: null,
    doc: null,
  };
  /** Untagged local patches skipped over the snapshot shown, reported since the last snapshot change. */
  private skipped = new Set<LocalPatch>();
  /** The last check, reused by the default commit and by undo and redo. */
  private lastCheck: { doc: Doc; patch: Patch; version: number; set: unknown; r: Check } | null = null;
  /** The base's snapshots that no longer apply, from its validation or its patch's. */
  private validity: { doc: Doc | null; invalid: ReadonlySet<string> } = { doc: null, invalid: new Set() };
  /** Commits sent through a host's port since its last echo, and the base they were made on. */
  private pending: { base: Doc; patches: Patch[] } | null = null;
  private baseRefs: { doc: Doc | null; refs: ReadonlySet<string> } = { doc: null, refs: new Set() };

  constructor(host: DocumentHost) {
    this.host = host;
    this.port = host.port ?? {
      route: (_p, _l, at) => (at.snapshot !== null || at.masked ? "local" : "commit"),
      commit: async (p) => this.applyBase(p),
      history: true,
    };
    this.history = new History(this, host.historyLimit);
  }

  private apply(doc: Doc, patch: Patch) {
    return applyPatch(doc, patch, this.host.svs(), { set: this.host.set(), base: this.host.base() });
  }
  /** `applyPatch` on the base, memoized by the last call's inputs. */
  private check(doc: Doc, patch: Patch): Check {
    const c = this.lastCheck;
    const version = this.host.svs().version;
    const set = this.host.set();
    if (c?.doc === doc && c.patch === patch && c.version === version && c.set === set) return c.r;
    const r = this.apply(doc, patch);
    this.lastCheck = { doc, patch, version, set, r };
    return r;
  }

  /** The base's snapshots that no longer apply over it (`SnapshotInfo.valid`). */
  invalidSnapshots(): ReadonlySet<string> {
    return this.validity.invalid;
  }
  /** Commits sent through a host's port and not yet echoed: undo and redo wait for them. */
  get pendingCommits(): boolean {
    return this.pending !== null;
  }

  /**
   * `setDocument` and `setData`: the first is the base; later ones through `rebase`, so unchanged nodes keep their
   * identity. `issues` are the document's validation, whose snapshot checks are kept. The history survives only an
   * echo of the commits pending.
   */
  replace(next: Doc, origin: Origin | undefined, issues: readonly Issue[]): boolean {
    const base = this.base.value;
    const doc = base ? rebase(base, next, this.host.svs()) : next;
    const foreign = this.echo(doc, true);
    this.setBase(doc, origin, invalidOf(doc, issues));
    return foreign;
  }

  /** Back to no document: layers and the history cleared. */
  reset() {
    this.local.set([]);
    this.preview.set(null);
    this.snapshot.set(null);
    this.history.clear();
    this.pending = null;
    this.setBase(null, undefined, new Set());
  }

  /** A base patch; `external`: the host's (`ctx.applyPatch`), which may echo commits pending. */
  applyBase(patch: Patch, origin?: Origin, external = false): Applied {
    const base = this.base.value;
    if (!base) return { ok: false, issues: [noDocument()] };
    const r = this.check(base, patch);
    if (!r.inverse) return { ok: false, issues: r.issues };
    // A snapshot already invalid is not reported again.
    for (const i of r.issues) if (!this.validity.invalid.has(invalidRef(r.doc, i) as string)) this.host.report(i);
    if (external) this.echo(r.doc, false);
    this.setBase(r.doc, origin, invalidOf(r.doc, r.issues));
    return { ok: true, issues: r.issues, inverse: r.inverse };
  }

  /**
   * A local patch, tagged with `snapshot` (by default the one shown): checked against the base's revision, applied
   * over the local layer, all or nothing; one tagged for a snapshot not shown is kept dormant, unchecked until then.
   */
  addLocal(patch: Patch, snapshot: string | null = this.snapshot.value): Applied {
    const base = this.base.value;
    if (!base) return { ok: false, issues: [noDocument()] };
    if (patch.base != null && patch.base !== (base.revision ?? undefined)) {
      const message = `for revision "${patch.base}"; the document is at "${base.revision ?? "none"}"`;
      return { ok: false, issues: [Issue.error("svs.patch.base", "/base", message)] };
    }
    const lp: LocalPatch = { patch: opsOf(patch), snapshot };
    let r: Check | null = null;
    if (snapshot === null || snapshot === this.snapshot.value) {
      r = this.apply(this.top as Doc, lp.patch);
      if (!r.inverse) return { ok: false, issues: r.issues };
      this.localDocs.push({ lp, doc: r.doc });
    }
    this.local.set([...this.local.value, lp]);
    this.recompute(false);
    return { ok: true, issues: r?.issues ?? [], inverse: r?.inverse ?? null, local: lp };
  }

  /** Applied over the current top first; an invalid one changes nothing. */
  previewPatch(patch: Patch | null): { ok: boolean; issues: readonly Issue[] } {
    if (patch !== null) {
      if (!this.top) return { ok: false, issues: [noDocument()] };
      const r = this.apply(this.top, opsOf(patch));
      if (!r.inverse) return { ok: false, issues: r.issues };
      this.previewDoc = { over: this.top, patch, doc: r.doc };
      this.preview.set(patch);
      this.recompute(false);
      return { ok: true, issues: r.issues };
    }
    this.preview.set(null);
    this.recompute(false);
    return { ok: true, issues: [] };
  }
  setPreview(patch: Patch | null): Promise<GraphResult> {
    const r = this.previewPatch(patch);
    return this.host.settle(r.issues, r.ok);
  }
  setSnapshot(ref: string | null, origin?: Origin) {
    this.snapshot.set(ref, origin);
    this.recompute(true);
  }
  /** Every local patch, and the history's local entries (seqspec:viewer#history). */
  clearLocal() {
    this.local.set([]);
    this.history.dropLayer("local");
    this.recompute(false);
  }
  /** By identity; false: it is no longer there. */
  removeLocal(lp: LocalPatch): boolean {
    const local = this.local.value;
    const i = local.indexOf(lp);
    if (i < 0) return false;
    this.local.set([...local.slice(0, i), ...local.slice(i + 1)]);
    this.recompute(false);
    return true;
  }

  private setBase(doc: Doc | null, origin: Origin | undefined, invalid: ReadonlySet<string>) {
    this.validity = { doc, invalid };
    this.base.set(doc, origin);
    this.revision.set(doc?.revision ?? null);
    this.recompute(false);
  }

  /**
   * The shown document. A snapshot that no longer applies is left once. A local patch tagged for another snapshot is
   * dormant; one that does not apply is dropped on a base or local change when it was made over what is shown
   * (`seqv.patch.dropped`), else skipped and kept, reported once per snapshot change (`seqv.patch.skipped`).
   */
  private recompute(snapshotChange: boolean) {
    const base = this.base.value;
    let ref = this.snapshot.value;
    if (this.snap.base !== base || this.snap.ref !== ref) {
      let doc = base;
      if (base && ref !== null) {
        const r = snapshotDoc(base, ref, this.host.svs(), { set: this.host.set(), base: this.host.base() });
        if (r.issues.length > 0) {
          // Left once; an invalidity the base's validation or patch gave is surfaced there.
          if (!this.validity.invalid.has(ref)) this.host.report(r.issues[r.issues.length - 1]);
          ref = null;
          this.snapshot.set(null, "graph");
        } else doc = r.doc;
      }
      if (this.snap.ref !== ref) this.skipped.clear();
      this.snap = { base, ref, doc };
      this.localDocs = [];
    }
    let doc = this.snap.doc;
    const local = this.local.value;
    const docs: { lp: LocalPatch; doc: Doc }[] = [];
    let dropped = false;
    for (let i = 0; doc && i < local.length; i++) {
      const lp = local[i];
      const held = this.localDocs[i];
      if (held?.lp === lp && docs.length === i) {
        docs.push(held);
        doc = held.doc;
        continue;
      }
      const r = lp.snapshot === null || lp.snapshot === ref ? this.apply(doc, lp.patch) : null;
      if (r?.inverse) doc = r.doc;
      else if (r && !snapshotChange && lp.snapshot === ref) {
        const why = r.issues[0]?.message ?? "";
        this.host.report(Issue.warning("seqv.patch.dropped", `/local/${i}`, `a local patch no longer applies: ${why}`));
        dropped = true;
        continue;
      } else if (r && !this.skipped.has(lp)) {
        this.skipped.add(lp);
        const message = `a local patch does not apply over snapshot "${ref}": ${r.issues[0]?.message ?? ""}`;
        this.host.report(Issue.warning("seqv.patch.skipped", `/local/${i}`, message));
      }
      docs.push({ lp, doc });
    }
    this.localDocs = docs;
    this.top = doc;
    if (dropped) this.local.set(docs.map((d) => d.lp));
    const preview = this.preview.value;
    const p = this.previewDoc;
    if (doc && preview && (p.over !== doc || p.patch !== preview)) {
      const r = this.apply(doc, opsOf(preview));
      if (!r.inverse) for (const i of r.issues) this.host.report(i);
      this.previewDoc = { over: doc, patch: preview, doc: r.inverse ? r.doc : doc };
    }
    this.shown.set(doc && preview ? this.previewDoc.doc : doc);
  }

  /** Decision 11: a node only a snapshot or a local patch adds, or a field an active local patch also sets. */
  private masked(patch: Patch): boolean {
    const base = this.base.value;
    if (this.baseRefs.doc !== base) this.baseRefs = { doc: base, refs: refsOf(base) };
    const refs = this.baseRefs.refs;
    const shown = this.snapshot.value;
    const sets = (o: Op, ref: string, keys: readonly string[]) =>
      o.op === "set" && o.ref === ref && keys.some((k) => Object.hasOwn(o.params ?? {}, k));
    for (const op of patch.ops) {
      const target = op.op === "insert" ? op.parent : op.ref;
      if (!refs.has(target) || (op.op === "move" && !refs.has(op.parent))) return true;
      const keys = op.op === "set" ? Object.keys(op.params ?? {}) : [];
      for (const lp of this.local.value) {
        if (lp.snapshot !== null && lp.snapshot !== shown) continue;
        if (keys.length > 0 && lp.patch.ops.some((o) => sets(o, target, keys))) return true;
      }
    }
    return false;
  }

  /** The base's revision on a patch that names none. */
  private stamp(patch: Patch, base: Doc): Patch {
    return base.revision == null ? opsOf(patch) : { format: FORMAT, base: base.revision, ops: patch.ops };
  }

  /** Through the port; with a host's port it is pending until the host echoes it. */
  private async commit(sent: Patch, label: string, base: Doc) {
    const echoed = this.host.port !== null && this.port.history;
    if (echoed) {
      this.pending ??= { base, patches: [] };
      this.pending.patches.push(sent);
    }
    const r = await this.port.commit(sent, label);
    const p = this.pending;
    if (!r.ok && p?.patches.includes(sent)) {
      p.patches.splice(p.patches.indexOf(sent), 1);
      if (p.patches.length === 0) this.pending = null;
    }
    return r;
  }

  /**
   * A base from outside: the commits pending, applied in order to the base they were made on, equal to it canonically
   * is the host's echo (the history kept); anything else clears the history. A `setDocument` with nothing pending is a
   * foreign edit; a host's base patch with nothing pending keeps it (its inverses are checked on undo).
   */
  /** True: not an echo of the commits pending (a new document, its history cleared). */
  private echo(next: Doc, replaced: boolean): boolean {
    const p = this.pending;
    this.pending = null;
    if (!p) {
      if (replaced) this.history.clear();
      return replaced;
    }
    let doc: Doc | null = p.base;
    for (let i = 0; doc && i < p.patches.length; i++) {
      const r: Check = this.apply(doc, opsOf(p.patches[i]));
      doc = r.inverse ? r.doc : null;
    }
    // `next` is rebased on the base (equal nodes shared), so this mostly compares by identity: no text of either.
    if (doc && equal(doc.root, next.root)) return false;
    this.history.clear();
    return true;
  }

  /**
   * A UI patch: the port routes it, given the snapshot shown and whether it is masked. A commit is checked against the
   * base first (a stale base or a patch that does not apply is refused here, never sent) and carries the base's
   * revision; local patches stay in this viewer, tagged; both enter the history.
   */
  async propose(patch: Patch, label: string, origin: Origin): Promise<ProposeResult> {
    const routed = this.port.route(patch, label, { snapshot: this.snapshot.value, masked: this.masked(patch) });
    this.host.patched({ patch, label, routed, origin });
    const fail = (issues: readonly Issue[]) => {
      for (const i of issues) this.host.report(i);
      return { routed, ok: false, issues };
    };
    if (routed === "refuse") return fail([Issue.warning("seqv.port.refused", "", `"${label}" was refused`)]);
    if (routed === "local") {
      const r = this.addLocal(patch);
      if (!r.ok) return fail(r.issues);
      const local = r.local as LocalPatch;
      if (r.inverse) this.history.record({ patch: local.patch, inverse: r.inverse, label, layer: "local", local });
      return { routed, ok: true, issues: r.issues };
    }
    const base = this.base.value;
    if (!base) return fail([noDocument()]);
    const sent = patch.base === undefined ? this.stamp(patch, base) : patch;
    const check = this.check(base, sent);
    if (!check.inverse) return fail(check.issues);
    const r = await this.commit(sent, label, base);
    if (r.ok && this.port.history)
      this.history.record({ patch: opsOf(patch), inverse: opsOf(check.inverse), label, layer: "base" });
    return { routed, ok: r.ok, issues: r.issues };
  }

  /** Applies an entry's patch (redo) or inverse (undo) to its layer; the entry to push on the other stack, or null. */
  async step(e: Entry, undo: boolean): Promise<Entry | null> {
    const refuse = (why: string) => {
      this.host.report(Issue.warning("seqv.history.refused", "", `${undo ? "undo" : "redo"} of "${e.label}": ${why}`));
      return null;
    };
    if (e.layer === "local") {
      if (undo) return this.removeLocal(e.local as LocalPatch) ? e : refuse("its local patch is gone");
      const r = this.addLocal(e.patch, e.local?.snapshot ?? null);
      if (!r.ok) return refuse(r.issues[0]?.message ?? "it no longer applies");
      return { ...e, local: r.local, inverse: r.inverse ?? e.inverse };
    }
    const base = this.base.value;
    if (!base) return refuse("there is no document");
    const sent = this.stamp(undo ? e.inverse : e.patch, base);
    const check = this.check(base, sent);
    if (!check.inverse) return refuse(check.issues[0]?.message ?? "it no longer applies");
    const r = await this.commit(sent, `${undo ? "Undo" : "Redo"} ${e.label}`, base);
    return r.ok ? (undo ? e : { ...e, inverse: opsOf(check.inverse) }) : null;
  }

  dispose() {
    for (const s of [this.base, this.local, this.preview, this.snapshot, this.shown, this.revision])
      (s as WritableSeqState<unknown>).complete();
    this.history.dispose();
  }
}

/** Entries of patches the context committed itself (base through a port with `history`, and local). */
export class History implements HistoryManager {
  readonly canUndo = state(false);
  readonly canRedo = state(false);
  private undos: Entry[] = [];
  private redos: Entry[] = [];
  private readonly docs: Documents;
  private readonly limit: () => number;

  constructor(docs: Documents, limit: () => number) {
    this.docs = docs;
    this.limit = limit;
  }

  record(e: Entry) {
    this.undos.push(e);
    if (this.undos.length > this.limit()) this.undos.shift();
    this.redos = [];
    this.sync();
  }
  undo() {
    return this.step(this.undos, this.redos, true);
  }
  redo() {
    return this.step(this.redos, this.undos, false);
  }
  clear() {
    this.undos = [];
    this.redos = [];
    this.sync();
  }
  dropLayer(layer: Entry["layer"]) {
    this.undos = this.undos.filter((e) => e.layer !== layer);
    this.redos = this.redos.filter((e) => e.layer !== layer);
    this.sync();
  }

  private async step(from: Entry[], to: Entry[], undo: boolean) {
    if (from.length > 0 && this.docs.pendingCommits) {
      const message = `${undo ? "undo" : "redo"} waits for the host to pass back its commits`;
      this.docs.host.report(Issue.warning("seqv.history.refused", "", message));
      return false;
    }
    const e = from.pop();
    const done = e ? await this.docs.step(e, undo) : null;
    if (done) to.push(done);
    this.sync();
    return done !== null;
  }
  private sync() {
    this.canUndo.set(this.undos.length > 0);
    this.canRedo.set(this.redos.length > 0);
  }
  dispose() {
    this.canUndo.complete();
    this.canRedo.complete();
  }
}

const sameViews = (a: readonly ViewInfo[], b: readonly ViewInfo[]) =>
  a.length === b.length && a.every((v, i) => v.ref === b[i].ref && v.label === b[i].label);

/** The shown document's views; the current one is kept while it exists, else `default_view`, else the first. */
export class Views implements ViewManager {
  readonly list = state<readonly ViewInfo[]>([], sameViews);
  readonly current = state<string | null>(null);

  update(doc: Doc | null) {
    const views = (doc?.root.children ?? []).filter((c) => c.kind === "view");
    this.list.set(
      views.map((v) => ({ ref: v.ref ?? "view", label: String(v.params?.label ?? v.ref) })),
      "graph",
    );
    const list = this.list.value;
    if (list.some((v) => v.ref === this.current.value)) return;
    const preferred = list.find((v) => v.ref === doc?.root.params?.default_view);
    this.current.set((preferred ?? list[0])?.ref ?? null, "graph");
  }
  setCurrent(ref: string, origin: Origin = "api") {
    if (this.list.value.some((v) => v.ref === ref)) this.current.set(ref, origin);
  }
  dispose() {
    this.list.complete();
    this.current.complete();
  }
}
