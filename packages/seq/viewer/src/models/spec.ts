/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The spec's model (_spec/README.md "SpecModel"): revisions are the truth, patches derived. It holds the head commit
// (the spec's immutable root, its revision, its parents), applies `repr` patches with the language, its slots and the
// `seq.lang` checks of a change (all or nothing; a stale one refused, or rebased when the host opts in), and makes a
// commit per applied change, keeping the recent roots by revision. Over the root it shows a preview (a patch, applied
// again when the root changes) and over that a story step (a stored patch: navigation, never an edit). It takes part
// in the host's history through Seq*'s history port and knows nothing of the host's stack.

import { batch, derive, type ReadonlyState, ReactiveModel, State } from '@molstar/core/reactive';
import {
  applyPatch,
  type ApplyOptions,
  type Changes,
  type Conflict,
  merge3,
  type Patch,
  type PatchResult,
  patchFromNode,
} from '@molstar/core/repr/patch';
import { hashNode, type Issue, issueOf, type Node, printAddress, resolveAddress } from '@molstar/core/repr/tree';
import { SeqLanguage, SeqPatchValidator, SeqValidateOptions } from '@molstar/seqx-core/lang';
import { createResolver, type ResolveOptions, type StaticModel } from '@molstar/seqx-core/resolve';
import { DefaultViewerConfig, type ViewerConfig } from '@molstar/seqx-viewer/config';
import { ViewerIssues } from '@molstar/seqx-viewer/issues';
import type {
  EditContext,
  EditOrigin,
  SeqHistory,
  SeqHistoryParticipant,
  SeqHistoryPort,
  SpecCommit,
} from '@molstar/seqx-viewer/ports';

/** A whole new root, or a function of the current one (a document loaded, a host's own setting). */
export type SpecChange = Node | ((root: Node) => Node);

/** Why `apply` refused a whole root. */
export type SpecRefusal = 'stale-revision';

export interface SpecApplyOptions {
  /** The revision the change was made against: refused when it is not the current one (EC-81). */
  readonly against?: string;
  /** The edit it is part of (`edit`'s context); without it, with a history attached, a change outside it. */
  readonly context?: EditContext;
}

export interface SpecPatchOptions {
  /** The edit it is part of (`edit` passes it); without it, with a history attached, a change outside the history. */
  readonly context?: EditContext;
  /** Refuse a patch without a `base` (an agent's or a remote channel's). */
  readonly requireBase?: boolean;
  readonly allowDangling?: boolean;
  readonly allowRebind?: boolean;
  /** A stale patch rebased through `merge3` from its base root (`edit.rebaseStale` by default: off). */
  readonly rebase?: boolean;
}

export type SpecPatchResult =
  | {
      readonly ok: true;
      readonly changed: boolean;
      /** The head after it (the same one when nothing changed). */
      readonly commit: SpecCommit;
      readonly changes: Changes;
      /** Warnings and infos: the core's, an edit masked by the step shown. */
      readonly issues: readonly Issue[];
      /** It was stale, and rebased through `merge3`. */
      readonly rebased?: true;
    }
  | {
      readonly ok: false;
      readonly issues: readonly Issue[];
      /** A rebase's conflicts, with its partial merge (ours and theirs' independent changes). */
      readonly conflicts?: readonly Conflict[];
      readonly partial?: Node;
    };

/** What an edit did: every patch it applied, whether any changed the spec, whether all applied. */
export interface EditOutcome {
  readonly ok: boolean;
  readonly changed: boolean;
  readonly results: readonly SpecPatchResult[];
  readonly issues: readonly Issue[];
}

/** An edit in progress: what it applies is one step of the host's history. */
export interface SpecEdit {
  readonly context: EditContext | undefined;
  apply(patch: Patch, options?: Omit<SpecPatchOptions, 'context'>): SpecPatchResult;
  /** Ends the preview in its one patch (applied in full); undefined when there is none. */
  commitPreview(): SpecPatchResult | undefined;
}

/** The core's options for a patch of the spec: the language, its slots, the `seq.lang` checks of a change. */
export const SpecPatchApply: ApplyOptions = {
  language: SeqLanguage,
  validate: SeqValidateOptions,
  validators: [SeqPatchValidator],
};

/** A preview at interactive rates: normalized and validated, not the language's own checks (they run at commit). */
const PreviewApply: ApplyOptions = { language: SeqLanguage, validate: SeqValidateOptions };

/** A story step over what is shown: the language and its slots, refs followed, nothing else refused. */
const StepApply: ApplyOptions = { language: SeqLanguage, validate: SeqValidateOptions, allowDangling: true };

const NoChanges: Changes = { made: [], params: [], children: [], removed: [], moved: [] };
const WholeRoot: Changes = { made: [[]], params: [], children: [], removed: [], moved: [] };

/** The step's own ops, by address and key, of `set` and `unset` (what masks an edit). */
function stepKeys(patch: Patch): Set<string> {
  const out = new Set<string>();
  for (const op of patch.ops)
    if (op.op === 'set' || op.op === 'unset') out.add(`${printAddress(op.at)}|${op.where ?? 'params'}|${op.key}`);
  return out;
}

/** The same issues in the same order, by code, address and message: a composition that says the same tells nothing. */
function sameIssues(a: readonly Issue[], b: readonly Issue[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++)
    if (a[i].code !== b[i].code || a[i].address !== b[i].address || a[i].message !== b[i].message) return false;
  return true;
}

const commitOf = (root: Node, parents: readonly string[]): SpecCommit => ({
  root,
  revision: hashNode(root),
  parents,
});

/** A shown story step: its id, its patch (undefined when its node no longer holds one), whether it applies now. */
export interface ShownStep {
  readonly id: string;
  readonly address: string;
  readonly label: string;
  readonly applies: boolean;
}

export class SpecModel extends ReactiveModel implements SeqHistoryParticipant {
  private readonly _head: State<SpecCommit>;
  private readonly _preview = new State<Patch | undefined>(undefined);
  private readonly _stepId = new State<string | undefined>(undefined);
  private readonly _step = new State<ShownStep | undefined>(undefined);
  private readonly _edited: State<Node>;
  private readonly _shown: State<Node>;
  private readonly _issues = new State<readonly Issue[]>([], sameIssues);
  private readonly _changed: State<SpecCommit>;
  private readonly _recent = new Map<string, Node>();
  private readonly _beforeCommit: ((changes: Changes) => void)[] = [];
  private _history: SeqHistory | undefined;
  private _contexts = 0;

  readonly config: ViewerConfig;
  /** The head commit: the root, its revision, its parents. */
  readonly head: ReadonlyState<SpecCommit>;
  /** The committed root (the base an edit patches): what the history checkpoints and a save writes. */
  readonly root: ReadonlyState<Node>;
  /** The root's hash: equal roots, equal revisions. */
  readonly revision: ReadonlyState<string>;
  /** A gesture's patch over the root (a drag resizing a track, a theme tried on hover): shown, never saved or in history. */
  readonly preview: ReadonlyState<Patch | undefined>;
  /** The root with the preview applied: what an edit made now would patch, plus the gesture in progress. */
  readonly edited: ReadonlyState<Node>;
  /** The story step shown over it, and whether it applies. */
  readonly step: ReadonlyState<ShownStep | undefined>;
  /** What is shown: `step(applyPatch(root, preview))`. */
  readonly shown: ReadonlyState<Node>;
  /** The static model of what is shown, by a resolver that keeps its entries (by hash) and paths between specs. */
  readonly model: ReadonlyState<StaticModel>;
  /** Why the preview was dropped, the step does not apply, an edit was masked: the last composition's, or edit's. */
  readonly issues: ReadonlyState<readonly Issue[]>;
  /** The last commit an edit made (the participant's face: a restore makes none). */
  readonly changed: ReadonlyState<SpecCommit>;

  constructor(root: Node, options: ResolveOptions = {}, config: ViewerConfig = DefaultViewerConfig) {
    super();
    this.config = config;
    const head = commitOf(root, []);
    this._remember(head);
    this._head = new State(head);
    this._changed = new State(head);
    this._edited = new State(root);
    this._shown = new State(root);
    this.head = this._head.asReadonly;
    this.root = derive(this._head, (h) => h.root);
    this.revision = derive(this._head, (h) => h.revision);
    this.preview = this._preview.asReadonly;
    this.edited = this._edited.asReadonly;
    this.step = this._step.asReadonly;
    this.shown = this._shown.asReadonly;
    this.issues = this._issues.asReadonly;
    this.changed = this._changed.asReadonly;
    const resolver = createResolver(options);
    this.model = derive(this.shown, (s) => resolver.resolve(s));
  }

  // ---- revisions ----

  /** A root by revision: the recent ones (`revisions.keepRecent`), then the history's (`root(revision)`). */
  rootOf(revision: string): Node | undefined {
    return this._recent.get(revision) ?? this._history?.root?.(revision);
  }

  private _remember(commit: SpecCommit): void {
    this._recent.delete(commit.revision);
    this._recent.set(commit.revision, commit.root);
    const keep = Math.max(1, this.config.revisions.keepRecent);
    while (this._recent.size > keep) this._recent.delete(this._recent.keys().next().value!);
  }

  /**
   * Tells `f` of each commit's changes before the new root is set (the session's addresses follow a move, as the
   * spec's refs do: `AddressModel.followMoves`); returns the detach.
   */
  beforeCommit(f: (changes: Changes) => void): () => void {
    this._beforeCommit.push(f);
    return () => {
      const i = this._beforeCommit.indexOf(f);
      if (i >= 0) this._beforeCommit.splice(i, 1);
    };
  }

  /** Makes a commit: the session follows, the head moves, the preview and the step apply again over it. */
  private _commit(
    root: Node,
    changes: Changes,
    context: EditContext | undefined,
    extra: readonly Issue[] = [],
  ): SpecCommit {
    const before = this._head.value;
    const commit = commitOf(root, [before.revision]);
    for (const f of this._beforeCommit) f(changes);
    this._remember(commit);
    batch(() => {
      this._head.set(commit);
      this._compose(extra);
      this._changed.set(commit);
    });
    if (this._history && !this._history.isActive(context))
      this._history.outside('The spec was changed outside the history');
    return commit;
  }

  // ---- patches ----

  /**
   * Applies a patch to the root (_spec/README.md "SpecModel"): with the language, its slots and the `seq.lang` checks
   * of what it changed; all or nothing; one commit when it changed something. A stale one (its `base` not the current
   * revision) is refused, or rebased through `merge3` when `rebase` (or `edit.rebaseStale`) says so. While a story step
   * is shown, the patch edits the root under it: one naming a node only the step makes is refused, one that sets a key
   * the step sets too is reported masked. Run it inside `edit` for one step of the host's history.
   */
  applyPatch(patch: Patch, options: SpecPatchOptions = {}): SpecPatchResult {
    const head = this._head.value;
    const stale = patch.base !== undefined && patch.base !== head.revision && patch.result !== head.revision;
    if (stale && (options.rebase ?? this.config.edit.rebaseStale)) return this._rebase(patch, options);
    const r = applyPatch(head.root, patch, this._options(options));
    if (!r.ok) return { ok: false, issues: this._stepHints(r.issues) };
    if (!r.changed) return { ok: true, changed: false, commit: head, changes: r.changes, issues: r.issues };
    const issues = [...r.issues, ...this._masked(r.changes)];
    const commit = this._commit(r.root, r.changes, options.context, issues);
    return { ok: true, changed: true, commit, changes: r.changes, issues };
  }

  private _options(options: SpecPatchOptions): ApplyOptions {
    return {
      ...SpecPatchApply,
      ...(options.requireBase ? { requireBase: true } : {}),
      ...(options.allowDangling ? { allowDangling: true } : {}),
      ...(options.allowRebind ? { allowRebind: true } : {}),
    };
  }

  /**
   * A stale patch rebased (merge.md's rule, the host's opt-in): `merge3(base, current, applyPatch(base, p))` with the
   * same options on both, the base root from the recent roots or the history's. Conflicts refuse it, with the partial
   * merge for the host to show or rescue.
   */
  private _rebase(patch: Patch, options: SpecPatchOptions): SpecPatchResult {
    const head = this._head.value;
    const base = this.rootOf(patch.base!);
    if (!base)
      return {
        ok: false,
        issues: [
          issueOf(ViewerIssues, 'seq.viewer.rebase.base', 'The revision the patch was made against is not known', {
            hint: ViewerIssues.codes['seq.viewer.rebase.base'].hint,
          }),
        ],
      };
    const o = this._options(options);
    const theirs = applyPatch(base, patch, o);
    if (!theirs.ok) return { ok: false, issues: theirs.issues };
    const m = merge3(base, head.root, theirs.root, o);
    if (!m.ok)
      return {
        ok: false,
        issues: [
          issueOf(ViewerIssues, 'seq.viewer.rebase.conflict', `${m.conflicts.length} conflict(s) with the spec`, {
            hint: ViewerIssues.codes['seq.viewer.rebase.conflict'].hint,
          }),
          ...m.issues,
        ],
        conflicts: m.conflicts,
        partial: m.partial,
      };
    if (m.root === head.root || hashNode(m.root) === head.revision)
      return { ok: true, changed: false, commit: head, changes: NoChanges, issues: m.issues, rebased: true };
    const issues = [...m.issues, ...this._masked(m.changes)];
    const commit = this._commit(m.root, m.changes, options.context, issues);
    return { ok: true, changed: true, commit, changes: m.changes, issues, rebased: true };
  }

  /** A target refused while a step is shown, that the step makes: the hint names the step. */
  private _stepHints(issues: readonly Issue[]): readonly Issue[] {
    const step = this._step.value;
    if (!step) return issues;
    return issues.map((i) =>
      i.code === 'repr.patch.target'
        ? { ...i, hint: `It may be only the story step ${step.address} that makes it: leave the step to edit it` }
        : i,
    );
  }

  /** An edit under a shown step that sets a key the step sets too: it lands in the root and shows nothing (an info). */
  private _masked(changes: Changes): Issue[] {
    const step = this._step.value;
    const patch = step && this._stepPatch(step.id);
    if (!patch) return [];
    const keys = stepKeys(patch);
    const out: Issue[] = [];
    for (const c of changes.params) {
      const at = printAddress(c.at);
      for (const k of c.keys)
        if (keys.has(`${at}|${c.where}|${k}`))
          out.push(
            issueOf(ViewerIssues, 'seq.viewer.edit.masked', `The step ${step!.address} sets ${k} too`, {
              address: at,
              hint: ViewerIssues.codes['seq.viewer.edit.masked'].hint,
            }),
          );
    }
    return out;
  }

  // ---- whole roots ----

  /**
   * Sets a whole new root (a document loaded, a host's own setting): one commit, its changes the whole root. An equal
   * tree is no change; one made against another revision is refused (EC-81). A preview is dropped.
   */
  apply(change: SpecChange, options: SpecApplyOptions = {}): boolean | SpecRefusal {
    const head = this._head.value;
    if (options.against !== undefined && options.against !== head.revision) return 'stale-revision';
    const next = typeof change === 'function' ? change(head.root) : change;
    this._preview.set(undefined);
    if (next === head.root || hashNode(next) === head.revision) {
      this._compose();
      return false;
    }
    this._commit(next, WholeRoot, options.context);
    return true;
  }

  // ---- edits and the history ----

  /** Binds the host's history (the port's factory, given this participant); returns the history, `detach` to unbind. */
  attachHistory(port: SeqHistoryPort): SeqHistory {
    const history = port(this);
    this._history = history;
    this.registerPermanentDispose(() => this.detachHistory(history));
    return history;
  }

  /** Unbinds a history (disposed unless the host owns it). */
  detachHistory(history: SeqHistory): void {
    if (this._history !== history) return;
    this._history = undefined;
    if (history.disposeWithViewer !== false) history.dispose?.();
  }

  /** The history bound, if any. */
  get history(): SeqHistory | undefined {
    return this._history;
  }

  /**
   * Runs an edit: everything `work` applies is one step of the host's history (with its label and origin), restored
   * there when a patch of it is refused. Without a history it runs at once, and a refusal moves the root back. Returns
   * the outcome (a promise when a history runs it).
   */
  edit(
    label: string,
    work: (edit: SpecEdit) => void,
    options: { readonly origin?: EditOrigin } = {},
  ): EditOutcome | Promise<EditOutcome> {
    const run = (context: EditContext | undefined): EditOutcome => {
      const results: SpecPatchResult[] = [];
      const edit: SpecEdit = {
        context,
        apply: (patch, o = {}) => {
          const r = this.applyPatch(patch, { ...o, context });
          results.push(r);
          return r;
        },
        commitPreview: () => {
          const r = this.commitPreview({ context });
          if (r) results.push(r);
          return r;
        },
      };
      work(edit);
      return {
        ok: results.every((r) => r.ok),
        changed: results.some((r) => r.ok && r.changed),
        results,
        issues: results.flatMap((r) => r.issues),
      };
    };
    const history = this._history;
    if (!history) {
      const before = this._head.value;
      const outcome = run({ id: ++this._contexts });
      // all or nothing: a refused patch takes back what the edit applied before it
      if (!outcome.ok && this._head.value !== before) this._move(before);
      return outcome;
    }
    return history.run({ label, kind: 'patch', origin: options.origin ?? { kind: 'person' } }, run, {
      changed: (o) => o.ok && o.changed,
      failed: (o) => !o.ok,
    });
  }

  /** Applies patches as one edit (one step), in order; all or nothing. */
  patch(
    label: string,
    patches: Patch | readonly Patch[],
    options: { readonly origin?: EditOrigin } & Omit<SpecPatchOptions, 'context'> = {},
  ): EditOutcome | Promise<EditOutcome> {
    const { origin, ...o } = options;
    const list = Array.isArray(patches) ? patches : [patches as Patch];
    return this.edit(
      label,
      (e) => {
        for (const p of list) if (!e.apply(p, o).ok) break;
      },
      origin ? { origin } : {},
    );
  }

  /** Opens an agent's turn (`origin.turn`) on the history: its edits join one step; returns its end. */
  beginTurn(label: string, origin: EditOrigin): () => void {
    return this._history?.beginTurn?.(label, origin) ?? (() => {});
  }

  // ---- previews ----

  /**
   * Shows a gesture's patch over the root without committing it (_spec/README.md "SpecModel", Previews): normalized
   * and validated, not the language's own checks. A later call replaces it (it is made from the root, not from the
   * last preview). Refused, nothing changes. When the root changes, it is applied again over the new root, and
   * dropped when it no longer applies.
   */
  setPreview(patch: Patch): PatchResult {
    const p: Patch = { ops: patch.ops };
    const r = applyPatch(this._head.value.root, p, PreviewApply);
    if (!r.ok) return r;
    batch(() => {
      this._preview.set(p);
      this._compose();
    });
    return r;
  }

  /** Drops the preview (Esc): the root is shown again, nothing recorded. */
  cancelPreview(): void {
    if (!this._preview.value) return;
    batch(() => {
      this._preview.set(undefined);
      this._compose();
    });
  }

  /**
   * Ends the preview in its one patch, applied in full (the `seq.lang` checks too): one commit, or nothing when
   * refused. Run it inside `edit` for one step. Undefined when there was none.
   */
  commitPreview(options: { readonly context?: EditContext } = {}): SpecPatchResult | undefined {
    const p = this._preview.value;
    if (!p) return undefined;
    this._preview.set(undefined);
    const r = this.applyPatch(p, options);
    if (!r.ok || !r.changed) this._compose();
    return r;
  }

  // ---- story steps ----

  /** The stored patch of a step, by id; undefined when the spec has none such. */
  private _stepPatch(id: string): Patch | undefined {
    const node = resolveAddress(this._head.value.root, ['snapshots', id]);
    const p = node?.kind === 'snapshot' ? node.params?.patch : undefined;
    return p && typeof p === 'object' && !Array.isArray(p) && 'kind' in p ? patchFromNode(p as Node).patch : undefined;
  }

  /**
   * Shows a story step over the root and the preview (navigation: no commit, no history step). Undefined when it was
   * entered (whether it applies is `step.applies`), else the issue why not.
   */
  enterStep(id: string): Issue | undefined {
    const node = resolveAddress(this._head.value.root, ['snapshots', id]);
    if (node?.kind !== 'snapshot')
      return issueOf(ViewerIssues, 'seq.viewer.step.unknown', `No story step ${id}`, {
        address: printAddress(['snapshots', id]),
        hint: ViewerIssues.codes['seq.viewer.step.unknown'].hint,
      });
    batch(() => {
      this._stepId.set(id);
      this._compose();
    });
    return undefined;
  }

  /** Leaves the step shown: the root (and the preview) shown again. */
  leaveStep(): void {
    if (this._stepId.value === undefined) return;
    batch(() => {
      this._stepId.set(undefined);
      this._compose();
    });
  }

  /**
   * What is shown, again: the preview over the root (dropped when it no longer applies), the step over that (shown as
   * not applying when it does not); the issues of this composition, after `extra` (an edit's).
   */
  private _compose(extra: readonly Issue[] = []): void {
    const root = this._head.value.root;
    const issues: Issue[] = [...extra];
    let edited = root;
    const preview = this._preview.value;
    if (preview) {
      const r = applyPatch(root, preview, PreviewApply);
      if (r.ok) edited = r.root;
      else {
        this._preview.set(undefined);
        if (this.config.preview.dropped === 'issue')
          issues.push(
            issueOf(ViewerIssues, 'seq.viewer.preview.dropped', 'The preview no longer applies over the spec', {
              hint: ViewerIssues.codes['seq.viewer.preview.dropped'].hint,
            }),
          );
      }
    }
    let shown = edited;
    let step: ShownStep | undefined;
    const id = this._stepId.value;
    if (id !== undefined) {
      const node = resolveAddress(root, ['snapshots', id]);
      const address = printAddress(['snapshots', id]);
      if (node?.kind !== 'snapshot') this._stepId.set(undefined);
      else {
        const patch = this._stepPatch(id);
        const r = patch ? applyPatch(edited, patch, StepApply) : undefined;
        const label = typeof node.params?.label === 'string' ? node.params.label : id;
        if (r?.ok) shown = r.root;
        else
          issues.push(
            issueOf(ViewerIssues, 'seq.viewer.step.stale', `The step ${label} no longer applies`, {
              address,
              hint: ViewerIssues.codes['seq.viewer.step.stale'].hint,
            }),
          );
        step = { id, address, label, applies: !!r?.ok };
      }
    }
    const prev = this._step.value;
    if (
      prev?.id !== step?.id ||
      prev?.applies !== step?.applies ||
      prev?.label !== step?.label ||
      prev?.address !== step?.address
    )
      this._step.set(step);
    this._edited.set(edited);
    this._shown.set(shown);
    this._issues.set(issues);
  }

  // ---- the participant ----

  checkpoint(): SpecCommit {
    return this._head.value;
  }

  /** Moves the root back to a commit (an undo, a redo): no new commit; a preview dropped, the step shown kept. */
  async restore(commit: SpecCommit): Promise<void> {
    this._preview.set(undefined);
    this._move(commit);
  }

  private _move(commit: SpecCommit): void {
    this._remember(commit);
    batch(() => {
      this._head.set(commit);
      this._compose();
    });
  }
}
