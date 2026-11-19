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
} from '@molstar/core/repr/patch';
import { hashNode, type Issue, issueOf, type Node, printAddress, resolveAddress } from '@molstar/core/repr/tree';
import { checkSpecBounds, SeqLanguage, SeqValidateOptions } from '@molstar/seqx-core/lang';
// the change validators are their own entry: reading a spec carries no editing (core _spec, "Checks of a change")
import { SeqPatchValidator, seqBoundsValidator, seqPatchValidator } from '@molstar/seqx-core/lang/changes';
import { createResolver, type Lengths, type ResolveOptions, type StaticModel } from '@molstar/seqx-core/resolve';
import type { StepState } from '@molstar/seqx-core/story';
import type { SeqStorySteps } from '@molstar/seqx-viewer/story';
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
export type SpecRefusal = 'stale-revision' | 'out-of-bounds' | 'frozen-step';

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

/** Who edits, and the host's run an edit is part of (a host's step whose work runs a Seq* command). */
export interface EditOptions {
  readonly origin?: EditOrigin;
  /** The host's run in progress: the edit joins it, no step (and no wait) of its own. */
  readonly within?: EditContext;
}

/** Patches made from the root and revision an edit runs on (`SpecModel.patch`). */
export type PatchMaker = (root: Node, revision: string) => Patch | readonly Patch[];

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
// (the bounds on nesting and size run on a preview and a step too: `_previewOptions`, `_stepOptions`)

/** A story step over what is shown: the language and its slots, refs followed, nothing else refused. */
const StepApply: ApplyOptions = { language: SeqLanguage, validate: SeqValidateOptions, allowDangling: true };
const StepIssueCodes = new Set(['seq.story.step.from', 'seq.story.step.form']);

type Move = Changes['moved'][number];

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
  /** A patch (edits patch the spec under it) or a full state (edits are refused: edit the step). */
  readonly form: 'patch' | 'full';
  /** Its state: the core's `stateOf` over the root and the preview (its view and timing among it). */
  readonly state: StepState;
}

/** No source-backed system's length: a spec model of its own (a test's) resolves without. */
const NoLengths: ReadonlyState<Lengths> = new State<Lengths>(new Map()).asReadonly;

export class SpecModel extends ReactiveModel implements SeqHistoryParticipant {
  private readonly _head: State<SpecCommit>;
  private readonly _boundsIssues = new State<readonly Issue[]>([]);
  private readonly _preview = new State<Patch | undefined>(undefined);
  private readonly _stepId = new State<string | undefined>(undefined);
  private readonly _step = new State<ShownStep | undefined>(undefined);
  private readonly _edited: State<Node>;
  private readonly _shown: State<Node>;
  private readonly _issues = new State<readonly Issue[]>([], sameIssues);
  private readonly _changed: State<SpecCommit>;
  private readonly _loaded: State<boolean>;
  private readonly _recent = new Map<string, Node>();
  /** Each recent commit's parent and moves, by revision: what the session is told when the head moves between them. */
  private readonly _lineage = new Map<string, { readonly parent: string; readonly moved: readonly Move[] }>();
  /** Edits running: their commits are told to `changed` when the edit ends well, never one it took back. */
  private _editing = 0;
  private readonly _beforeCommit: ((changes: Changes) => void)[] = [];
  private _history: SeqHistory | undefined;
  private _contexts = 0;
  /** The steps' states, by the revision shown under them and the step (_spec/README.md "SpecModel", Story steps). */
  /** The story steps' reading and its cache (`attachStory`): none, no step is entered. */
  private _stories: ReturnType<SeqStorySteps['createStoryCache']> | undefined;
  private _story: SeqStorySteps | undefined;
  /** What `_compose` showed: the frame slot shows over it, and leaves back to it. */
  private _composed: Node;
  private _frame: Node | undefined;
  /** The preview's last result over a root: `_compose` reuses it (one apply a move). */
  private _previewApplied: { readonly root: Node; readonly patch: Patch; readonly result: Node } | undefined;

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
  /** What is shown: `frame(step(applyPatch(root, preview)))`, the frame an animation's, transient. */
  readonly shown: ReadonlyState<Node>;
  /** The static model of what is shown, by a resolver that keeps its entries (by hash) and paths between specs. */
  readonly model: ReadonlyState<StaticModel>;
  /** Why the preview was dropped, the step does not apply, an edit was masked: the last composition's, or edit's. */
  readonly issues: ReadonlyState<readonly Issue[]>;
  /** The issues of the last whole root `apply` refused for its bounds (`seq.lang.bound.*`); empty after one it took. */
  readonly boundsIssues: ReadonlyState<readonly Issue[]>;
  /** The last commit an edit made (the participant's face: a restore makes none). */
  readonly changed: ReadonlyState<SpecCommit>;
  /**
   * Whether a spec was ever set: false for the empty spec a viewer starts with, true after the first `apply` (a
   * document loaded, even an equal one). Until then no address names anything, so an address waits (_spec/README.md
   * "AddressModel").
   */
  readonly loaded: ReadonlyState<boolean>;

  /**
   * `lengths`: the lengths the opened sources' coverages give source-backed systems (the viewer's, from its store): the
   * model is resolved again with them (_spec/README.md "SpecModel").
   */
  constructor(
    root: Node,
    options: ResolveOptions = {},
    config: ViewerConfig = DefaultViewerConfig,
    loaded = true,
    lengths: ReadonlyState<Lengths> = NoLengths,
  ) {
    super();
    this.config = config;
    this._loaded = new State(loaded);
    this.loaded = this._loaded.asReadonly;
    const head = commitOf(root, []);
    this._remember(head);
    this._head = new State(head);
    this._changed = new State(head);
    this._edited = new State(root);
    this._shown = new State(root);
    this._composed = root;
    this.head = this._head.asReadonly;
    this.root = derive(this._head, (h) => h.root);
    this.revision = derive(this._head, (h) => h.revision);
    this.preview = this._preview.asReadonly;
    this.edited = this._edited.asReadonly;
    this.step = this._step.asReadonly;
    this.shown = this._shown.asReadonly;
    this.issues = this._issues.asReadonly;
    this.changed = this._changed.asReadonly;
    this.boundsIssues = this._boundsIssues.asReadonly;
    const resolver = createResolver(options);
    this.model = derive([this.shown, lengths], (s, l) => resolver.resolve(s, l));
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
    while (this._recent.size > keep) {
      const old = this._recent.keys().next().value!;
      this._recent.delete(old);
      this._lineage.delete(old);
    }
  }

  /**
   * The moves from one commit to another along the recent lineage: a commit's own going forward (a redo), their
   * inverses going back (an undo, a rollback); none when the two are not on one line of recent commits.
   */
  private _movesBetween(from: SpecCommit, to: SpecCommit): readonly Move[] {
    if (from.revision === to.revision) return [];
    const up = (start: string, stop: string) => {
      const path: (readonly Move[])[] = [];
      for (let r = start, n = 0; n <= this._lineage.size; n++) {
        if (r === stop) return path;
        const l = this._lineage.get(r);
        if (!l) return undefined;
        path.push(l.moved);
        r = l.parent;
      }
      return undefined;
    };
    const back = up(from.revision, to.revision);
    if (back) return back.flatMap((m) => [...m].reverse().map((x) => ({ from: x.to, to: x.from })));
    const forward = up(to.revision, from.revision);
    return forward ? forward.reverse().flat() : [];
  }

  /**
   * Tells `f` of each new head's changes before its root is set (the session's addresses follow a move, as the spec's
   * refs do: `AddressModel.followMoves`): a commit's own; on an undo, a redo or a rollback the moves between the two
   * commits (only `moved` is set). Returns the detach.
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
    this._lineage.set(commit.revision, { parent: before.revision, moved: changes.moved });
    batch(() => {
      this._head.set(commit);
      this._compose(extra);
      if (!this._editing) this._changed.set(commit);
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
    const frozen = this._frozen();
    if (frozen) return { ok: false, issues: [frozen] };
    const stale = patch.base !== undefined && patch.base !== head.revision && patch.result !== head.revision;
    if (stale && (options.rebase ?? this.config.edit.rebaseStale)) return this._rebase(patch, options);
    const r = applyPatch(head.root, patch, this._options(options));
    if (!r.ok) return { ok: false, issues: this._stepHints(patch, r.issues) };
    if (!r.changed) return { ok: true, changed: false, commit: head, changes: r.changes, issues: r.issues };
    const issues = [...r.issues, ...this._masked(r.changes)];
    const commit = this._commit(r.root, r.changes, options.context, issues);
    return { ok: true, changed: true, commit, changes: r.changes, issues };
  }

  /** While a full-state step is shown, an edit is refused: the step does not follow the spec, so it would show nothing. */
  private _frozen(): Issue | undefined {
    const step = this._step.value;
    if (step?.form !== 'full') return undefined;
    return issueOf(ViewerIssues, 'seq.viewer.step.frozen', `The step ${step.label} is a full state`, {
      address: step.address,
      hint: ViewerIssues.codes['seq.viewer.step.frozen'].hint,
    });
  }

  /**
   * Shows a frame over what is shown (_spec/README.md "SpecModel", The frame slot): an animation's transient spec, no
   * patch, no revision, no history step; undefined leaves it, what was composed shown again.
   */
  setFrame(root: Node | undefined): void {
    if (root === this._frame) return;
    this._frame = root;
    this._shown.set(root ?? this._composed);
  }

  /** The frame shown over the state, if any. */
  get frame(): Node | undefined {
    return this._frame;
  }

  /** A preview's options: no language checks at interactive rates, but the bounds of what it makes and sets. */
  private _previewOptions(): ApplyOptions {
    return { ...PreviewApply, validators: [seqBoundsValidator(this.config.bounds)] };
  }

  /** A story step's: the bounds too, so a stored step cannot build past them. */
  private _stepOptions(): ApplyOptions {
    return { ...StepApply, validators: [seqBoundsValidator(this.config.bounds)] };
  }

  private _options(options: SpecPatchOptions): ApplyOptions {
    return {
      ...SpecPatchApply,
      // the configuration's bounds on nesting, as the whole spec's check holds them
      validators: [seqPatchValidator(this.config.bounds)],
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
    // `merge3` infers no move: theirs' own, where they still stand (the session follows them as the refs did)
    const moved = theirs.changes.moved.filter(
      (x) => resolveAddress(head.root, x.from) !== undefined && resolveAddress(m.root, x.to) !== undefined,
    );
    const changes: Changes = { ...m.changes, moved };
    const issues = [...m.issues, ...this._masked(changes)];
    const commit = this._commit(m.root, changes, options.context, issues);
    return { ok: true, changed: true, commit, changes, issues, rebased: true };
  }

  /** A target refused while a step is shown, that only the step makes (it is in what is shown): the hint names it. */
  private _stepHints(patch: Patch, issues: readonly Issue[]): readonly Issue[] {
    const step = this._step.value;
    if (!step) return issues;
    const shown = this._shown.value;
    return issues.map((i) => {
      if (i.code !== 'repr.patch.target') return i;
      // the op's index is the issue's path in the patch tree (`['children', i, …]`)
      const index = i.path?.[0] === 'children' ? i.path[1] : undefined;
      const op = typeof index === 'number' ? patch.ops[index] : undefined;
      const at = op ? ('at' in op ? op.at : op.parent) : undefined;
      if (!at || !resolveAddress(shown, at)) return i;
      return { ...i, hint: `Only the story step ${step.address} makes it: leave the step to edit it` };
    });
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
    // a full-state step shown refuses every edit, a whole root's too (leave it first: a new document does)
    if (this._frozen()) return 'frozen-step';
    const next = typeof change === 'function' ? change(head.root) : change;
    this._preview.set(undefined);
    // the bounds first, before anything hashes or walks the root (a tree too deep for the stack): refused with the
    // issues, the spec as it was
    const bounded = next === head.root ? [] : checkSpecBounds(next, this.config.bounds);
    this._boundsIssues.set(bounded);
    if (bounded.length) return 'out-of-bounds';
    // before the commit, so what listens to the new root already sees a loaded spec
    this._loaded.set(true);
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
    if (this._history) this.detachHistory(this._history);
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
  edit(label: string, work: (edit: SpecEdit) => void, options: EditOptions = {}): EditOutcome | Promise<EditOutcome> {
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
      this._editing++;
      try {
        work(edit);
      } finally {
        this._editing--;
      }
      const outcome = {
        ok: results.every((r) => r.ok),
        changed: results.some((r) => r.ok && r.changed),
        results,
        issues: results.flatMap((r) => r.issues),
      };
      // the edit's commits are told once it ended well; a refused one is taken back, never told
      if (outcome.ok && outcome.changed && !this._editing) this._changed.set(this._head.value);
      return outcome;
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
      ...(options.within ? { within: options.within } : {}),
    });
  }

  /**
   * Applies patches as one edit (one step), in order; all or nothing. `patches` may be a function of the root and
   * revision the edit runs on: a command's patch is made there, when the history runs it, so edits issued one after
   * another without waiting each apply to the one before (made earlier, the second would be stale).
   */
  patch(
    label: string,
    patches: Patch | readonly Patch[] | PatchMaker,
    options: EditOptions & Omit<SpecPatchOptions, 'context'> = {},
  ): EditOutcome | Promise<EditOutcome> {
    const { origin, within, ...o } = options;
    return this.edit(
      label,
      (e) => {
        const made = typeof patches === 'function' ? patches(this.root.value, this.revision.value) : patches;
        const list: readonly Patch[] = Array.isArray(made) ? made : [made as Patch];
        for (const p of list) if (!e.apply(p, o).ok) break;
      },
      { ...(origin ? { origin } : {}), ...(within ? { within } : {}) },
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
    const frozen = this._frozen();
    if (frozen) return { ok: false, issues: [frozen] };
    const p: Patch = { ops: patch.ops };
    const root = this._head.value.root;
    const r = applyPatch(root, p, this._previewOptions());
    if (!r.ok) return r;
    // composed from this result: one apply a move
    this._previewApplied = { root, patch: p, result: r.root };
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
    // one batch: the preview gone and the new root shown are told together
    let r: SpecPatchResult | undefined;
    batch(() => {
      this._preview.set(undefined);
      r = this.applyPatch(p, options);
      if (!r.ok || !r.changed) this._compose();
    });
    return r;
  }

  // ---- story steps ----

  /**
   * Reads story steps with the core's (the viewer's story part, `StoryPart`): before, `enterStep` refuses every step
   * (`seq.viewer.step.none`), so a viewer that lists no story bundles none of it.
   */
  attachStory(story: SeqStorySteps): void {
    if (this._story === story) return;
    this._story = story;
    this._stories = story.createStoryCache();
  }

  /** The stored patch of a step, by id; undefined when the spec has none such. */
  private _stepPatch(id: string): Patch | undefined {
    const node = resolveAddress(this._head.value.root, ['snapshots', id]);
    return node?.kind === 'snapshot' ? this._story?.stepPatchOf(node) : undefined;
  }

  /**
   * Shows a story step over the root and the preview (navigation: no commit, no history step). Undefined when it was
   * entered (whether it applies is `step.applies`), else the issue why not.
   */
  enterStep(id: string): Issue | undefined {
    if (!this._story)
      return issueOf(ViewerIssues, 'seq.viewer.step.none', `This viewer enters no story step (${id})`, {
        hint: ViewerIssues.codes['seq.viewer.step.none'].hint,
      });
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
      const cached = this._previewApplied;
      const r =
        cached?.root === root && cached.patch === preview
          ? ({ ok: true, root: cached.result } as const)
          : applyPatch(root, preview, this._previewOptions());
      if (r.ok) {
        edited = r.root;
        this._previewApplied = { root, patch: preview, result: r.root };
      } else {
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
      if (node?.kind !== 'snapshot' || !this._story) this._stepId.set(undefined);
      else {
        // the core's state of the step over the edited root: a full state, the patch over it, or a delta (the ops that
        // still apply, flagged)
        const st = this._story.stateOf(edited, id, { cache: this._stories!, apply: this._stepOptions() })!;
        shown = st.root;
        for (const i of st.issues) if (StepIssueCodes.has(i.code)) issues.push(i);
        if (!st.applies)
          issues.push(
            issueOf(ViewerIssues, 'seq.viewer.step.stale', `The step ${st.label} no longer applies`, {
              address: st.address,
              hint: ViewerIssues.codes['seq.viewer.step.stale'].hint,
            }),
          );
        step = { id, address: st.address, label: st.label, applies: st.applies, form: st.form, state: st };
      }
    }
    const prev = this._step.value;
    if (
      prev?.id !== step?.id ||
      prev?.applies !== step?.applies ||
      prev?.label !== step?.label ||
      prev?.address !== step?.address ||
      prev?.state !== step?.state
    )
      this._step.set(step);
    this._edited.set(edited);
    this._composed = shown;
    this._shown.set(this._frame ?? shown);
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
    const moved = this._movesBetween(this._head.value, commit);
    for (const f of this._beforeCommit) f({ ...NoChanges, moved });
    this._remember(commit);
    batch(() => {
      this._head.set(commit);
      this._compose();
    });
  }
}
