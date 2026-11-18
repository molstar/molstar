/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Seq*'s history port (_spec/README.md "Ports", History): Seq* owns no undo stack; the host does (a structure viewer,
// Prism, a Mol* plugin keeps one stack for everything). The port has two faces, in Seq*'s own terms:
//
// - the participant (the host drives Seq*): `checkpoint()` is a commit, the spec's immutable root with its revision
//   and parents (cheap by structural sharing); `restore(c)` moves the root back (no new commit) and drops a preview;
//   `changed` tells of each commit;
// - the history (Seq* asks the host): an edit runs as one host step with its label, change kind and origin, so the
//   host's stack interleaves Seq*'s edits with its own; a change made outside a run is reported (`outside`); `root`
//   resolves a revision older than the spec's recent roots, for a rebase.
//
// Only spec patches enter it. Adapters live in the composition roots: `@molstar/ui/history` (`uiHistory`,
// `hostHistory`) for the standalone viewer and Prism; none headless.

import type { ReadonlyState } from '@molstar/core/reactive';
import type { Node } from '@molstar/core/repr/tree';

/** A revision as a history records it: the spec's root, its revision (the root's hash) and its parents' revisions. */
export interface SpecCommit {
  readonly root: Node;
  readonly revision: string;
  /** One (an edit), two (a merge), none (a document loaded). */
  readonly parents: readonly string[];
}

/** Who made an edit, as a host's history groups steps: `person`, `agent` (with its turn), `script`, `link`. */
export interface EditOrigin {
  readonly kind: string;
  /** An agent's turn: its edits join one step while it is the latest. */
  readonly turn?: string;
}

/** What a step is, for the host's list and its Undo button. */
export interface EditStep {
  readonly label: string;
  /** Always `patch`: only spec patches enter a history (_spec/README.md "Commands"). */
  readonly kind: 'patch';
  readonly origin: EditOrigin;
}

/** An edit in progress: what work passes back so its changes are part of the step. Opaque. */
export interface EditContext {
  readonly id: number;
}

/** How a run's value says what happened: no step when nothing changed, the checkpoint restored when it failed. */
export interface EditRunOptions<T> {
  readonly changed: (value: T) => boolean;
  readonly failed: (value: T) => boolean;
  /** A run of the host's in progress this edit is part of (its work runs a Seq* command): no step of its own. */
  readonly within?: EditContext;
}

/** The participant face: what the spec gives the host (`SpecModel` is one). */
export interface SeqHistoryParticipant {
  /** The current commit: synchronous and cheap. */
  checkpoint(): SpecCommit;
  /** Moves the root back to a commit (no new commit) and drops a preview. */
  restore(commit: SpecCommit): Promise<void>;
  /** The last commit an edit made (a restore makes none). */
  readonly changed: ReadonlyState<SpecCommit>;
}

/** The history face: the host's stack as Seq* asks it. */
export interface SeqHistory {
  /** Runs an edit as one step of the host's stack: a checkpoint before it, restored when it fails. */
  run<T>(step: EditStep, work: (context: EditContext) => T, options: EditRunOptions<T>): Promise<T>;
  /** Opens an agent's turn (`origin.turn`): its edits join one step while it is the latest; returns its end. */
  beginTurn?(label: string, origin: EditOrigin): () => void;
  /** Whether a context is an edit still in progress. */
  isActive(context: EditContext | undefined): boolean;
  /** The spec changed outside a run: whether that invalidates the earlier steps is the host's rule. */
  outside(reason: string): void;
  /** A root by its revision, when the host's stack holds it (older than the spec's recent roots). */
  root?(revision: string): Node | undefined;
  /** False: the host owns it (a stack shared with other participants), and the viewer leaves it. */
  readonly disposeWithViewer?: boolean;
  dispose?(): void;
}

/** The port in a viewer's spec: a factory given the participant, called once per viewer. */
export type SeqHistoryPort = (participant: SeqHistoryParticipant) => SeqHistory;
