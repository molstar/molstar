/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `@molstar/ui/history` adapter of Seq*'s history port (_spec/README.md "Ports", History), for the composition
// roots that use that stack: the standalone viewer (`uiHistory`, a stack of its own) and a host with one stack for
// everything (`hostHistory`: Seq*'s spec is one participant of `HistoryParticipants`, beside the host's own). An edit
// runs as one `HistoryModel.run` step, an agent's turn is a `beginGroup`, a change outside a run invalidates the earlier
// steps, and the commits the stack checkpointed resolve by revision (`root`).

import type { Node } from '@molstar/core/repr/tree';
import { HistoryModel, type HistoryOptions, type HistoryPort } from '@molstar/ui/history/history-model';
import type { SeqHistory, SeqHistoryParticipant, SeqHistoryPort, SpecCommit } from './history.js';

/** The roots the stack's checkpoints hold, by revision, the latest `keep` of them. */
class CheckpointedRoots {
  private readonly roots = new Map<string, Node>();
  constructor(private readonly keep: () => number) {}
  note(commit: SpecCommit): void {
    this.roots.delete(commit.revision);
    this.roots.set(commit.revision, commit.root);
    while (this.roots.size > Math.max(1, this.keep())) this.roots.delete(this.roots.keys().next().value!);
  }
  get(revision: string): Node | undefined {
    return this.roots.get(revision);
  }
}

/** The participant as a `HistoryPort`, noting each commit it checkpoints. */
function recording(participant: SeqHistoryParticipant, roots: CheckpointedRoots): HistoryPort<SpecCommit> {
  return {
    checkpoint: () => {
      const c = participant.checkpoint();
      roots.note(c);
      return c;
    },
    restore: (c) => participant.restore(c),
  };
}

/** How many commits a stack's checkpoints may hold: its steps and its redo steps. */
const heldBy = (stack: HistoryModel<any>) => () => stack.limits.capacity + stack.limits.maxRedo + 1;

function historyOver(
  stack: HistoryModel<any>,
  roots: CheckpointedRoots,
  dispose: () => void,
): Omit<SeqHistory, 'disposeWithViewer'> {
  return {
    run: (step, work, options) =>
      stack
        .run(step.label, step.origin, async (context) => work(context), {
          changed: options.changed,
          failed: options.failed,
          ...(options.within ? { within: options.within } : {}),
        })
        .then((o) => o.value),
    beginTurn: (label, origin) => stack.beginGroup(label, origin),
    isActive: (context) => stack.isActive(context),
    // `@molstar/ui/history`'s rule: a change outside the stack makes its earlier steps unsafe
    outside: (reason) => stack.invalidate(reason),
    root: (revision) => roots.get(revision),
    dispose,
  };
}

/** Seq*'s history over a stack of its own: the standalone viewer's. `stack` is its Undo and Redo. */
export interface UiHistory extends SeqHistory {
  readonly stack: HistoryModel<SpecCommit>;
}

/** The standalone viewer's history: a `@molstar/ui/history` stack of the spec alone, disposed with the viewer. */
export function uiHistory(options: HistoryOptions = {}): SeqHistoryPort {
  return (participant): UiHistory => {
    let stack: HistoryModel<SpecCommit> | undefined;
    const roots = new CheckpointedRoots(() => (stack ? heldBy(stack)() : 1));
    stack = new HistoryModel<SpecCommit>(recording(participant, roots), options);
    const s = stack;
    return { ...historyOver(s, roots, () => s.dispose()), stack: s };
  };
}

/**
 * A host's stack of several participants under one `HistoryPort` (Seq*'s spec, a structure's state): a checkpoint is
 * each one's, by key; a restore restores each. Give it to the host's `HistoryModel`, and the participants join as
 * they come (`add`).
 */
export class HistoryParticipants implements HistoryPort<Readonly<Record<string, unknown>>> {
  private readonly parts = new Map<string, HistoryPort<unknown>>();

  /** Adds a participant under a key; returns its removal. */
  add<C>(key: string, port: HistoryPort<C>): () => void {
    if (this.parts.has(key)) throw new Error(`A participant ${key} is already in the history`);
    this.parts.set(key, port as HistoryPort<unknown>);
    return () => {
      if (this.parts.get(key) === port) this.parts.delete(key);
    };
  }

  checkpoint(): Readonly<Record<string, unknown>> {
    const out: Record<string, unknown> = {};
    for (const [key, port] of this.parts) out[key] = port.checkpoint();
    return out;
  }

  async restore(checkpoint: Readonly<Record<string, unknown>>): Promise<void> {
    for (const [key, port] of this.parts) if (Object.hasOwn(checkpoint, key)) await port.restore(checkpoint[key]);
  }
}

/**
 * Seq*'s history on a host's stack (`stack`, over `participants`): the spec joins as the participant `key` while the
 * viewer lives; the stack is the host's and stays when the viewer goes.
 */
export function hostHistory(stack: HistoryModel<any>, participants: HistoryParticipants, key = 'seq'): SeqHistoryPort {
  return (participant) => {
    const roots = new CheckpointedRoots(heldBy(stack));
    const remove = participants.add(key, recording(participant, roots));
    return historyOver(stack, roots, remove);
  };
}
