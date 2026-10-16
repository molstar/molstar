/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Message } from '../../messages/index.js';
import { HarnessError, type HeadEdit } from '../contracts/harness.js';
import { type HeadId, MAIN, type ScopeId } from '../contracts/records.js';
import type { ToolJournal } from '../contracts/runner.js';
import { append, applyEdit, openPointer, type Run, readPointer, type Writer } from '../loop/run.js';
import type { Staged } from '../loop/state.js';

/**
 * A tool's own context for one attempt (harness.md §8): writes are staged in the step, the first takes a
 * hold on the pointer, and a failed or aborted attempt's writes and holds are dropped.
 */
export class Journal implements ToolJournal, Writer {
  private readonly run: Run;
  private readonly holder: string;
  private readonly saved = new Map<Staged, { entries: string[]; appended: string[] | undefined }>();
  private chain: Promise<void> = Promise.resolve();
  private failure: unknown;

  constructor(run: Run, holder: string) {
    this.run = run;
    this.holder = holder;
  }

  async read(scope: ScopeId): Promise<{ head?: HeadId; messages: Message[] }> {
    await this.chain;
    return readPointer(this.run, scope);
  }

  append(scope: ScopeId, messages: readonly Message[]): void {
    this.hold(scope);
    // Queued behind an edit (which sizes the context, so it is async): writes apply in order.
    this.queue(async () => append(this.run, this.open(scope), messages));
  }

  edit(scope: ScopeId, edit: HeadEdit): void {
    this.hold(scope);
    this.queue(async () => {
      if (!(await applyEdit(this.run, this.open(scope), edit)))
        throw new HarnessError('invalid', `the edit of ${scope} leaves an unprojectable head`);
    });
  }

  /** Holds a pointer for this writer, or fails `pointer-busy` (harness.md §8). */
  hold(scope: ScopeId): void {
    const holders = this.run.collector.holders;
    const holder = holders.get(scope);
    if (scope === MAIN || (holder !== undefined && holder !== this.holder))
      throw new HarnessError('pointer-busy', `${scope} is written by another writer in this step`);
    holders.set(scope, this.holder);
  }

  /** The pointer's transient head, held, with its state kept so a drop can restore it. */
  open(scope: ScopeId): Staged {
    this.hold(scope);
    const staged = openPointer(this.run, scope);
    if (!this.saved.has(staged))
      this.saved.set(staged, { entries: [...staged.entries], appended: staged.appended && [...staged.appended] });
    return staged;
  }

  /** Waits for queued writes; a write that failed fails the attempt. */
  async settled(): Promise<void> {
    await this.chain;
    if (this.failure !== undefined) throw this.failure;
  }

  /** The attempt failed or was aborted: its writes and holds go. */
  drop(): void {
    for (const [staged, saved] of this.saved) {
      staged.entries = saved.entries;
      staged.appended = saved.appended;
      if (staged.scope === undefined) continue;
      this.run.collector.compaction.delete(staged.scope);
      if (this.run.collector.holders.get(staged.scope) === this.holder) this.run.collector.holders.delete(staged.scope);
    }
    for (const [scope, holder] of this.run.collector.holders)
      if (holder === this.holder) this.run.collector.holders.delete(scope);
    this.saved.clear();
  }

  private queue(write: () => Promise<void>): void {
    this.chain = this.chain.then(write).catch((error: unknown) => {
      this.failure ??= error;
    });
  }
}
