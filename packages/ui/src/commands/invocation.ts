/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// One open command (commands §3): its values (`InvocationModel`), and typing into its fields, taking choices, the keys
// the menu forwards, the run, the assistant. DOM-free: the menu draws its states and calls its methods.

import { batch } from '@molstar/core/reactive';
import type { InvocationError, InvocationInit } from './invocation-state.js';
import { checkFormat } from './format.js';
import { errorMessage, InvocationModel } from './invocation-model.js';
import { isNumeric, isReadable, nextEmpty, readNumber } from './sentence.js';
import type { CommandFormat, CommandResult } from './types.js';

export class CommandInvocation<Ctx> extends InvocationModel<Ctx> {
  /** Stops a running assist. */
  readonly #abort = new AbortController();
  /** The run's signal: aborted by an explicit cancel, or a close. */
  #run: AbortController | undefined;

  constructor(init: InvocationInit<Ctx>) {
    super(init);
    // What remains of the tail is typed into the first empty field (and taken when its format reads it).
    if (this.tailRest.length) {
      this.type(this.tailRest.join(' '));
      const s = this.state.value;
      const f = this.fieldOf(s.field);
      if (f?.format && checkOk(this.formatOf(f), s.text, init.context.value)) this.accept();
    }
    this.updatePreview();
  }

  // --- Typing and choosing (commands §3.3) ---------------------------------------------------------

  /** Text typed into the active field: its choices follow; nothing is set until taken. */
  type(text: string): void {
    const s = this.state.value;
    const f = this.fieldOf(s.field);
    if (!f || !this.editing()) return;
    // One change: the text, then what its choices make of it (the derived `choices`, listed once for this text).
    batch(() => {
      this.patch({
        text,
        error: undefined,
        escalated: undefined,
        ...(s.proposed.includes(f.id) ? { proposed: s.proposed.filter((x) => x !== f.id) } : {}),
      });
      const items = this.choices.value?.items.length ?? 0;
      const typed = !!text.trim();
      this.patch({
        // A typed number is taken as typed: a preset is highlighted only by the arrows.
        pick: typed && items && !isNumeric(f) ? 0 : -1,
        canAssist:
          typed && !items && this.#assistable() && !isReadable(f, text, this.formatOf(f), this.init.context.value),
      });
    });
  }

  /** Moves the highlight through the choices, wrapping. */
  move(delta: number): void {
    const n = this.choices.value?.items.length ?? 0;
    if (!n || !this.editing()) return;
    const { pick } = this.state.value;
    this.patch({ pick: pick < 0 ? (delta > 0 ? 0 : n - 1) : (((pick + delta) % n) + n) % n });
  }

  /** Highlights a choice (under the pointer); -1 for none. */
  point(index: number): void {
    const n = this.choices.value?.items.length ?? 0;
    if (index >= -1 && index < n && index !== this.state.value.pick) this.patch({ pick: index });
  }

  /** Takes the highlighted choice (or `pick`), else the typed text when it reads as the field's value. */
  accept(pick?: number): boolean {
    const s = this.state.value;
    const f = this.fieldOf(s.field);
    if (!f || !this.editing()) return false;
    const chosen = this.choices.value?.items[pick ?? s.pick];
    if (chosen?.fill) {
      // An example to edit: into the text, not the value.
      this.type(chosen.value);
      return false;
    }
    if (chosen) return this.set(f.id, isNumeric(f) ? Number(chosen.value) : chosen.value);
    const text = s.text.trim();
    if (!text) return false;
    if (isNumeric(f)) {
      const v = readNumber(text);
      return v === undefined ? this.fail(f.id, `"${text}" is not a number.`) : this.set(f.id, v);
    }
    if (f.options) {
      const w = text.toLowerCase();
      const o = f.options.find((x) => x.value.toLowerCase() === w || x.label.toLowerCase() === w);
      if (o) return this.set(f.id, o.value);
      return this.fail(
        f.id,
        `"${text}" is not one of the choices.`,
        `Choose: ${f.options.map((x) => x.label).join(', ')}.`,
      );
    }
    const format = this.formatOf(f);
    if (!format) return this.set(f.id, text);
    const r = checkFormat(format, text, this.init.context.value);
    if (r.ok) return this.set(f.id, r.value);
    this.patch({
      error: { field: f.id, message: r.message, ...(r.hint ? { hint: r.hint } : {}) },
      canAssist: !!r.unreadable && this.#assistable(),
    });
    return false;
  }

  // --- Keys the menu forwards ----------------------------------------------------------------------

  /** Enter: take what is typed or highlighted (running too when that leaves no field empty), else run. */
  async enter(): Promise<CommandResult | undefined> {
    const s = this.state.value;
    if (s.status !== 'editing') return undefined;
    if (s.text.trim() || s.pick >= 0) {
      if (!this.accept()) return undefined;
      const values = this.state.value.values;
      if (!this.fields.every((f) => values[f.id] !== undefined)) return undefined;
    }
    return this.confirm();
  }

  /** Tab: let the assistant read unreadable text, else take a choice, else the next field. */
  async tab(delta: 1 | -1 = 1): Promise<void> {
    const s = this.state.value;
    if (s.canAssist) return this.assist();
    if ((s.text.trim() && this.choices.value?.items.length) || s.pick >= 0) {
      this.accept();
      return;
    }
    this.step(delta);
  }

  // --- Running (commands §3.5) ---------------------------------------------------------------------

  /** Runs the command with the current values: one outcome per started run (done, failed and still open, cancelled). */
  async confirm(): Promise<CommandResult | undefined> {
    const s = this.state.value;
    if (s.status !== 'editing') return undefined;
    if (s.blocked) {
      this.patch({ error: { message: s.blocked } });
      return undefined;
    }
    if (s.missing.length) {
      this.patch({ field: s.missing[0], text: '', pick: -1 });
      return undefined;
    }
    this.patch({ status: 'running', error: undefined });
    const run = new AbortController();
    this.#run = run;
    const ctx = this.init.context.value;
    const problem = this.#prepare(ctx);
    if (problem) {
      this.#run = undefined;
      this.patch({ status: 'editing', error: problem });
      return undefined;
    }
    let result: CommandResult | undefined;
    let failure: unknown;
    let failed = false;
    try {
      this.init.onRun?.();
      result = (await this.command.run(this.args(), ctx, run.signal)) ?? {};
    } catch (e) {
      failed = true;
      failure = e;
    }
    this.#run = undefined;
    // Closed by its owner (its registration went): it reports nothing.
    if (this.state.value.status === 'closed') return undefined;
    if (!failed) {
      // Done even when cancelled: the effect landed, so the person sees it (and its Undo).
      const done = result as CommandResult;
      this.clearPreview(true);
      this.patch({ status: 'done', result: done });
      this.init.onSettled?.({ kind: 'done', command: this.command, result: done });
      return done;
    }
    if (run.signal.aborted) {
      this.init.onSettled?.({ kind: 'cancelled', command: this.command });
      return undefined;
    }
    const text = errorMessage(failure);
    this.patch({ status: 'editing', error: { message: text } });
    this.init.onSettled?.({ kind: 'failed', command: this.command, message: text, open: true });
    return undefined;
  }

  /** Cancels the run in flight (an explicit act): its signal aborts; the outcome follows when it settles. */
  cancel(): boolean {
    if (this.state.value.status !== 'running' || !this.#run) return false;
    this.#run.abort();
    this.patch({ status: 'cancelling' });
    return true;
  }

  /** Lets the assistant read the typed words into values (proposed; nothing runs). */
  async assist(): Promise<void> {
    const s = this.state.value;
    const assistant = this.init.assistant;
    if (!assistant || !this.command.assist || s.status !== 'editing') return;
    this.patch({ status: 'resolving', canAssist: false, error: undefined, note: undefined, escalated: undefined });
    const signal = this.#abort.signal;
    const ctx = this.init.context.value;
    try {
      const r = await assistant.propose(
        {
          command: this.command,
          params: this.params,
          values: this.values(),
          text: s.text,
          ...(s.field ? { field: s.field } : {}),
          ctx,
        },
        signal,
      );
      if (signal.aborted || this.state.value.status === 'closed') return;
      if (r.kind === 'answer') {
        const result: CommandResult = { detail: r.text };
        this.patch({ status: 'done', result });
        this.init.onSettled?.({ kind: 'done', command: this.command, result });
        return;
      }
      if (r.kind === 'escalate') {
        this.patch({ status: 'editing', escalated: s.text, error: { message: r.reason } });
        return;
      }
      this.patch({ status: 'editing', text: '', pick: -1 });
      const rejected: string[] = [];
      for (const [k, v] of Object.entries(r.values)) {
        if (v === undefined) continue;
        const f = this.fieldOf(k);
        const checked = f && f.format && typeof v === 'string' ? checkWith(this.formatOf(f), v, ctx) : undefined;
        if (checked && !checked.ok) rejected.push(`${k}: ${checked.message}`);
        else if (!this.set(k, checked?.ok ? checked.value : v, { proposed: true })) rejected.push(k);
      }
      this.patch({
        ...(r.note ? { note: r.note } : {}),
        field: nextEmpty(this.fields, this.state.value.values) ?? this.state.value.field,
        ...(rejected.length ? { error: { message: `The proposal did not fit: ${rejected.join('; ')}.` } } : {}),
      });
    } catch (e) {
      if (signal.aborted) return;
      this.patch({ status: 'editing', error: { message: errorMessage(e) } });
    }
  }

  /**
   * Closes it at once: drops the preview, stops a running assist. A run in flight is aborted and reports nothing (its
   * owner went); hosts go through the manager, which refuses to close a running invocation and cancels explicitly.
   */
  close(): void {
    if (this.state.value.status === 'closed') return;
    this.#abort.abort();
    this.#run?.abort();
    this.patch({ status: 'closed' });
    this.clearPreview();
    this.retire();
  }

  #assistable(): boolean {
    return !!this.command.assist && !!this.init.assistant;
  }

  /** Why the run cannot start: `when`, a value its format no longer reads or its field no longer takes, `check`. */
  #prepare(ctx: Ctx): InvocationError | undefined {
    const w = this.blockedReason();
    if (w) return { message: w };
    const values = this.state.value.values;
    for (const f of this.fields) {
      const v = values[f.id];
      const format = this.formatOf(f);
      if (format && typeof v === 'string') {
        const r = checkFormat(format, v, ctx);
        if (!r.ok) return { field: f.id, message: r.message, ...(r.hint ? { hint: r.hint } : {}) };
      }
    }
    try {
      const problem = this.command.check?.(this.args(), ctx);
      return problem ? { message: problem } : undefined;
    } catch (e) {
      return { message: errorMessage(e) };
    }
  }
}

function checkWith<Ctx>(format: CommandFormat<Ctx> | undefined, text: string, ctx: Ctx) {
  return format ? checkFormat(format, text, ctx) : undefined;
}

function checkOk<Ctx>(format: CommandFormat<Ctx> | undefined, text: string, ctx: Ctx): boolean {
  return !!checkWith(format, text, ctx)?.ok;
}
