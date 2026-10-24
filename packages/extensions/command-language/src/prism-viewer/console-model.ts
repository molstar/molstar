/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/core';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerApi } from '@molstar/prism-viewer/api/area';
import type { PrismStore } from '@molstar/prism-viewer/spec';
import { applyCompletion, type CompletionItem, complete, type Signature } from '../language/completion.js';
import { describeCommand } from '../language/describe.js';
import { formatDiagnostic } from '../language/diagnostics.js';
import type { Dialect } from '../language/ir.js';
import { dialectOfFile, dialectPrefix, rootScripts } from '../language/script.js';
import { defaultDialect } from './commands.js';
import { type ScriptRunValue, scriptTexts } from './operations.js';

// The console's model (_spec §9), DOM-free: a dialect, a draft, the submissions with what each printed and said, and
// the person's history of lines (kept with the preferences, the last `HistoryMax`). A submission runs as `script.run`
// (one step) or, `undo`/`redo` alone, as the viewer's history. A script file opened with the files it loads runs whole.

export const HistoryPreference = 'command-language.history';
const HistoryMax = 200;

export type ConsoleLineKind = 'did' | 'value' | 'warning' | 'note' | 'error';

export interface ConsoleLine {
  readonly kind: ConsoleLineKind;
  readonly text: string;
  /** What to do instead (an error's). */
  readonly hint?: string;
}

export interface ConsoleEntry {
  readonly id: number;
  readonly dialect: Dialect;
  /** What ran: the text, or the file's name. */
  readonly input: string;
  readonly state: 'running' | 'done' | 'failed';
  readonly lines: readonly ConsoleLine[];
  /** How long it ran (ms), once done. */
  readonly ms?: number;
  /** It ran files (a script, a drop), not a typed line. */
  readonly files?: true;
}

/** What the input will do, read as it is typed (the composer's line under it). */
export type ConsoleInspection =
  | { readonly kind: 'empty' }
  | { readonly kind: 'ready'; readonly dialect: Dialect; readonly does: string; readonly warnings: readonly string[] }
  | { readonly kind: 'error'; readonly dialect: Dialect; readonly message: string; readonly hint?: string };

/** A completion of the input: its word, and the input with it. */
export interface ConsoleCompletion {
  readonly label: string;
  readonly detail?: string;
  /** What it is: a command, a representation, a color, a keyword, an object, … */
  readonly kind: CompletionItem['kind'];
  /** A color's RGB, for its swatch. */
  readonly rgb?: number;
  /** The input with it. */
  readonly text: string;
}

/** What helps the person type the input: the completions of the word at its end, and the command it is in. */
export interface ConsoleAssist {
  readonly items: readonly ConsoleCompletion[];
  /** The part of the word typed so far (the items' labels start with it, or hold it). */
  readonly typed: string;
  readonly signature?: Signature;
}

/** Lines that show what each dialect reads (the empty console's examples). */
export const ConsoleExamples: Readonly<Record<Dialect, readonly string[]>> = {
  pymol: [
    'fetch 1cbs',
    'as cartoon',
    'show sticks, resn REA',
    'util.cbag resn REA',
    'color slate, chain A',
    'spectrum b, blue_white_red',
    'zoom resn REA, 5',
    'help',
  ],
  vmd: [
    'mol pdbload 1cbs',
    'mol delrep 0 top',
    'mol representation NewCartoon',
    'mol color Structure',
    'mol addrep top',
    'set sel [atomselect top "resname REA"]',
    '$sel num',
    'help',
  ],
};

/** The script files a console runs, by extension; the other files are what they load. */
export const ScriptExtensions = ['pml', 'tcl', 'vmd'];

export class ConsoleModel {
  readonly #dialect = observableState<Dialect>('pymol');
  readonly #entries = observableState<readonly ConsoleEntry[]>([]);
  readonly #draft = observableState('');
  readonly #history: string[];
  /** Where ↑/↓ are in the history; `history.length`: the draft. */
  #at: number;
  #saved = '';
  #next = 0;
  /** The runs, one after another: a submission made while one runs reads the scene that one leaves. */
  #queue: Promise<unknown> = Promise.resolve();

  constructor(
    readonly api: ViewerApi,
    readonly preferences?: PrismStore,
  ) {
    this.#history = readHistory(preferences);
    this.#at = this.#history.length;
    const viewer = api.script;
    const preferred = viewer ? defaultDialect(preferences, viewer) : undefined;
    if (preferred) this.#dialect.set(preferred);
  }

  get dialect(): ReadonlyState<Dialect> {
    return this.#dialect.asReadonly;
  }
  get entries(): ReadonlyState<readonly ConsoleEntry[]> {
    return this.#entries.asReadonly;
  }
  get draft(): ReadonlyState<string> {
    return this.#draft.asReadonly;
  }
  get history(): readonly string[] {
    return this.#history;
  }

  setDialect(d: Dialect): void {
    this.#dialect.set(d);
  }

  setDraft(text: string): void {
    this.#draft.set(text);
    this.#at = this.#history.length;
  }

  /** The line before in the history (↑); the draft is kept to come back to. */
  previous(): void {
    if (this.#at === 0) return;
    if (this.#at === this.#history.length) this.#saved = this.#draft.value;
    this.#at--;
    this.#draft.set(this.#history[this.#at]);
  }

  /** The line after (↓), back to the draft at the end. */
  next(): void {
    if (this.#at >= this.#history.length) return;
    this.#at++;
    this.#draft.set(this.#at === this.#history.length ? this.#saved : this.#history[this.#at]);
  }

  clear(): void {
    this.#entries.set([]);
  }

  /** Runs the draft (or the text given); a `pymol ` or `vmd ` first word runs it in that dialect. */
  async submit(text = this.#draft.value, files: readonly File[] = []): Promise<ConsoleEntry | undefined> {
    if (!text.trim()) return undefined;
    this.#remember(text);
    this.#draft.set('');
    this.#at = this.#history.length;
    const prefixed = dialectPrefix(text);
    const dialect = prefixed?.dialect ?? this.#dialect.value;
    return this.#enqueue(() => this.#run(dialect, prefixed ? prefixed.rest : text, text, files));
  }

  /** Runs after what was submitted before; the entry shows at once, as running. */
  #enqueue<T>(work: () => Promise<T>): Promise<T> {
    const next = this.#queue.then(work, work);
    this.#queue = next.catch(() => undefined);
    return next;
  }

  /**
   * Opens script files with the files they load: each script no other one includes runs (the others run where they
   * are included), in the dialect its extension says. Files without a script are loaded, as the line that loads them
   * (`load 1cbs.cif`, `mol new 1cbs.cif`): the person sees the command, and can run it again.
   */
  async openFiles(files: readonly File[]): Promise<void> {
    const texts = await scriptTexts(files);
    const roots = new Set(rootScripts(texts));
    // Scripts that all include each other have no root: each runs (and says it includes itself).
    const scripts = files.filter((f) => (roots.size ? roots.has(f.name) : texts.has(f.name)));
    if (!scripts.length) {
      const dialect = this.#dialect.value;
      const quoted = (name: string) => (/[\s,;"'{}]/.test(name) ? `"${name.replace(/"/g, '')}"` : name);
      const text = files
        .map((f) => (dialect === 'pymol' ? `load ${quoted(f.name)}` : `mol new ${quoted(f.name)}`))
        .join('\n');
      await this.#enqueue(() => this.#run(dialect, text, text, files));
      return;
    }
    for (const script of scripts) {
      await this.#enqueue(() =>
        this.#run(dialectOfFile(script.name)!, texts.get(script.name) ?? '', script.name, files, true),
      );
    }
  }

  async #run(
    dialect: Dialect,
    text: string,
    input: string,
    files: readonly File[],
    file = false,
  ): Promise<ConsoleEntry> {
    const area = this.api.script;
    const id = this.#push({ dialect, input, state: 'running', lines: [], ...(file ? { files: true as const } : {}) });
    if (!area) return this.#finish(id, 'failed', [{ kind: 'error', text: 'This viewer has no command language.' }]);
    // Everything that does not read, at once (the run would say only the first).
    const preview = area.preview(
      text,
      dialect,
      files.map((f) => f.name),
      await scriptTexts(files),
    );
    const errors = preview.diagnostics.filter((d) => d.severity === 'error');
    if (errors.length) {
      return this.#finish(
        id,
        'failed',
        errors.map((d) => ({
          kind: 'error',
          text: formatDiagnostic(text, { ...d, hint: undefined } as typeof d),
          ...(d.hint ? { hint: d.hint } : {}),
        })),
      );
    }
    const lone = preview.commands.length === 1 ? preview.commands[0] : undefined;
    if (lone?.kind === 'history') {
      const r = await (lone.action === 'undo' ? this.api.history.undo() : this.api.history.redo());
      return this.#finish(id, r.ok ? 'done' : 'failed', [
        { kind: r.ok ? 'did' : 'error', text: r.ok ? r.did : r.message },
      ]);
    }
    const r = await area.run({ dialect, text, files: [...files] });
    if (!r.ok) {
      return this.#finish(id, 'failed', [{ kind: 'error', text: r.message, ...(r.hint ? { hint: r.hint } : {}) }]);
    }
    const value = r.value as ScriptRunValue | undefined;
    const lines: ConsoleLine[] = (value?.outputs ?? []).map((o) => ({ kind: o.kind, text: o.text }));
    return this.#finish(id, 'done', lines.length ? lines : [{ kind: 'did', text: r.did }]);
  }

  #push(e: Omit<ConsoleEntry, 'id'>): number {
    const id = this.#next++;
    this.#entries.set([...this.#entries.value, { ...e, id }]);
    this.#started.set(id, performance.now());
    return id;
  }

  readonly #started = new Map<number, number>();

  #finish(id: number, state: ConsoleEntry['state'], lines: readonly ConsoleLine[]): ConsoleEntry {
    let done: ConsoleEntry | undefined;
    const started = this.#started.get(id);
    this.#started.delete(id);
    const ms = started === undefined ? undefined : Math.round(performance.now() - started);
    this.#entries.set(
      this.#entries.value.map((e) => {
        if (e.id !== id) return e;
        done = { ...e, state, lines, ...(ms !== undefined ? { ms } : {}) };
        return done;
      }),
    );
    return done!;
  }

  /** What a line will do, or why it will not (nothing runs). */
  inspect(text: string): ConsoleInspection {
    const area = this.api.script;
    if (!text.trim() || !area) return { kind: 'empty' };
    const prefixed = dialectPrefix(text);
    const dialect = prefixed?.dialect ?? this.#dialect.value;
    const body = prefixed ? prefixed.rest : text;
    if (!body.trim()) return { kind: 'empty' };
    const compiled = area.preview(body, dialect);
    const error = compiled.diagnostics.find((d) => d.severity === 'error');
    if (error) {
      return {
        kind: 'error',
        dialect,
        message: formatDiagnostic(body, { ...error, hint: undefined } as typeof error),
        ...(error.hint ? { hint: error.hint } : {}),
      };
    }
    const does = compiled.commands.map(describeCommand).filter(Boolean);
    return {
      kind: 'ready',
      dialect,
      does: does.length > 3 ? `${does.slice(0, 3).join('; ')}; and ${does.length - 3} more` : does.join('; '),
      warnings: compiled.diagnostics.filter((d) => d.severity !== 'error').map((d) => d.message),
    };
  }

  /** The words that may come at the input's end, each with the input completed by it. */
  completions(text: string, max = 6): ConsoleCompletion[] {
    return this.assist(text, max).items.slice();
  }

  /**
   * The completions of the word at the input's end (at most `max`), and the command the input is in. An empty input has
   * none unless `asked` (Tab on it: every command).
   */
  assist(text: string, max = 12, asked = false): ConsoleAssist {
    const area = this.api.script;
    if (!area || (!text.trim() && !asked)) return { items: [], typed: '' };
    const prefixed = dialectPrefix(text);
    const dialect = prefixed?.dialect ?? this.#dialect.value;
    const head = prefixed ? text.slice(0, prefixed.offset) : '';
    const body = prefixed ? prefixed.rest : text;
    const c = complete(body, dialect, area.env());
    const typed = body.slice(c.from, c.to);
    const items = c.items
      // What is typed whole already is no completion.
      .filter((i) => i.label.toLowerCase() !== typed.toLowerCase())
      .slice(0, max)
      .map(
        (i): ConsoleCompletion => ({
          label: i.label,
          kind: i.kind,
          ...(i.detail ? { detail: i.detail } : {}),
          ...(i.rgb !== undefined ? { rgb: i.rgb } : {}),
          text: head + applyCompletion(body, c, i),
        }),
      );
    return { items, typed, ...(c.signature ? { signature: c.signature } : {}) };
  }

  /** An entry's text back in the input, to change and run again. */
  edit(entry: ConsoleEntry): void {
    if (!this.editable(entry)) return;
    this.setDialect(entry.dialect);
    this.setDraft(entry.input);
  }

  /** Whether an entry's input is a line to edit (not a script file's name, nor the files of a drop). */
  editable(entry: ConsoleEntry): boolean {
    return !entry.files;
  }

  #remember(text: string): void {
    if (this.#history[this.#history.length - 1] !== text) this.#history.push(text);
    if (this.#history.length > HistoryMax) this.#history.splice(0, this.#history.length - HistoryMax);
    try {
      this.preferences?.set(HistoryPreference, this.#history);
    } catch {
      // A full or unavailable store: the history lasts for the session.
    }
  }

  /** What the input says it reads, for its placeholder. */
  placeholder(d: Dialect): string {
    return d === 'pymol' ? 'A PyMOL command: show sticks, resn REA' : 'A VMD command: mol addrep top';
  }
}

function readHistory(preferences: PrismStore | undefined): string[] {
  try {
    const v = preferences?.get(HistoryPreference);
    return Array.isArray(v) ? v.filter((x): x is string => typeof x === 'string').slice(-HistoryMax) : [];
  } catch {
    return [];
  }
}
