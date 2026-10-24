/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/core';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import type { ViewerApi } from '@molstar/prism-viewer/api/area';
import type { PrismStore } from '@molstar/prism-viewer/spec';
import { formatDiagnostic } from '../language/diagnostics.js';
import { type Dialect, DialectNames } from '../language/ir.js';
import { dialectOfFile, dialectPrefix } from '../language/script.js';
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
}

export interface ConsoleEntry {
  readonly id: number;
  readonly dialect: Dialect;
  /** What ran: the text, or the file's name. */
  readonly input: string;
  readonly state: 'running' | 'done' | 'failed';
  readonly lines: readonly ConsoleLine[];
}

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
    return this.#run(dialect, prefixed ? prefixed.rest : text, text, files);
  }

  /** Opens script files with the files they load: each script runs, in the dialect its extension says. */
  async openFiles(files: readonly File[]): Promise<void> {
    const scripts = files.filter((f) => dialectOfFile(f.name));
    if (!scripts.length) {
      this.#push({
        dialect: this.#dialect.value,
        input: files.map((f) => f.name).join(', '),
        state: 'failed',
        lines: [
          { kind: 'error', text: 'No script among these: open a .pml, .tcl or .vmd file (with the files it loads).' },
        ],
      });
      return;
    }
    for (const script of scripts) {
      const text = await script.text();
      await this.#run(dialectOfFile(script.name)!, text, script.name, files);
    }
  }

  async #run(dialect: Dialect, text: string, input: string, files: readonly File[]): Promise<ConsoleEntry> {
    const area = this.api.script;
    const id = this.#push({ dialect, input, state: 'running', lines: [] });
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
        errors.map((d) => ({ kind: 'error', text: formatDiagnostic(text, d) })),
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
      return this.#finish(id, 'failed', [{ kind: 'error', text: r.hint ? `${r.message} ${r.hint}` : r.message }]);
    }
    const value = r.value as ScriptRunValue | undefined;
    const lines: ConsoleLine[] = (value?.outputs ?? []).map((o) => ({ kind: o.kind, text: o.text }));
    return this.#finish(id, 'done', lines.length ? lines : [{ kind: 'did', text: r.did }]);
  }

  #push(e: Omit<ConsoleEntry, 'id'>): number {
    const id = this.#next++;
    this.#entries.set([...this.#entries.value, { ...e, id }]);
    return id;
  }

  #finish(id: number, state: ConsoleEntry['state'], lines: readonly ConsoleLine[]): ConsoleEntry {
    let done: ConsoleEntry | undefined;
    this.#entries.set(
      this.#entries.value.map((e) => {
        if (e.id !== id) return e;
        done = { ...e, state, lines };
        return done;
      }),
    );
    return done!;
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
    return d === 'pymol'
      ? `${DialectNames.pymol}: fetch 1cbs; show sticks, resn REA · help`
      : `${DialectNames.vmd}: mol pdbload 1cbs; mol addrep top · help`;
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
