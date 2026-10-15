/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The command manager: what can run, what matches, what is open, what runs, and what came of it.

import { type Observable, Subject } from 'rxjs';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';
import { CommandInvocation } from './invocation.js';
import { fuzzyScore } from './match.js';
import type {
  CommandAssistant,
  CommandDefinition,
  CommandFormats,
  CommandGroup,
  CommandOutcome,
  CommandPreset,
  CommandResult,
  KeyChord,
  PreviewSink,
  RecentsStore,
} from './types.js';

export interface CommandManagerOptions<Ctx> {
  /** What commands see: the plugin, the selection, the focused view… */
  context: ReadonlyState<Ctx>;
  formats: CommandFormats<Ctx>;
  preview?: PreviewSink;
  assistant?: CommandAssistant<Ctx>;
  recents?: RecentsStore;
  maxRecents?: number;
  /**
   * What in the context previews read (the scene, the selection): an open form's preview is redrawn when one of them
   * changes (`Object.is`, item by item). Without it, every context emission redraws it.
   */
  previewDeps?: (ctx: Ctx) => readonly unknown[];
}

export interface CommandMatch<Ctx> {
  command: CommandDefinition<Ctx>;
  /** Why it cannot run now, if it cannot. */
  blocked?: string;
  /** The way in it was found by, if one: opened with `open(id, tail, preset.id)`. */
  preset?: CommandPreset;
}

export interface CommandQuery<Ctx> {
  /** Recently used first when nothing is typed; best match first otherwise. */
  recent: CommandMatch<Ctx>[];
  matches: CommandMatch<Ctx>[];
  /** The text starts with a command's alias and goes on: that command, and the rest as arguments. */
  direct?: CommandMatch<Ctx> & { tail: string };
}

export type OpenResult<Ctx> =
  | { kind: 'form'; invocation: CommandInvocation<Ctx> }
  | { kind: 'ran'; result: Promise<CommandResult | undefined> }
  | { kind: 'blocked'; reason: string };

export interface KeyEventLike {
  key: string;
  altKey?: boolean;
  shiftKey?: boolean;
  metaKey?: boolean;
  ctrlKey?: boolean;
  /** Layout-independent key, e.g. "KeyT": used when `key` is not a plain letter (⌥T on macOS is "†"). */
  code?: string;
  /** Where the key went (an element in a page); a text-entry target keeps keys that type. */
  target?: unknown;
  /** The composed path of a DOM event: across a shadow root `target` is its host, the path's first item the real one. */
  composedPath?(): EventTarget[];
  isComposing?: boolean;
  repeat?: boolean;
  defaultPrevented?: boolean;
}

/** The run in flight. */
export interface CommandRun<Ctx> {
  command: CommandDefinition<Ctx>;
  /** Cancel was asked for; the run has not settled yet. */
  cancelling: boolean;
}

interface Run<Ctx> {
  command: CommandDefinition<Ctx>;
  /** The registration it ran under: its outcome reports only while that lasts. */
  owner: object | undefined;
  invocation?: CommandInvocation<Ctx>;
  cancel(): void;
}

const message = (e: unknown) => (e instanceof Error ? e.message : String(e));

export class CommandManager<Ctx> {
  readonly #options: CommandManagerOptions<Ctx>;
  readonly #commands = observableState<readonly CommandDefinition<Ctx>[]>([]);
  readonly #groups = observableState<readonly CommandGroup[]>([]);
  /** Each command's registration (a token per `register` call). */
  readonly #owners = new Map<CommandDefinition<Ctx>, object>();
  readonly #active = observableState<CommandInvocation<Ctx> | undefined>(undefined);
  readonly #running = observableState<CommandRun<Ctx> | undefined>(undefined);
  #run: Run<Ctx> | undefined;
  readonly #recents;
  readonly #outcomes = new Subject<CommandOutcome<Ctx>>();

  constructor(options: CommandManagerOptions<Ctx>) {
    this.#options = options;
    let saved: string[] = [];
    try {
      saved = options.recents?.load() ?? [];
    } catch {
      saved = [];
    }
    this.#recents = observableState<readonly string[]>(saved);
  }

  get context(): ReadonlyState<Ctx> {
    return this.#options.context;
  }
  get commands(): ReadonlyState<readonly CommandDefinition<Ctx>[]> {
    return this.#commands.asReadonly;
  }
  /** The open form, if any. */
  get active(): ReadonlyState<CommandInvocation<Ctx> | undefined> {
    return this.#active.asReadonly;
  }
  get recents(): ReadonlyState<readonly string[]> {
    return this.#recents.asReadonly;
  }
  /** The run in flight, if any: nothing else opens until it settles. */
  get running(): ReadonlyState<CommandRun<Ctx> | undefined> {
    return this.#running.asReadonly;
  }
  /** What came of each run, once, as it settles; never replayed to a later subscriber. */
  get outcomes(): Observable<CommandOutcome<Ctx>> {
    return this.#outcomes.asObservable();
  }

  /**
   * Adds commands (a screen, while it is mounted); the returned function takes them away again,
   * cancelling a run of theirs and closing their form. A run that settles after that reports to
   * no one: its screen is gone.
   */
  register(defs: readonly CommandDefinition<Ctx>[]): () => void {
    const owner = {};
    const ids = new Set(defs.map((d) => d.id));
    for (const d of defs) this.#owners.set(d, owner);
    this.#commands.set([...this.#commands.asReadonly.value.filter((c) => !ids.has(c.id)), ...defs]);
    return () => {
      const mine = (c: CommandDefinition<Ctx>) => defs.includes(c) && this.#owners.get(c) === owner;
      const run = this.#run;
      if (run && mine(run.command)) {
        run.cancel();
        this.#endRun(run);
      }
      const active = this.#active.asReadonly.value;
      if (active && mine(active.command)) {
        active.close();
        this.#active.set(undefined);
      }
      for (const d of defs) if (this.#owners.get(d) === owner) this.#owners.delete(d);
      this.#commands.set(this.#commands.asReadonly.value.filter((c) => !defs.includes(c)));
    };
  }

  /** The command tree's groups, for browsing. */
  get groups(): ReadonlyState<readonly CommandGroup[]> {
    return this.#groups.asReadonly;
  }

  /** Adds groups to the command tree (a contribution's, while active); the returned function takes them away. */
  registerGroups(groups: readonly CommandGroup[]): () => void {
    const ids = new Set(groups.map((g) => g.id));
    this.#groups.set([...this.#groups.asReadonly.value.filter((g) => !ids.has(g.id)), ...groups]);
    return () => this.#groups.set(this.#groups.asReadonly.value.filter((g) => !groups.includes(g)));
  }

  /** A group's path from the top level to it (itself last); empty for none or an unknown id. */
  groupPath(id: string | undefined): CommandGroup[] {
    const all = this.#groups.asReadonly.value;
    const path: CommandGroup[] = [];
    let g = id === undefined ? undefined : all.find((x) => x.id === id);
    while (g && !path.includes(g)) {
      path.unshift(g);
      g = g.parent === undefined ? undefined : all.find((x) => x.id === g!.parent);
    }
    return path;
  }

  /** What a group holds (or the top level, for none): its subgroups and its commands, in order. */
  children(id: string | undefined): { groups: CommandGroup[]; commands: CommandMatch<Ctx>[] } {
    const known = new Set(this.#groups.asReadonly.value.map((g) => g.id));
    const byOrder = (a: CommandGroup, b: CommandGroup) => (a.order ?? 0) - (b.order ?? 0);
    const groups = this.#groups.asReadonly.value
      .filter((g) => (g.parent && known.has(g.parent) ? g.parent : undefined) === id)
      .sort(byOrder);
    const commands = this.#commands.asReadonly.value
      .filter((c) => (c.group && known.has(c.group) ? c.group : undefined) === id)
      .map((c) => this.#match(c));
    return { groups, commands };
  }

  /** Is the command inside the group (at any depth)? */
  inGroup(command: CommandDefinition<Ctx>, id: string): boolean {
    return this.groupPath(command.group).some((g) => g.id === id);
  }

  get(id: string): CommandDefinition<Ctx> | undefined {
    return this.#commands.asReadonly.value.find((c) => c.id === id);
  }

  /** Why the command cannot run now (its `when`, or the error it threw), if it cannot. */
  blocked(command: CommandDefinition<Ctx>): string | undefined {
    try {
      const w = command.when?.(this.#options.context.value);
      return typeof w === 'string' ? w : undefined;
    } catch (e) {
      return message(e);
    }
  }

  #match(command: CommandDefinition<Ctx>): CommandMatch<Ctx> {
    const blocked = this.blocked(command);
    return { command, ...(blocked ? { blocked } : {}) };
  }

  /** What the typed text (after "/") matches. */
  query(text: string): CommandQuery<Ctx> {
    const all = this.#commands.asReadonly.value;
    const q = text.trim();
    const lower = text.toLowerCase().trimStart().replace(/\s+/g, ' ');
    // Still the start of a title ("focus on" of "Focus on…") or of an alias ("add scene" of "add scene",
    // not "add" with "scene"): a query, not a command with arguments.
    const typed = lower.trimEnd();
    const titled = all.some(
      (c) =>
        c.title.toLowerCase().replace(/…$/, '').startsWith(typed) ||
        aliases(c).some((a) => a.toLowerCase().startsWith(typed)),
    );
    // The longest alias that fits; among commands sharing it, one that can run.
    const fitting = titled
      ? []
      : all
          .flatMap((c) => (c.params ? aliases(c).map((a) => ({ c, a: a.toLowerCase() })) : []))
          .filter(({ a }) => lower.startsWith(`${a} `) && lower.trim().length > a.length)
          .sort((x, y) => y.a.length - x.a.length);
    const longest = fitting.filter((x) => x.a.length === fitting[0]?.a.length);
    const direct = longest.find((x) => !this.blocked(x.c)) ?? longest[0];
    if (direct) {
      const tail = text
        .trimStart()
        .replace(/\s+/g, ' ')
        .slice(direct.a.length + 1)
        .trim();
      const match = this.#match(direct.c);
      // A blocked one keeps the other ways in sight: what the whole text matches.
      const others = match.blocked ? this.#ranked(all, q).filter((m) => m.command !== direct.c) : [];
      return { recent: [], matches: [match, ...others], direct: { ...match, tail } };
    }
    if (!q) {
      const recent = this.#recents.asReadonly.value
        .map((id) => all.find((c) => c.id === id))
        .filter((c): c is CommandDefinition<Ctx> => !!c)
        .map((c) => this.#match(c));
      const seen = new Set(recent.map((m) => m.command.id));
      return { recent, matches: all.filter((c) => !seen.has(c.id)).map((c) => this.#match(c)) };
    }
    return { recent: [], matches: this.#ranked(all, q) };
  }

  /**
   * Commands whose title or an alias fuzzy-matches the text, best first; and their presets that the
   * text says more about than the command alone (by the preset's title, its words, or both titles).
   */
  #ranked(all: readonly CommandDefinition<Ctx>[], q: string): CommandMatch<Ctx>[] {
    const score = (text: string, less = 0) => (fuzzyScore(text, q) ?? -1) - less;
    const hits = all.flatMap((c) => {
      const s = Math.max(score(c.title), ...aliases(c).map((a) => score(a, 0.5)));
      const presets = this.#presets(c).flatMap((p) => {
        const own = Math.max(score(p.title), ...(p.keywords ?? []).map((k) => score(k, 0.5)));
        const both = score(`${c.title} ${p.title}`);
        return own >= 0 || (s < 0 && both >= 0) ? [{ c, p, s: Math.max(own, both) }] : [];
      });
      return [...(s >= 0 ? [{ c, p: undefined, s }] : []), ...presets];
    });
    // Better matches first (by whole points: a title that starts with the text, a word that does, then looser ones);
    // equal matches keep their category together, in the order categories were registered, then go by title, so
    // "Select …" commands sit together and "Select with …" ones next to each other.
    const categories = new Map<string, number>();
    for (const c of all) if (!categories.has(c.category)) categories.set(c.category, categories.size);
    const title = (x: (typeof hits)[number]) => (x.p ? `${x.c.title} ${x.p.title}` : x.c.title);
    return hits
      .sort(
        (a, b) =>
          Math.floor(b.s) - Math.floor(a.s) ||
          (categories.get(a.c.category) ?? 0) - (categories.get(b.c.category) ?? 0) ||
          title(a).localeCompare(title(b)),
      )
      .map((x) => ({ ...this.#match(x.c), ...(x.p ? { preset: x.p } : {}) }));
  }

  #presets(command: CommandDefinition<Ctx>): readonly CommandPreset[] {
    try {
      return command.presets?.(this.#options.context.value) ?? [];
    } catch {
      return [];
    }
  }

  /**
   * Opens a command: its form (with a preset's values, when one is named), or runs it at once
   * when it has none. Refused while a run is in flight, and when the command cannot run; never
   * throws.
   */
  open(id: string, tail?: string, preset?: string): OpenResult<Ctx> {
    const command = this.get(id);
    if (!command) return { kind: 'blocked', reason: `No command "${id}".` };
    const values = preset === undefined ? undefined : this.#presets(command).find((p) => p.id === preset)?.values;
    if (preset !== undefined && !values) return { kind: 'blocked', reason: `${command.title} has no "${preset}".` };
    if (this.#run) return { kind: 'blocked', reason: `${this.#run.command.title.replace(/…$/, '')} is still running.` };
    const blocked = this.blocked(command);
    if (blocked) return { kind: 'blocked', reason: blocked };
    this.close();
    if (!command.params) return { kind: 'ran', result: this.#runNow(command) };
    let invocation: CommandInvocation<Ctx>;
    try {
      invocation = new CommandInvocation<Ctx>({
        command,
        context: this.#options.context,
        formats: this.#options.formats,
        ...(this.#options.preview ? { preview: this.#options.preview } : {}),
        ...(this.#options.assistant ? { assistant: this.#options.assistant } : {}),
        ...(this.#options.previewDeps ? { previewDeps: this.#options.previewDeps } : {}),
        ...(tail ? { tail } : {}),
        ...(values ? { preset: values } : {}),
        onRun: () =>
          this.#startRun({ command, owner: this.#owners.get(command), invocation, cancel: () => invocation.cancel() }),
        onSettled: (outcome) => this.#settled(invocation, outcome),
      });
    } catch (e) {
      return { kind: 'blocked', reason: message(e) };
    }
    this.#active.set(invocation);
    return { kind: 'form', invocation };
  }

  /** Closes the open form. False while it runs: cancel it instead. */
  close(): boolean {
    const active = this.#active.asReadonly.value;
    if (!active) return true;
    if (this.#run?.invocation === active) return false;
    active.close();
    this.#active.set(undefined);
    return true;
  }

  /** Cancels the run in flight (an explicit act); its outcome follows when it settles. */
  cancel(): boolean {
    const run = this.#run;
    if (!run) return false;
    if (!this.#running.asReadonly.value?.cancelling) this.#running.set({ command: run.command, cancelling: true });
    run.cancel();
    return true;
  }

  /**
   * Runs a command's keybinding, if the event is one it may take. True when it handled the
   * event; a refusal is a `blocked` outcome.
   */
  handleKey(e: KeyEventLike, mac = true): boolean {
    if (e.defaultPrevented || e.repeat || e.isComposing || e.key === 'Process') return false;
    const typing = typesText(e.composedPath?.()[0] ?? e.target);
    const hits = this.#commands.asReadonly.value.filter(
      (c) => c.shortcut && (!typing || typesNothing(c.shortcut)) && chordMatches(c.shortcut, e, mac),
    );
    const command = hits.find((c) => !this.blocked(c)) ?? hits[0];
    if (!command) return false;
    const r = this.open(command.id);
    if (r.kind === 'blocked') this.#outcomes.next({ kind: 'blocked', command, message: r.reason });
    return true;
  }

  async #runNow(command: CommandDefinition<Ctx>): Promise<CommandResult | undefined> {
    const ctx = this.#options.context.value;
    const abort = new AbortController();
    const run: Run<Ctx> = { command, owner: this.#owners.get(command), cancel: () => abort.abort() };
    this.#startRun(run);
    try {
      const values = command.defaults?.(ctx) ?? {};
      const problem = command.check?.(values, ctx);
      if (problem) {
        this.#endRun(run);
        this.#report(run, { kind: 'failed', command, message: problem });
        return undefined;
      }
      const result = await command.run(values, ctx, abort.signal);
      this.#endRun(run);
      this.#report(run, { kind: 'done', command, result });
      return result;
    } catch (e) {
      this.#endRun(run);
      this.#report(
        run,
        abort.signal.aborted ? { kind: 'cancelled', command } : { kind: 'failed', command, message: message(e) },
      );
      return undefined;
    }
  }

  #startRun(run: Run<Ctx>): void {
    this.#run = run;
    this.#running.set({ command: run.command, cancelling: false });
  }

  #endRun(run: Run<Ctx>): void {
    if (this.#run !== run) return;
    this.#run = undefined;
    this.#running.set(undefined);
  }

  /** A form settled: the run ends, its outcome is told, and it closes unless it stays to show an error. */
  #settled(invocation: CommandInvocation<Ctx>, outcome: CommandOutcome<Ctx>): void {
    const run = this.#run?.invocation === invocation ? this.#run : undefined;
    if (run) this.#endRun(run);
    if (outcome.kind !== 'failed' && this.#active.asReadonly.value === invocation) {
      invocation.close();
      this.#active.set(undefined);
    }
    this.#report(
      run ?? { command: outcome.command, owner: this.#owners.get(outcome.command), cancel: () => {} },
      outcome,
    );
  }

  #report(run: Run<Ctx>, outcome: CommandOutcome<Ctx>): void {
    if (outcome.kind === 'done') this.#remember(outcome.command);
    // Its screen went while it ran: there is no one to tell.
    if (!run.owner || this.#owners.get(run.command) !== run.owner) return;
    this.#outcomes.next(outcome);
  }

  #remember(command: CommandDefinition<Ctx>): void {
    const recents = [command.id, ...this.#recents.asReadonly.value.filter((id) => id !== command.id)].slice(
      0,
      this.#options.maxRecents ?? 5,
    );
    this.#recents.set(recents);
    try {
      this.#options.recents?.save(recents);
    } catch {
      // Recents are a convenience; a store that cannot save is not an error.
    }
  }
}

function aliases(c: { alias?: string | readonly string[] }): readonly string[] {
  return c.alias === undefined ? [] : typeof c.alias === 'string' ? [c.alias] : c.alias;
}

const NOT_TEXT = new Set(['button', 'checkbox', 'radio', 'range', 'color', 'file', 'image', 'reset', 'submit']);

/** Does the target take typed text (a text input, a textarea, a select, an editable element)? */
function typesText(target: unknown): boolean {
  if (!target || typeof target !== 'object') return false;
  const t = target as {
    tagName?: string;
    type?: string;
    isContentEditable?: boolean;
    getAttribute?(name: string): string | null;
  };
  if (t.isContentEditable) return true;
  const editable = typeof t.getAttribute === 'function' ? t.getAttribute('contenteditable') : null;
  if (editable !== null && editable !== 'false') return true;
  const tag = t.tagName?.toUpperCase();
  if (tag === 'TEXTAREA' || tag === 'SELECT') return true;
  return tag === 'INPUT' && !NOT_TEXT.has((t.type ?? 'text').toLowerCase());
}

/** ⌘ or Ctrl without Alt never types a character (Alt does on macOS; Ctrl+Alt is AltGr). */
function typesNothing(chord: KeyChord): boolean {
  return !!chord.mod && !chord.alt;
}

/** A single character that is not a letter, a digit or a space. */
function isSymbol(key: string): boolean {
  return /^[^\p{L}\p{N}\s]$/u.test(key);
}

function chordMatches(chord: KeyChord, e: KeyEventLike, mac: boolean): boolean {
  const mod = mac ? !!e.metaKey : !!e.ctrlKey;
  const other = mac ? !!e.ctrlKey : !!e.metaKey;
  if (!!chord.mod !== mod || other || !!chord.alt !== !!e.altKey) return false;
  // A symbol (`?` is Shift+/ on most layouts) already implies its modifiers; letters and digits match Shift strictly.
  if (!!chord.shift !== !!e.shiftKey && (chord.shift || !isSymbol(chord.key))) return false;
  const want = chord.key.toLowerCase();
  const key = e.key.toLowerCase();
  if (key === want) return true;
  // A key that is not a plain letter or digit (⌥F on macOS is "ƒ"; a non-Latin layout): its code.
  if (/^[a-z0-9]$/.test(key) || !e.code) return false;
  const m = /^(?:Key([A-Z])|Digit([0-9]))$/.exec(e.code);
  return !!m && (m[1] ?? m[2] ?? '').toLowerCase() === want;
}
