/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The command manager (commands §4): what is registered, what matches, what is open, what runs and what came of it,
// the menu's highlighted row, and the keys.

import { type Observable, type ReadonlyState, ReactiveModel, State, Subject } from '@molstar/core/reactive';
import { type CommandBrowse, CommandBrowser } from './browse.js';
import { CommandInvocation } from './invocation.js';
import { errorMessage } from './invocation-model.js';
import { isMac, isPlainChord, type KeyEventLike } from './keys.js';
import { KeyIndex, Keymap, type KeymapOptions } from './keymap.js';
import { type CommandMatch, type CommandQuery, queryCommands } from './query.js';
import { CommandRegistry, type LazyCommands, type LazyRegistration } from './registry.js';
import type {
  CommandAssistant,
  CommandDefinition,
  CommandFormats,
  CommandGroup,
  CommandOutcome,
  CommandPreset,
  CommandResult,
  CommandValues,
  KeyChord,
  PreviewSink,
  RecentsStore,
} from './types.js';

export interface CommandManagerOptions<Ctx> {
  /** What commands see: the app, the selection, the focused view… */
  context: ReadonlyState<Ctx>;
  formats: CommandFormats<Ctx>;
  preview?: PreviewSink;
  assistant?: CommandAssistant<Ctx>;
  recents?: RecentsStore;
  /** How many recent commands are kept (default 5). */
  maxRecents?: number;
  /**
   * What in the context previews read: an open invocation's preview is drawn again when one of them changes
   * (`Object.is`, item by item). Without it, every change of the context draws it again.
   */
  previewDeps?: (ctx: Ctx) => readonly unknown[];
  /** The app's keys and the store of the person's remaps (commands §6.3). */
  keymap?: KeymapOptions;
  /** Ms a highlighted row stays before it is previewed and counted (`browse`, default 80). */
  browseDelay?: number;
}

export type OpenResult<Ctx> =
  | { kind: 'form'; invocation: CommandInvocation<Ctx> }
  | { kind: 'ran'; result: Promise<CommandResult | undefined> }
  | { kind: 'blocked'; reason: string };

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

export class CommandManager<Ctx> extends ReactiveModel {
  readonly #options: CommandManagerOptions<Ctx>;
  readonly #registry: CommandRegistry<Ctx>;
  /** Which key runs which command: the commands' own, the app's, the person's (commands §6.3). */
  readonly keymap: Keymap;
  readonly #keys: KeyIndex<Ctx>;
  readonly #active = new State<CommandInvocation<Ctx> | undefined>(undefined);
  readonly #running = new State<CommandRun<Ctx> | undefined>(undefined);
  #run: Run<Ctx> | undefined;
  readonly #recents: State<readonly string[]>;
  readonly #outcomes = new Subject<CommandOutcome<Ctx>>();
  /** The groups' shared values the person picked, by group id. */
  readonly #shared = new State<Readonly<Record<string, string>>>({});
  readonly #browser: CommandBrowser<Ctx>;
  readonly #blocked = (c: CommandDefinition<Ctx>) => this.blocked(c);

  constructor(options: CommandManagerOptions<Ctx>) {
    super();
    this.#options = options;
    this.keymap = new Keymap(options.keymap);
    this.#keys = new KeyIndex(this.keymap);
    this.#registry = new CommandRegistry<Ctx>((gone) => this.#removed(gone));
    this.#browser = new CommandBrowser<Ctx>({
      context: options.context,
      preview: options.preview,
      delay: options.browseDelay ?? 80,
      blocked: this.#blocked,
      args: (c, ctx) => this.#argsFor(c, ctx),
    });
    let saved: string[] = [];
    try {
      saved = options.recents?.load() ?? [];
    } catch {
      saved = [];
    }
    this.#recents = new State<readonly string[]>(saved);
  }

  get context(): ReadonlyState<Ctx> {
    return this.#options.context;
  }
  get commands(): ReadonlyState<readonly CommandDefinition<Ctx>[]> {
    return this.#registry.commands;
  }
  get groups(): ReadonlyState<readonly CommandGroup[]> {
    return this.#registry.groups;
  }
  /** The open invocation, if any. */
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

  // --- Registration (commands §4.1–4.2) ----------------------------------------------------------------

  /**
   * Adds commands (a screen's, while it is mounted) as one change; the returned function takes them away, cancelling a
   * run of theirs and closing their invocation. A run that settles after that reports to no one. An id already
   * registered (or repeated in `defs`) throws, naming it, and nothing is added.
   */
  register(defs: readonly CommandDefinition<Ctx>[]): () => void {
    return this.#registry.register(defs);
  }

  /** Adds groups to the command tree; the returned function takes them away. A taken id throws, adding nothing. */
  registerGroups(groups: readonly CommandGroup[]): () => void {
    return this.#registry.registerGroups(groups);
  }

  /** Adds commands made when wanted (commands §4.2). */
  registerLazy(source: LazyCommands<Ctx>): LazyRegistration {
    return this.#registry.registerLazy(source);
  }

  /** The lists are wanted (the menu, while open) until the returned release: lazy commands are made in a microtask. */
  want(): () => void {
    return this.#registry.want();
  }

  get(id: string): CommandDefinition<Ctx> | undefined {
    return this.#registry.get(id);
  }

  /** A group's path from the top level to it (itself last); empty for none or an unknown id. */
  groupPath(id: string | undefined): CommandGroup[] {
    return this.#registry.groupPath(id);
  }

  /** What a group holds (or the top level, for none): its subgroups and its commands, in order. */
  children(id: string | undefined): { groups: CommandGroup[]; commands: CommandMatch<Ctx>[] } {
    const { groups, commands } = this.#registry.children(id);
    return { groups, commands: commands.map((c) => this.#match(c)) };
  }

  /** Is the command inside the group (at any depth)? */
  inGroup(command: CommandDefinition<Ctx>, id: string): boolean {
    return this.groupPath(command.group).some((g) => g.id === id);
  }

  /** Why the command cannot run now (its `when`, or the error it threw), if it cannot. */
  blocked(command: CommandDefinition<Ctx>): string | undefined {
    try {
      const w = command.when?.(this.#options.context.value);
      return typeof w === 'string' ? w : undefined;
    } catch (e) {
      return errorMessage(e);
    }
  }

  /** The command's key now, after the app's and the person's keymap (commands §6.3). */
  shortcut(command: CommandDefinition<Ctx>): KeyChord | undefined {
    return this.keymap.shortcut(command);
  }

  // --- Shared arguments (commands §2.3) --------------------------------------------------------------

  /** The groups' shared values the person picked, by group id (the defaults are the groups'). */
  get shared(): ReadonlyState<Readonly<Record<string, string>>> {
    return this.#shared.asReadonly;
  }

  /** A group's shared value now: the person's pick, else its default; none for a group without one. */
  sharedValue(groupId: string): string | undefined {
    const g = this.#registry.groupsById().get(groupId);
    if (!g?.shared) return undefined;
    const v = this.#shared.value[groupId];
    return v !== undefined && g.shared.options.some(([o]) => o === v) ? v : g.shared.default;
  }

  /** Sets a group's shared value, for every command in it from now on. */
  setShared(groupId: string, value: string): void {
    if (this.#shared.value[groupId] !== value) this.#shared.set({ ...this.#shared.value, [groupId]: value });
  }

  /** The shared values a command gets from the groups it is in (the nearest group's wins for a name). */
  #sharedFor(command: CommandDefinition<Ctx>): Record<string, string> {
    const out: Record<string, string> = {};
    for (const g of this.groupPath(command.group)) {
      const v = g.shared ? this.sharedValue(g.id) : undefined;
      if (g.shared && v !== undefined) out[g.shared.name] = v;
    }
    return out;
  }

  /** What a command without params runs with: its defaults and its groups' shared values. */
  #argsFor(command: CommandDefinition<Ctx>, ctx: Ctx): Record<string, unknown> {
    return { ...(command.defaults?.(ctx) ?? {}), ...this.#sharedFor(command) };
  }

  // --- Finding (commands §4.3) ----------------------------------------------------------------------

  /** What the typed text matches: recent commands (nothing typed), matches with presets, a direct alias. */
  query(text: string): CommandQuery<Ctx> {
    this.#registry.fresh();
    return queryCommands(this.#registry.commands.value, text, {
      recents: this.#recents.value,
      byId: this.#registry.byId(),
      match: (c, preset) => this.#match(c, preset),
      presets: (c) => this.#presets(c),
    });
  }

  #match(command: CommandDefinition<Ctx>, preset?: CommandPreset): CommandMatch<Ctx> {
    const blocked = this.blocked(command);
    return { command, ...(blocked ? { blocked } : {}), ...(preset ? { preset } : {}) };
  }

  #presets(command: CommandDefinition<Ctx>): readonly CommandPreset[] {
    try {
      return command.presets?.(this.#options.context.value) ?? [];
    } catch {
      return [];
    }
  }

  // --- Browsing (commands §4.5) ---------------------------------------------------------------------

  /** The menu's highlighted command (`browse`), and its count once there is one. */
  get browsing(): ReadonlyState<CommandBrowse | undefined> {
    return this.#browser.browsing;
  }

  /**
   * The menu's highlighted row (commands §4.5): a command without params that can run is previewed and counted after
   * `browseDelay`, for that row only; a move of the highlight (or none) aborts that and clears the preview.
   */
  browse(id: string | undefined): void {
    this.#browser.browse(id === undefined ? undefined : this.get(id));
  }

  // --- Opening and running (commands §4.4) ------------------------------------------------------------

  /**
   * Opens a command: its invocation (with a preset's values when one is named, or the values given: a row it is for),
   * or runs it at once when it has no params. Refused while a run is in flight and when the command cannot run; never
   * throws.
   */
  open(id: string, tail?: string, preset?: string | CommandValues): OpenResult<Ctx> {
    const command = this.get(id);
    if (!command) return { kind: 'blocked', reason: `No command "${id}".` };
    const values =
      preset === undefined
        ? undefined
        : typeof preset === 'string'
          ? this.#presets(command).find((p) => p.id === preset)?.values
          : preset;
    if (typeof preset === 'string' && !values)
      return { kind: 'blocked', reason: `${command.title} has no "${preset}".` };
    if (this.#run) return { kind: 'blocked', reason: `${this.#run.command.title.replace(/…$/, '')} is still running.` };
    const blocked = this.blocked(command);
    if (blocked) return { kind: 'blocked', reason: blocked };
    this.#browser.end();
    this.close();
    if (!command.params) return { kind: 'ran', result: this.#runNow(command) };
    const options = this.#options;
    const start = this.#sharedStart(command, values);
    let invocation: CommandInvocation<Ctx>;
    try {
      invocation = new CommandInvocation<Ctx>({
        command,
        context: options.context,
        formats: options.formats,
        ...(options.preview ? { preview: options.preview } : {}),
        ...(options.assistant ? { assistant: options.assistant } : {}),
        ...(options.previewDeps ? { previewDeps: options.previewDeps } : {}),
        ...(tail ? { tail } : {}),
        ...(start ? { preset: start } : {}),
        onRun: () =>
          this.#startRun({
            command,
            owner: this.#registry.ownerOf(command),
            invocation,
            cancel: () => invocation.cancel(),
          }),
        onSettled: (outcome) => this.#settled(invocation, outcome),
      });
    } catch (e) {
      return { kind: 'blocked', reason: errorMessage(e) };
    }
    this.#active.set(invocation);
    return { kind: 'form', invocation };
  }

  /** An invocation's start: the groups' shared values its params take, under the values it is opened with. */
  #sharedStart(command: CommandDefinition<Ctx>, values: CommandValues | undefined): CommandValues | undefined {
    const shared = this.#sharedFor(command);
    if (!Object.keys(shared).length) return values;
    let params: object;
    try {
      params = command.params?.(this.#options.context.value, values ?? {}) ?? {};
    } catch {
      return values;
    }
    const fit = Object.fromEntries(Object.entries(shared).filter(([k]) => k in params));
    return Object.keys(fit).length ? { ...fit, ...values } : values;
  }

  /** Closes the open invocation. False while it runs: cancel it instead. */
  close(): boolean {
    const active = this.#active.value;
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
    if (!this.#running.value?.cancelling) this.#running.set({ command: run.command, cancelling: true });
    run.cancel();
    return true;
  }

  // --- Keys (commands §6.1) ---------------------------------------------------------------------------

  /**
   * Runs a command's key from a keydown the host forwards while the menu is closed; true when it took the event. A
   * refused command is a `blocked` outcome, except on a plain key, which is then not taken.
   */
  handleKey(e: KeyEventLike, mac = isMac()): boolean {
    if (e.defaultPrevented || e.repeat || e.isComposing || e.key === 'Process') return false;
    const registry = this.#registry;
    // A key the app or the person bound may name a lazy command that is not made yet: made first. An unbound key
    // makes nothing (this runs at every keystroke).
    if (registry.anyStale && this.#keys.bound(registry.commands.value, e, mac)) registry.fresh();
    const hit = this.#keys.find(registry.commands.value, e, mac, this.#blocked);
    if (!hit) return false;
    const { command } = hit;
    // A plain key (an arrow: the next snapshot) is ambient: where its command cannot run, the key is not taken.
    if (isPlainChord(hit.chord) && this.blocked(command)) return false;
    const r = this.open(command.id);
    if (r.kind === 'blocked') this.#outcomes.next({ kind: 'blocked', command, message: r.reason });
    return true;
  }

  // --- Internals -------------------------------------------------------------------------------------

  /** Commands going away: a run of theirs is cancelled (it reports to no one), their invocation closed. */
  #removed(gone: ReadonlySet<CommandDefinition<Ctx>>): void {
    const run = this.#run;
    if (run && gone.has(run.command)) {
      run.cancel();
      this.#endRun(run);
    }
    const active = this.#active.value;
    if (active && gone.has(active.command)) {
      active.close();
      this.#active.set(undefined);
    }
    const browsed = this.#browser.command;
    if (browsed && gone.has(browsed)) this.#browser.end();
  }

  async #runNow(command: CommandDefinition<Ctx>): Promise<CommandResult | undefined> {
    const ctx = this.#options.context.value;
    const abort = new AbortController();
    const run: Run<Ctx> = { command, owner: this.#registry.ownerOf(command), cancel: () => abort.abort() };
    this.#startRun(run);
    try {
      const args = this.#argsFor(command, ctx);
      const problem = command.check?.(args, ctx);
      if (problem) {
        this.#endRun(run);
        this.#report(run, { kind: 'failed', command, message: problem });
        return undefined;
      }
      const result = await command.run(args, ctx, abort.signal);
      this.#endRun(run);
      this.#report(run, { kind: 'done', command, result });
      return result;
    } catch (e) {
      this.#endRun(run);
      this.#report(
        run,
        abort.signal.aborted ? { kind: 'cancelled', command } : { kind: 'failed', command, message: errorMessage(e) },
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

  /** An invocation settled: the run ends, its outcome is told, and it closes unless it stays to show an error. */
  #settled(invocation: CommandInvocation<Ctx>, outcome: CommandOutcome<Ctx>): void {
    const run = this.#run?.invocation === invocation ? this.#run : undefined;
    if (run) this.#endRun(run);
    if (outcome.kind !== 'failed' && this.#active.value === invocation) {
      invocation.close();
      this.#active.set(undefined);
    }
    this.#report(
      run ?? { command: outcome.command, owner: this.#registry.ownerOf(outcome.command), cancel: () => {} },
      outcome,
    );
  }

  #report(run: Run<Ctx>, outcome: CommandOutcome<Ctx>): void {
    if (outcome.kind === 'done') this.#remember(outcome.command);
    // Its registration went while it ran: there is no one to tell.
    if (!run.owner || this.#registry.ownerOf(run.command) !== run.owner) return;
    this.#outcomes.next(outcome);
  }

  #remember(command: CommandDefinition<Ctx>): void {
    const recents = [command.id, ...this.#recents.value.filter((id) => id !== command.id)].slice(
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

  /** Retired for good: the open invocation closes, a run is cancelled, browsing ends, `outcomes` completes. */
  protected override onRetire(): void {
    const run = this.#run;
    if (run) {
      run.cancel();
      this.#endRun(run);
    }
    this.#active.value?.close();
    this.#active.set(undefined);
    this.#browser.end();
    this.#outcomes.complete();
  }
}
