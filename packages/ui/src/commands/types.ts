/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a command is and the ports the layer reads through (commands §2). DOM-free.

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { IconName } from '../components/icon/icon.js';
import type { FieldHint, FieldValue, FieldValues } from '../forms/flat-fields.js';

/** A command's inline argument values by param name; primitive, as inline fields hold them. */
export type CommandValue = FieldValue;
export type CommandValues = FieldValues;
/** What a run gets: the inline values, the options set, and the params' defaults for the rest. */
export type CommandArgs = Readonly<Record<string, unknown>>;

/** What a command reports when it ran. */
export interface CommandResult {
  /** One line for the person: what happened. */
  message?: string;
  /** More text (an answer, a detail). */
  detail?: string;
  /** The equivalent operation, in mono, when there is one. */
  ops?: string;
  /** It changed something Undo can take back. */
  undoable?: boolean;
  /** Anything else for the host. */
  data?: unknown;
}

/**
 * What happened to a command: one event per settled run, never replayed. `blocked` is a keybinding the manager refused
 * (it has no caller to tell).
 */
export type CommandOutcome<Ctx> =
  | { kind: 'done'; command: CommandDefinition<Ctx>; result: CommandResult }
  | {
      kind: 'failed';
      command: CommandDefinition<Ctx>;
      message: string;
      /** The invocation stays open and shows the error itself. */
      open?: boolean;
    }
  | { kind: 'cancelled'; command: CommandDefinition<Ctx> }
  | { kind: 'blocked'; command: CommandDefinition<Ctx>; message: string };

/** A key chord that runs a command without the menu, e.g. { key: 't', alt: true }. */
export interface KeyChord {
  key: string;
  alt?: boolean;
  shift?: boolean;
  /** ⌘ on macOS, Ctrl elsewhere. */
  mod?: boolean;
}

/** How a param reads inline (commands §2.2): a before-word and a format, and for a number the values it offers. */
export interface CommandFieldHint extends FieldHint {
  /** A number's ready values, listed while the field is active: a value, or [value, label]. */
  presets?: readonly (number | readonly [value: number, label: string])[];
}

export interface CommandDefinition<Ctx> {
  id: string;
  title: string;
  category: string;
  /** Shown before the title in the command menu and on buttons that run it. */
  icon?: IconName;
  /**
   * The group it sits in, by id (`select.expression`), for browsing: the menu opens groups like folders and shows a
   * command's path in search results. Top level without one; a search finds it either way.
   */
  group?: string;
  /** Words that start this command when typed with arguments ("show as" in "show as sticks"); the longest that fits wins. */
  alias?: string | readonly string[];
  /** Its default key; the app's and the person's keymap may replace it (commands §6). */
  shortcut?: KeyChord;
  /**
   * Where its shortcut works (commands §6.2): `global` (the default: anywhere) or a scope id (`canvas`, `panel:<id>`, …)
   * that an element around the focus names with `data-msui-scope`. Inside its scope it wins over a global command.
   */
  scope?: string;
  /** Can it run now? A string says why not; the command stays listed, greyed. */
  when?(ctx: Ctx): true | string;
  /**
   * The arguments, as params. Without them, picking runs the command. Text, string selects, booleans and numbers read
   * inline as a sentence; the rest are its options (commands §3.4). They may follow the values: see `reshapes`.
   */
  params?(ctx: Ctx, values: CommandValues): PD.Params;
  /** How params read inline: before-words, formats and a number's presets, by param name. */
  fields?: Readonly<Record<string, CommandFieldHint>>;
  /** The label of the action that runs it: "Apply". */
  action?: string;
  /** Params whose value changes the params. Set or cleared, they are built again, keeping the values that still fit. */
  reshapes?: readonly string[];
  /** Ways in that a search finds (never listed unasked), each opening the invocation with its values. */
  presets?(ctx: Ctx): readonly CommandPreset[];
  /** Initial values, e.g. from the selection. */
  defaults?(ctx: Ctx): CommandValues;
  /**
   * What it would do, handed to the `PreviewSink` as it is: while its invocation is open and its sentence complete, and,
   * for a command without params, while its row is highlighted in the menu (commands §4.5).
   */
  preview?(args: CommandArgs, ctx: Ctx): unknown;
  /**
   * How much it would do, in a few words ("412 residues"), shown on its highlighted row: computed for that row only, and
   * `signal` aborts when the highlight moves (commands §4.5). Commands without params only.
   */
  count?(args: CommandArgs, ctx: Ctx, signal: AbortSignal): string | undefined | Promise<string | undefined>;
  /** A last check across fields before running; a message blocks the run. */
  check?(args: CommandArgs, ctx: Ctx): string | undefined;
  /** Runs it. Once `signal` aborts (an explicit cancel), settle promptly: reject, or finish if the effect landed. */
  run(args: CommandArgs, ctx: Ctx, signal: AbortSignal): CommandResult | Promise<CommandResult>;
  /** An assistant may fill this command in from words; the instruction says what the values mean. */
  assist?: { instruction: string; answers?: boolean };
  /**
   * Where it is in the UI (a card's section, a tab): the menu offers to go there instead (⌥↵), for what needs a panel.
   * `title` names the place: "Export › Models".
   */
  view?: { title: string; open(ctx: Ctx): void };
}

/** A node of the command tree for browsing (Select › Expression); groups nest by `parent`. */
export interface CommandGroup {
  id: string;
  title: string;
  icon?: IconName;
  /** The group it sits in; top level without one. */
  parent?: string;
  /** Position among its siblings (lower first; registration order otherwise). */
  order?: number;
  /** An argument set once for every command in the group, at any depth (commands §2.3). */
  shared?: SharedArgument;
}

/** A group's shared argument: one of a few named values. */
export interface SharedArgument {
  /** The arg's name (`mode`). */
  name: string;
  /** What it is called in the menu's header ("Mode"). */
  title: string;
  /** Its values, [value, label]. */
  options: readonly (readonly [value: string, label: string])[];
  /** The value until the person picks another. */
  default: string;
}

/** A ready start of a command: "Load structure… from AlphaFold DB". */
export interface CommandPreset {
  /** Unique within its command. */
  id: string;
  /** Said after the command's title: "from AlphaFold DB". */
  title: string;
  values: CommandValues;
  /** Other words it is found by (a source's id). */
  keywords?: readonly string[];
}

// --- Ports ------------------------------------------------------------------------------------

export type FormatCheck =
  | { ok: true; value: string }
  | { ok: false; message: string; hint?: string; unreadable?: boolean };

export interface FormatCompletion {
  value: string;
  label: string;
  detail?: string;
  /**
   * Taking it fills the field's text for editing instead of choosing it (an example expression to adapt): nothing runs
   * until the text is taken as it is.
   */
  fill?: boolean;
}

/** A named input type ("selection"): checks typed text into a value and completes it (commands §2.2). */
export interface CommandFormat<Ctx> {
  /** What can be typed, in a few words (shown while the field is active). */
  hint?(ctx: Ctx): string | undefined;
  check(value: string, ctx: Ctx): FormatCheck;
  complete(text: string, ctx: Ctx): FormatCompletion[];
}

export interface CommandFormats<Ctx> {
  get(name: string): CommandFormat<Ctx> | undefined;
}

/** Where previews go (the host draws them); cleared when the invocation closes or runs. */
export interface PreviewSink {
  set(preview: unknown): void;
  clear(): void;
}

export interface AssistRequest<Ctx> {
  command: CommandDefinition<Ctx>;
  params: PD.Params;
  /** Values already set (they stay unless the assistant replaces them). */
  values: CommandValues;
  /** What the person typed, in words. */
  text: string;
  /** The field it was typed into. */
  field?: string;
  ctx: Ctx;
}

export type AssistResult =
  | { kind: 'values'; values: CommandValues; note?: string }
  | { kind: 'answer'; text: string }
  | { kind: 'escalate'; reason: string };

/** Reads words into a command's values. It proposes; nothing runs until confirmed. */
export interface CommandAssistant<Ctx> {
  propose(request: AssistRequest<Ctx>, signal: AbortSignal): Promise<AssistResult>;
}

/** Keys by command id (commands §6.3): a chord replaces the command's own, `null` takes it away. */
export type KeyBindings = Readonly<Record<string, KeyChord | null>>;

/** Where the person's remapped keys are kept (a host's preferences); a store that throws keeps nothing. */
export interface KeymapStore {
  load(): KeyBindings | undefined;
  save(bindings: KeyBindings): void;
}

/** Recently used command ids, most recent first (e.g. kept in the browser). */
export interface RecentsStore {
  load(): string[];
  save(ids: readonly string[]): void;
}

/** One thing a command-menu prefix offers (a scene node, a snapshot, a completion, …). */
export interface CommandMenuItem {
  id: string;
  label: string;
  detail?: string;
  icon?: IconName;
  /**
   * Does it; return true to keep the menu open (it opened an invocation, or typed another prefix), or `{ text }` to put
   * that text in the menu and keep it open (a completion that edits the line rather than runs it).
   */
  run(): boolean | void | { readonly text: string };
}

/**
 * A prefix that turns the menu into another list (`@` a scene node, `#` a snapshot, `pymol ` a PyMOL command); the host
 * gives them. `char` is matched as the text's start, so a word with its space works.
 */
export interface CommandMenuPrefix {
  char: string;
  title: string;
  /** What to type after it. */
  placeholder?: string;
  /** The rows for the text after the prefix (leading spaces trimmed). */
  items(text: string): readonly CommandMenuItem[];
}
