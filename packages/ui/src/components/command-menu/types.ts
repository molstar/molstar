/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the command menu takes (commands §5): its props, a command's own form, a field editor by format, the assistant.

import type { ReadonlyState } from '@molstar/core/reactive';
import type { CommandInvocation } from '../../commands/invocation.js';
import type { CommandManager } from '../../commands/manager.js';
import type { FieldChoices } from '../../commands/sentence.js';
import type { CommandMenuPrefix } from '../../commands/types.js';
import type { Child } from '../../dom/index.js';

export interface CommandMenuProps {
  manager: CommandManager<any>;
  /** Closes the menu (Esc on the list, a press on the scrim, a command that ran). */
  onClose(): void;
  /** Text to start with ("select "). */
  initial?: string;
  /** Prefixes that turn the menu into another list (`@`, `#`, `pymol `); `?` help is built in. */
  prefixes?: readonly CommandMenuPrefix[];
  /** Opens inside this group. */
  group?: string;
  /** The assistant: typed words get an Ask row. */
  assistant?: CommandMenuAssistant;
  /** Opens the app's help (a `HelpView`): the first row of `?`, then the menu closes. */
  help?: () => void;
  /** Opened for one command from elsewhere: its sentence is the menu, so Esc (and ⌫ past its start) closes the menu. */
  direct?: boolean;
  /** Commands' own forms, by command id (commands §5.4). */
  forms?: Readonly<Record<string, CommandFormView>>;
  /** Field editors by format name (commands §5.4): drawn in place of the choices while a field of that format is active. */
  fieldViews?: Readonly<Record<string, CommandFieldView>>;
}

/** What a command's own form and a field editor get: the invocation to edit, and the input's keys first. */
export interface CommandViewProps {
  manager: CommandManager<any>;
  invocation: CommandInvocation<any>;
  /** Runs the command as Enter does. */
  run(): void;
  /** A key in the menu's input, before the menu acts on it: true takes it. */
  keys(handler: (e: KeyboardEvent) => boolean): void;
}

/** A command's own form (commands §5.4), under the sentence's choices. */
export interface CommandFormView {
  /** The fields it draws itself: while one is active, the menu lists no choices for it. */
  fields?: readonly string[];
  /** `false`: it draws the params past the sentence itself; the Options fold is not shown. */
  options?: false;
  render(props: CommandViewProps): Child;
}

/** What a field editor gets: the field it edits and its choices, besides the invocation. */
export interface CommandFieldViewProps extends CommandViewProps {
  /** The param name of the field. */
  field: string;
  /** The field's choices as the layer lists them (filtered by what is typed). */
  choices: ReadonlyState<FieldChoices | undefined>;
  /** Sets the field's value (and moves on, as taking a choice does). */
  set(value: string | number | boolean): void;
}

/** An editor for a kind of field many commands share (a color's swatches), by format name (commands §5.4). */
export interface CommandFieldView {
  render(props: CommandFieldViewProps): Child;
}

/**
 * What the menu asks of an assistant: words read into one command (its sentence opens filled; nothing runs until Enter),
 * or the request handed on; without a model, the setup.
 */
export interface CommandMenuAssistant {
  ready(): boolean;
  /** A command these words fit, with them read as its arguments, if one does. */
  propose(text: string): { command: string; tail?: string } | undefined;
  /** Hands the request on (a turn of the assistant's own). */
  ask(text: string): void;
  /** Opens the model setup. */
  setup(): void;
}
