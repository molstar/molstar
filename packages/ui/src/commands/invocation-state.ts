/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An invocation's state and how it is made (commands §3.1).

import type { ReadonlyState } from '@molstar/core/reactive';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { FormFieldView } from '../forms/flat-fields.js';
import type {
  CommandAssistant,
  CommandDefinition,
  CommandFormats,
  CommandOutcome,
  CommandResult,
  CommandValues,
  PreviewSink,
} from './types.js';

export type InvocationStatus = 'editing' | 'resolving' | 'running' | 'cancelling' | 'done' | 'closed';

export interface InvocationError {
  /** The field it is about, if one. */
  field?: string;
  message: string;
  hint?: string;
}

export interface InvocationState {
  status: InvocationStatus;
  /** All the params now (they follow the fields that reshape them). */
  params: PD.Params;
  /** The sentence's fields, in order. */
  fields: readonly FormFieldView[];
  /** The params the sentence does not hold (its options), in order. */
  formOnly: readonly string[];
  /** The run action's label: "Apply". */
  action: string;
  /** The sentence's values set. */
  values: CommandValues;
  /** The options set (form-only and hidden params): what a run gets over their defaults. */
  optionValues: Readonly<Record<string, unknown>>;
  /** The field being typed into. */
  field?: string;
  /** What is typed into it (not yet a value). */
  text: string;
  /** The highlighted choice (`choices.items`); -1 for none. */
  pick: number;
  /** Fields whose values the assistant proposed and nobody has touched since. */
  proposed: readonly string[];
  /** Required fields still empty. */
  missing: readonly string[];
  error?: InvocationError;
  /** The typed text is not this field's grammar, and an assistant could read it. */
  canAssist: boolean;
  /** What the assistant said about its proposal. */
  note?: string;
  /** Words the assistant could not express as this command. */
  escalated?: string;
  /** The command cannot run now (its `when`). */
  blocked?: string;
  result?: CommandResult;
}

export interface InvocationInit<Ctx> {
  command: CommandDefinition<Ctx>;
  context: ReadonlyState<Ctx>;
  formats: CommandFormats<Ctx>;
  preview?: PreviewSink;
  assistant?: CommandAssistant<Ctx>;
  /** Text typed after the command's alias: "sticks red". */
  tail?: string;
  /** Values it opens with, over the defaults (a preset's, a shared argument's). */
  preset?: CommandValues;
  /** What in the context the preview reads; compared item by item. */
  previewDeps?: (ctx: Ctx) => readonly unknown[];
  /** The run started (its effect may land from now on). */
  onRun?: () => void;
  /** How it ended: once per started run, and for an answer from the assistant. */
  onSettled?: (outcome: CommandOutcome<Ctx>) => void;
}
