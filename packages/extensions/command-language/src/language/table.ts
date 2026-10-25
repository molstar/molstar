/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Diagnostic, Span } from './diagnostics.js';
import type { CompileState } from './env.js';
import type { CommandBody, Dialect } from './ir.js';

// The command tables (_spec §2): one per dialect, each command's name, signature, help and examples, and how it
// compiles. Parsing, completion, `help`, the docs and the spec's scope tables read them.

/** What an argument is: completion offers by it, and the compiler reads it. */
export type ArgKind =
  | 'selection'
  | 'rep'
  | 'color'
  | 'number'
  | 'int'
  | 'bool'
  /** A new name (a selection's, an object's). */
  | 'name'
  | 'object'
  | 'molecule'
  | 'file'
  | 'setting'
  /** A setting's value, read by the setting. */
  | 'value'
  | 'enum'
  | 'text';

export interface ParamSpec {
  readonly name: string;
  readonly kind: ArgKind;
  readonly optional?: boolean;
  /** An `enum`'s words (also offered for others where a few words are the usual ones). */
  readonly values?: readonly string[];
}

/** A command's help: what `help`, completion and the docs show. */
export interface CommandHelp {
  /** As typed (`show`, `mol addrep`). */
  readonly name: string;
  readonly aliases?: readonly string[];
  /** One line: what it does in Mol*. */
  readonly summary: string;
  /** The forms it takes (`show [rep [, selection]]`). */
  readonly usage: string;
  readonly params?: readonly ParamSpec[];
  /** Lines that compile (a test runs them all). */
  readonly examples: readonly string[];
}

/** A real command of a dialect that Mol* does not do: why, and what to do instead. */
export interface Unsupported {
  readonly name: string;
  readonly why: string;
  readonly hint?: string;
}

/** What a command's `compile` gets besides its arguments. */
export interface Compiler {
  readonly dialect: Dialect;
  readonly state: CompileState;
  /** The command's place and text. */
  readonly at: Span;
  readonly source: string;
  /** A warning or note about this command (it still runs). */
  report(d: Diagnostic): void;
}

/** What a command compiles to: IR commands without their place (the compiler adds it). */
export type Compiled = readonly CommandBody[];
