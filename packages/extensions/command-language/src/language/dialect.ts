/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Diagnostic, Span } from './diagnostics.js';
import type { Dialect } from './ir.js';
import type { CommandHelp, Compiled, Compiler } from './table.js';

// A dialect (_spec §3): how its text splits into commands, and how one compiles. PyMOL's and VMD's implement it; the
// driver (`script.ts`) and completion are the same for both.

/** One command of a submission, as split: its name, place and text (the dialect keeps its own parse beside). */
export interface Unit {
  readonly name: string;
  readonly at: Span;
  readonly source: string;
}

export interface DialectImpl {
  readonly id: Dialect;
  /** The table, for help and completion. */
  readonly commands: readonly CommandHelp[];
  split(text: string): { readonly units: readonly Unit[]; readonly diagnostics: readonly Diagnostic[] };
  /** Compiles one unit; throws `ScriptError` for an error, reports warnings through `c`. */
  compile(unit: Unit, c: Compiler): Compiled;
  /** A unit that runs another script given with the submission (`@file.pml`, `source file.tcl`): its name. */
  includeOf?(unit: Unit): { readonly name: string; readonly at: Span } | undefined;
}
