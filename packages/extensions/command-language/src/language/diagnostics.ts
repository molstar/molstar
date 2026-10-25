/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What went wrong, and where (_spec §7): every message names its place in the submission, so a console marks the line
// and the menu shows the part that did not read. Errors before anything runs (syntax, unknown, unsupported, argument,
// selection) are all reported at once; a runtime failure stops the run at its command.

/** A place in the submission: character offsets, `to` exclusive. */
export interface Span {
  readonly from: number;
  readonly to: number;
}

export type DiagnosticCode =
  /** The text does not read: an unclosed quote or brace, a missing argument. */
  | 'syntax'
  /** Not a command of the dialect. */
  | 'unknown-command'
  /** A real command of the dialect that Mol* does not do (the hint says what it does instead). */
  | 'unsupported'
  /** An option the command takes that Mol* does not (an error when the result would differ, else a warning). */
  | 'unsupported-argument'
  /** An argument of the wrong kind or out of range. */
  | 'argument'
  /** A selection that does not read (the transpiler's message, placed in the line). */
  | 'selection'
  /** A name no object, selection, molecule or variable has. */
  | 'no-name'
  /** The host failed the command (an operation's failure: its code, message and hint). */
  | 'runtime';

export interface Diagnostic {
  readonly severity: 'error' | 'warning' | 'note';
  readonly code: DiagnosticCode;
  readonly message: string;
  /** What to do instead. */
  readonly hint?: string;
  readonly at: Span;
}

/** A diagnostic thrown out of the compiler to stop the command it is about. */
export class ScriptError extends Error {
  constructor(readonly diagnostic: Diagnostic) {
    super(diagnostic.message);
  }
}

export function error(code: DiagnosticCode, message: string, at: Span, hint?: string): Diagnostic {
  return hint ? { severity: 'error', code, message, hint, at } : { severity: 'error', code, message, at };
}

export function warning(code: DiagnosticCode, message: string, at: Span, hint?: string): Diagnostic {
  return hint ? { severity: 'warning', code, message, hint, at } : { severity: 'warning', code, message, at };
}

export function note(message: string, at: Span): Diagnostic {
  return { severity: 'note', code: 'unsupported-argument', message, at };
}

/** Throws an error diagnostic: the compiler's way out of a command. */
export function fail(code: DiagnosticCode, message: string, at: Span, hint?: string): never {
  throw new ScriptError(error(code, message, at, hint));
}

/** The 1-based line and column of an offset. */
export function positionOf(text: string, offset: number): { readonly line: number; readonly column: number } {
  let line = 1;
  let start = 0;
  for (let i = 0; i < offset && i < text.length; i++) {
    if (text.charCodeAt(i) === 10) {
      line++;
      start = i + 1;
    }
  }
  return { line, column: offset - start + 1 };
}

/** One line for a person: "line 2: no command `shw` (did you mean `show`?)". Single-line text omits the line. */
export function formatDiagnostic(text: string, d: Diagnostic): string {
  const multiline = text.includes('\n');
  const where = multiline ? `line ${positionOf(text, d.at.from).line}: ` : '';
  const kind = d.severity === 'error' ? '' : `${d.severity}: `;
  return `${where}${kind}${d.message}${d.hint ? ` ${d.hint}` : ''}`;
}

/** The edit distance of two short words, for "did you mean". */
export function editDistance(a: string, b: string): number {
  const m = a.length;
  const n = b.length;
  let prev = Array.from({ length: n + 1 }, (_, j) => j);
  for (let i = 1; i <= m; i++) {
    const cur = [i];
    for (let j = 1; j <= n; j++) {
      cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1));
    }
    prev = cur;
  }
  return prev[n];
}

/** The nearest of the words, when it is near enough to be a typo. */
export function nearest(word: string, words: Iterable<string>): string | undefined {
  const w = word.toLowerCase();
  let best: string | undefined;
  let bestD = Infinity;
  for (const c of words) {
    const d = editDistance(w, c.toLowerCase());
    if (d < bestD) {
      bestD = d;
      best = c;
    }
  }
  return best !== undefined && bestD <= Math.max(1, Math.floor(w.length / 3)) ? best : undefined;
}
