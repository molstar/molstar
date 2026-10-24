/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Diagnostic, error, type Span } from '../diagnostics.js';

// Tcl's word rules (_spec §3.2), for VMD's commands: words split on whitespace; `{…}` is literal (nested braces),
// `"…"` substitutes; `$name`, `${name}` and `[…]` substitute in quoted and bare words; `\` escapes; `;` and newlines
// end a command; `#` at a command's start is a comment; `\` before a newline joins the lines. What the words mean is
// the dialect's (`dialect.ts`): only a closed subset of Tcl is evaluated.

export type TclPart =
  | { readonly kind: 'text'; readonly text: string }
  | { readonly kind: 'var'; readonly name: string; readonly at: Span }
  /** A command substitution: its script, from `from` in the submission. */
  | { readonly kind: 'command'; readonly script: string; readonly from: number; readonly at: Span };

export interface TclWord {
  readonly at: Span;
  /** `{…}`: literal, no substitution. */
  readonly braced: boolean;
  readonly parts: readonly TclPart[];
}

export interface TclCommand {
  readonly words: readonly TclWord[];
  readonly at: Span;
  readonly source: string;
  /** Ends with a space: completion offers the next word. */
  readonly open: boolean;
}

const Escapes: Readonly<Record<string, string>> = {
  n: '\n',
  t: '\t',
  r: '\r',
  '\\': '\\',
  '"': '"',
  $: '$',
  '[': '[',
  ']': ']',
  '{': '{',
  '}': '}',
  ' ': ' ',
  ';': ';',
};

class Reader {
  i: number;
  constructor(
    readonly text: string,
    start: number,
    readonly end: number,
    readonly diagnostics: Diagnostic[],
  ) {
    this.i = start;
  }
  get done() {
    return this.i >= this.end;
  }
  peek(k = 0) {
    return this.text[this.i + k];
  }
  /** Skips spaces and tabs (and `\` newline joins) inside a command. */
  skipBlanks() {
    for (;;) {
      const c = this.peek();
      if (c === ' ' || c === '\t' || c === '\r') this.i++;
      else if (c === '\\' && this.peek(1) === '\n') this.i += 2;
      else break;
    }
  }
  /** The index of the `]` matching a `[` at `i`, or -1. */
  matchBracket(from: number): number {
    let depth = 0;
    let brace = 0;
    let quote = false;
    for (let k = from; k < this.end; k++) {
      const c = this.text[k];
      if (c === '\\') {
        k++;
        continue;
      }
      if (quote) {
        if (c === '"') quote = false;
        continue;
      }
      if (brace) {
        if (c === '{') brace++;
        else if (c === '}') brace--;
        continue;
      }
      if (c === '"') quote = true;
      else if (c === '{') brace++;
      else if (c === '[') depth++;
      else if (c === ']') {
        depth--;
        if (depth === 0) return k;
      }
    }
    return -1;
  }
}

/** `$name`, `${name}` or `[…]` at the reader: a part, or undefined when it is plain text. */
function substitution(r: Reader): TclPart | undefined {
  const c = r.peek();
  const from = r.i;
  if (c === '$') {
    if (r.peek(1) === '{') {
      const close = r.text.indexOf('}', r.i + 2);
      if (close < 0 || close >= r.end) {
        r.diagnostics.push(error('syntax', 'The variable name is not closed.', { from, to: r.end }));
        r.i = r.end;
        return { kind: 'text', text: '' };
      }
      const name = r.text.slice(r.i + 2, close);
      r.i = close + 1;
      return { kind: 'var', name, at: { from, to: r.i } };
    }
    const m = r.text.slice(r.i + 1, r.end).match(/^[A-Za-z0-9_]+(?:\([^)]*\))?/);
    if (!m) return undefined;
    r.i += 1 + m[0].length;
    return { kind: 'var', name: m[0], at: { from, to: r.i } };
  }
  if (c === '[') {
    const close = r.matchBracket(r.i);
    if (close < 0) {
      r.diagnostics.push(error('syntax', 'The bracket is not closed.', { from, to: r.end }));
      r.i = r.end;
      return { kind: 'text', text: '' };
    }
    const part: TclPart = {
      kind: 'command',
      script: r.text.slice(r.i + 1, close),
      from: r.i + 1,
      at: { from, to: close + 1 },
    };
    r.i = close + 1;
    return part;
  }
  return undefined;
}

/** Reads one word at the reader (it is at the word's first character). */
function readWord(r: Reader): TclWord {
  const from = r.i;
  const parts: TclPart[] = [];
  let text = '';
  const flush = () => {
    if (text) parts.push({ kind: 'text', text });
    text = '';
  };
  const c = r.peek();
  if (c === '{') {
    let depth = 0;
    let k = r.i;
    for (; k < r.end; k++) {
      const ch = r.text[k];
      if (ch === '\\') {
        k++;
        continue;
      }
      if (ch === '{') depth++;
      else if (ch === '}') {
        depth--;
        if (depth === 0) break;
      }
    }
    if (k >= r.end) {
      r.diagnostics.push(error('syntax', 'The brace is not closed.', { from, to: r.end }));
      const body = r.text.slice(r.i + 1, r.end);
      r.i = r.end;
      return { at: { from, to: r.end }, braced: true, parts: [{ kind: 'text', text: body }] };
    }
    const body = r.text.slice(r.i + 1, k).replace(/\\\n\s*/g, ' ');
    r.i = k + 1;
    const next = r.peek();
    if (!r.done && next !== ' ' && next !== '\t' && next !== '\n' && next !== ';' && next !== '\r') {
      r.diagnostics.push(error('syntax', 'Extra characters after the closing brace.', { from: r.i, to: r.i + 1 }));
    }
    return { at: { from, to: r.i }, braced: true, parts: [{ kind: 'text', text: body }] };
  }
  const quoted = c === '"';
  if (quoted) r.i++;
  while (!r.done) {
    const ch = r.peek();
    if (quoted && ch === '"') {
      r.i++;
      flush();
      return { at: { from, to: r.i }, braced: false, parts };
    }
    if (!quoted && (ch === ' ' || ch === '\t' || ch === '\n' || ch === ';' || ch === '\r')) break;
    if (ch === '\\') {
      const nx = r.peek(1);
      if (nx === '\n') {
        text += ' ';
        r.i += 2;
        while (r.peek() === ' ' || r.peek() === '\t') r.i++;
      } else {
        text += nx === undefined ? '\\' : (Escapes[nx] ?? nx);
        r.i += 2;
      }
      continue;
    }
    if (ch === '$' || ch === '[') {
      const part = substitution(r);
      if (part) {
        flush();
        parts.push(part);
        continue;
      }
    }
    text += ch;
    r.i++;
  }
  if (quoted) r.diagnostics.push(error('syntax', 'The quote is not closed.', { from, to: r.i }));
  flush();
  return { at: { from, to: r.i }, braced: false, parts };
}

/** Splits Tcl text (a submission, or a substitution's script) into commands of words. */
export function splitTcl(
  text: string,
  from = 0,
  to = text.length,
): { commands: TclCommand[]; diagnostics: Diagnostic[] } {
  const diagnostics: Diagnostic[] = [];
  const commands: TclCommand[] = [];
  const r = new Reader(text, from, to, diagnostics);
  while (!r.done) {
    // Between commands: blanks, newlines, semicolons.
    for (;;) {
      r.skipBlanks();
      const c = r.peek();
      if (c === '\n' || c === ';') r.i++;
      else break;
    }
    if (r.done) break;
    if (r.peek() === '#') {
      // A comment runs to the end of its line (a `\` newline continues it).
      while (!r.done && r.peek() !== '\n') r.i += r.peek() === '\\' ? 2 : 1;
      continue;
    }
    const start = r.i;
    const words: TclWord[] = [];
    while (!r.done) {
      r.skipBlanks();
      const c = r.peek();
      if (r.done || c === '\n' || c === ';') break;
      words.push(readWord(r));
    }
    let end = r.i;
    while (end > start && /\s/.test(text[end - 1])) end--;
    if (words.length) {
      commands.push({
        words,
        at: { from: start, to: end },
        source: text.slice(start, end),
        open: r.i > (words[words.length - 1]?.at.to ?? 0),
      });
    }
  }
  return { commands, diagnostics };
}

/** A word's text when it has no substitutions (a command name, a keyword); undefined otherwise. */
export function literal(w: TclWord): string | undefined {
  if (w.parts.every((p) => p.kind === 'text')) return w.parts.map((p) => (p as { text: string }).text).join('');
  return undefined;
}
