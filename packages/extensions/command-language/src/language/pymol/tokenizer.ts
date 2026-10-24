/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Diagnostic, error, type Span } from '../diagnostics.js';

// PyMOL's command lines (_spec §3.1), by hand so every piece keeps its place and an unfinished line still reads:
// `command arg, arg, key=value`. The first space (or comma) ends the command's name; commas at depth 0 (outside
// parentheses, brackets and quotes) split the arguments; `;` and newlines end a command; `#` starts a comment outside
// quotes; `\` at a line's end continues it. A `python … python end` block is one unsupported command.

export interface PymolArg {
  /** The argument's text, trimmed (quotes kept: the command's parameter kind decides). */
  readonly text: string;
  readonly at: Span;
  /** `key=value`: the key, as written; `text` is then the value. */
  readonly key?: string;
  readonly keyAt?: Span;
}

export interface PymolCall {
  readonly name: string;
  readonly nameAt: Span;
  readonly args: readonly PymolArg[];
  /** The whole command, from its name to its end. */
  readonly at: Span;
  /** Its text, as written. */
  readonly source: string;
  /** A trailing comma or space after the last argument: completion offers the next one. */
  readonly open: boolean;
}

export interface PymolSplit {
  readonly calls: readonly PymolCall[];
  readonly diagnostics: readonly Diagnostic[];
}

/** A command's text between `from` and `to` in the submission. */
interface Segment {
  readonly from: number;
  readonly to: number;
}

const isSpace = (c: string) => c === ' ' || c === '\t' || c === '\r';

/** Splits a submission into command segments: `;` and newlines outside quotes, `\` continuations, `#` comments. */
function segments(text: string, diagnostics: Diagnostic[]): Segment[] {
  const out: Segment[] = [];
  let start = 0;
  let quote: string | undefined;
  let quoteAt = 0;
  let i = 0;
  const end = (to: number) => {
    out.push({ from: start, to });
  };
  while (i < text.length) {
    const c = text[i];
    if (quote) {
      if (c === '\\' && i + 1 < text.length && text[i + 1] !== '\n') i += 2;
      else if (c === quote) {
        quote = undefined;
        i++;
      } else if (c === '\n') {
        diagnostics.push(error('syntax', 'The quote is not closed on its line.', { from: quoteAt, to: i }));
        quote = undefined;
        end(i);
        start = ++i;
      } else i++;
      continue;
    }
    if (c === '"' || c === "'") {
      quote = c;
      quoteAt = i;
      i++;
    } else if (c === '#') {
      end(i);
      const nl = text.indexOf('\n', i);
      i = nl < 0 ? text.length : nl;
      start = i;
    } else if (c === ';' || c === '\n') {
      end(i);
      start = ++i;
    } else i++;
  }
  if (quote) diagnostics.push(error('syntax', 'The quote is not closed.', { from: quoteAt, to: text.length }));
  end(text.length);
  return out;
}

/** The argument texts of a command's rest: commas at depth 0 outside quotes. */
function splitArgs(text: string, from: number, to: number): { args: PymolArg[]; open: boolean } {
  const args: PymolArg[] = [];
  if (from >= to) return { args, open: false };
  let depth = 0;
  let quote: string | undefined;
  let start = from;
  const push = (end: number) => {
    let a = start;
    let b = end;
    while (a < b && /\s/.test(text[a])) a++;
    while (b > a && /\s/.test(text[b - 1])) b--;
    const raw = text.slice(a, b);
    const kw = raw.match(/^([A-Za-z_]\w*)\s*=(?!=)\s*/);
    if (kw) {
      args.push({
        text: raw.slice(kw[0].length),
        at: { from: a + kw[0].length, to: b },
        key: kw[1],
        keyAt: { from: a, to: a + kw[1].length },
      });
    } else args.push({ text: raw, at: { from: a, to: b } });
  };
  for (let i = from; i < to; i++) {
    const c = text[i];
    if (quote) {
      if (c === '\\') i++;
      else if (c === quote) quote = undefined;
      continue;
    }
    if (c === '"' || c === "'") quote = c;
    else if (c === '(' || c === '[' || c === '{') depth++;
    else if ((c === ')' || c === ']' || c === '}') && depth > 0) depth--;
    else if (c === ',' && depth === 0) {
      push(i);
      start = i + 1;
    }
  }
  push(to);
  // `show sticks,` and `show ` leave an empty last argument: the person is about to type it.
  const last = args[args.length - 1];
  const open = !!last && !last.text && !last.key;
  if (open || (args.length === 1 && !args[0].text)) args.pop();
  return { args, open: open || /[\s,]$/.test(text.slice(from, to)) };
}

/** Splits a PyMOL submission into commands. */
export function splitPymol(source: string): PymolSplit {
  // A `\` before a newline joins the lines: blanked out (same length, so every place stays the source's).
  const text = source.replace(/\\[ \t]*\r?\n/g, (m) => ' '.repeat(m.length));
  const diagnostics: Diagnostic[] = [];
  const calls: PymolCall[] = [];
  const segs = segments(text, diagnostics);
  for (let k = 0; k < segs.length; k++) {
    let { from } = segs[k];
    const { to } = segs[k];
    while (from < to && /\s/.test(text[from])) from++;
    if (from >= to) continue;
    let end = to;
    while (end > from && /\s/.test(text[end - 1])) end--;
    // The name: `@file`, `util.cbag`, `cmd.show` read as one word.
    const m = text.slice(from, end).match(/^(@\S*|[A-Za-z_][\w.]*)/);
    if (!m) {
      diagnostics.push(error('syntax', `'${text.slice(from, end)}' does not start with a command.`, { from, to: end }));
      continue;
    }
    const name = m[1];
    const nameAt = { from, to: from + name.length };
    // `python` … `python end`: one block, skipped (no Python here).
    if (name === 'python' && !/\bend\s*$/.test(text.slice(nameAt.to, end))) {
      let close = k + 1;
      while (close < segs.length && !/^\s*python\s+end\s*$/.test(text.slice(segs[close].from, segs[close].to))) close++;
      const blockEnd = close < segs.length ? segs[close].to : text.length;
      calls.push({
        name,
        nameAt,
        args: [],
        at: { from, to: blockEnd },
        source: text.slice(from, blockEnd),
        open: false,
      });
      k = close;
      continue;
    }
    let restFrom = nameAt.to;
    while (restFrom < to && isSpace(text[restFrom])) restFrom++;
    if (text[restFrom] === ',') restFrom++;
    const { args, open } = splitArgs(text, restFrom, to);
    calls.push({
      name,
      nameAt,
      args,
      at: { from, to: end },
      source: text.slice(from, end),
      open: open || (args.length === 0 && to > nameAt.to),
    });
  }
  return { calls, diagnostics };
}

/** An argument without its quotes (`'my file.pdb'` → `my file.pdb`). */
export function unquote(text: string): string {
  const t = text.trim();
  if (t.length >= 2 && (t[0] === '"' || t[0] === "'") && t[t.length - 1] === t[0]) return t.slice(1, -1);
  return t;
}
