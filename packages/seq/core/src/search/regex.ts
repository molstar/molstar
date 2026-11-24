/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A regular expression's search (_spec/README.md "Search", Regular expressions): a subset of POSIX ERE read here (letters,
// `.`, classes, groups, alternation, `?`, `*`, `+`, counted repeats, `^` and `$` at the ends, `\` to take a character
// literally), compiled to a Glushkov automaton (one state a letter position, no ε-moves) and simulated bit-parallel
// (Navarro and Raffinot's tables: a table a byte of the state). Never backtracking: a pattern is untrusted, and the cost
// is the letters times the words whatever the pattern. Every bound is checked before a table is allocated.

import type { Issue } from '@molstar/core/repr/tree';
import { issueOf } from '@molstar/core/repr/tree';
import { SearchIssues } from './issues.js';
import { EdgeEnd, EdgeStart, type HitVisitor, type SearchBounds, type SearchKernel, type StepBudget } from './types.js';

/** A parsed expression: a set of codes, a concatenation, an alternation, a repeat (`max` −1 unbounded), nothing. */
type Ast =
  | { readonly t: 'set'; readonly set: Uint8Array; readonly letter?: number }
  | { readonly t: 'cat'; readonly xs: readonly Ast[] }
  | { readonly t: 'alt'; readonly xs: readonly Ast[] }
  | { readonly t: 'rep'; readonly x: Ast; readonly min: number; readonly max: number }
  | { readonly t: 'eps' };

/** A regular expression read: its tree, its anchors, its literal letters (the alphabet check), its positions. */
export interface RegexParsed {
  readonly ast: Ast;
  readonly anchorStart: boolean;
  readonly anchorEnd: boolean;
  readonly letters: readonly number[];
  readonly positions: number;
  /** The longest match it can make (−1: unbounded). */
  readonly maxMatch: number;
}

const upper = (c: number) => (c >= 97 && c <= 122 ? c - 32 : c);
const Unsupported = 'seq.search.pattern.unsupported';

class Refused {
  constructor(readonly issue: Issue) {}
}

/** Reads a pattern, or says why it does not read (each issue with its offset in the pattern). */
export function parseRegex(pattern: string, bounds: SearchBounds): RegexParsed | readonly Issue[] {
  if (pattern.length > bounds.maxRegexChars)
    return [
      issueOf(
        SearchIssues,
        'seq.search.pattern.size',
        `The pattern has ${pattern.length} characters; a regular expression has at most ${bounds.maxRegexChars}`,
      ),
    ];
  let i = 0;
  let depth = 0;
  const letters: number[] = [];
  const fail = (code: 'seq.search.pattern.unsupported' | 'seq.search.pattern.syntax', message: string, at = i) => {
    throw new Refused(issueOf(SearchIssues, code, message, { offset: at }));
  };
  for (let k = 0; k < pattern.length; k++)
    if (pattern.charCodeAt(k) > 255) fail(Unsupported, 'Only Latin letters are searched', k);
  const anchorStart = pattern.startsWith('^');
  const end = pattern.length - (pattern.endsWith('$') && !pattern.endsWith('\\$') ? 1 : 0);
  if (anchorStart) i = 1;

  const literal = (c: number): Ast => {
    const set = new Uint8Array(256);
    set[upper(c)] = 1;
    letters.push(upper(c));
    return { t: 'set', set, letter: upper(c) };
  };
  const escaped = (): Ast => {
    const c = pattern[i + 1];
    if (c === undefined || i + 1 >= end) fail('seq.search.pattern.syntax', 'A pattern ends with \\');
    if (/[0-9]/.test(c)) fail(Unsupported, 'Backreferences are not searched');
    if (/[a-zA-Z]/.test(c)) fail(Unsupported, `\\${c} is not searched: write the letters, or a class`);
    i += 2;
    return literal(c.charCodeAt(0));
  };
  const klass = (): Ast => {
    const start = i;
    i++;
    const neg = pattern[i] === '^';
    if (neg) i++;
    const set = new Uint8Array(256);
    let first = true;
    while (i < end && (pattern[i] !== ']' || first)) {
      first = false;
      let a = pattern.charCodeAt(i);
      if (pattern[i] === '\\') {
        const e = pattern[i + 1];
        if (e === undefined) break;
        if (/[a-zA-Z0-9]/.test(e)) fail(Unsupported, `\\${e} is not searched in a class`);
        a = e.charCodeAt(0);
        i++;
      } else if (pattern[i] === '[' && /^\[[:=.]/.test(pattern.slice(i, i + 2)))
        fail(Unsupported, 'Named classes are not searched');
      i++;
      let b = a;
      if (pattern[i] === '-' && i + 1 < end && pattern[i + 1] !== ']') {
        b = pattern.charCodeAt(i + 1);
        if (pattern[i + 1] === '\\') fail(Unsupported, 'An escaped range end is not searched');
        if (b < a) fail('seq.search.pattern.syntax', 'A class range runs backward');
        i += 2;
      }
      for (let c = a; c <= b; c++) set[upper(c)] = 1;
    }
    if (pattern[i] !== ']' || i >= end) fail('seq.search.pattern.syntax', 'A class is not closed', start);
    i++;
    if (neg) for (let c = 0; c < 256; c++) set[c] = set[c] ? 0 : 1;
    return { t: 'set', set };
  };
  const count = (): { min: number; max: number } | undefined => {
    const m = /^\{(\d{1,6})(,(\d{0,6}))?\}/.exec(pattern.slice(i, end));
    if (!m) return undefined;
    const min = Number(m[1]);
    const max = m[2] === undefined ? min : m[3] === '' ? -1 : Number(m[3]);
    if (max !== -1 && max < min) fail('seq.search.pattern.syntax', 'A count’s maximum is below its minimum');
    i += m[0].length;
    return { min, max };
  };
  const atom = (): Ast | undefined => {
    const c = pattern[i];
    if (c === '(') {
      const start = i;
      if (pattern[i + 1] === '?') {
        if (pattern[i + 2] !== ':') fail(Unsupported, 'Lookaround, named groups and flags are not searched');
        i += 3;
      } else i++;
      if (++depth > bounds.maxRegexDepth)
        throw new Refused(
          issueOf(SearchIssues, 'seq.search.pattern.size', `Groups nest deeper than ${bounds.maxRegexDepth}`, {
            offset: start,
          }),
        );
      const x = alt();
      if (pattern[i] !== ')' || i >= end) fail('seq.search.pattern.syntax', 'A group is not closed', start);
      i++;
      depth--;
      return x;
    }
    if (c === '[') return klass();
    if (c === '\\') return escaped();
    if (c === '.') {
      i++;
      return { t: 'set', set: new Uint8Array(256).fill(1) };
    }
    if (c === '^' || c === '$')
      fail(Unsupported, `${c} is searched only at the pattern’s ${c === '^' ? 'start' : 'end'}`);
    if (c === '*' || c === '+' || c === '?') fail('seq.search.pattern.syntax', `${c} repeats nothing`);
    if (c === '{' && count()) fail('seq.search.pattern.syntax', 'A count repeats nothing');
    if (c === ')' || c === '|' || c === undefined) return undefined;
    i++;
    return literal(c.charCodeAt(0));
  };
  const repeat = (): Ast | undefined => {
    let x = atom();
    if (!x) return undefined;
    for (;;) {
      const c = pattern[i];
      if (i >= end) break;
      let q: { min: number; max: number } | undefined;
      if (c === '*') q = { min: 0, max: -1 };
      else if (c === '+') q = { min: 1, max: -1 };
      else if (c === '?') q = { min: 0, max: 1 };
      if (q) i++;
      else if (c === '{') q = count();
      if (!q) break;
      if (i < end && (pattern[i] === '?' || pattern[i] === '+'))
        fail(Unsupported, 'Lazy and possessive repeats are not searched: a hit is the longest match from its start');
      x = { t: 'rep', x, min: q.min, max: q.max };
    }
    return x;
  };
  const cat = (): Ast => {
    const xs: Ast[] = [];
    while (i < end && pattern[i] !== '|' && pattern[i] !== ')') {
      const x = repeat();
      if (!x) break;
      xs.push(x);
    }
    return xs.length === 1 ? xs[0] : xs.length ? { t: 'cat', xs } : { t: 'eps' };
  };
  function alt(): Ast {
    const xs = [cat()];
    while (i < end && pattern[i] === '|') {
      i++;
      xs.push(cat());
    }
    return xs.length === 1 ? xs[0] : { t: 'alt', xs };
  }
  try {
    const ast = alt();
    if (i < end)
      fail('seq.search.pattern.syntax', pattern[i] === ')' ? 'A ) closes no group' : 'The pattern does not read');
    const positions = positionsOf(ast, bounds.maxRegexPositions + 1);
    if (positions > bounds.maxRegexPositions)
      return [
        issueOf(
          SearchIssues,
          'seq.search.pattern.size',
          `The pattern takes more than ${bounds.maxRegexPositions} letter positions once its repeats are counted out`,
        ),
      ];
    if (positions === 0 || maxMatchOf(ast) === 0)
      return [issueOf(SearchIssues, 'seq.search.pattern.empty', 'The pattern matches only the empty string')];
    return { ast, anchorStart, anchorEnd: end < pattern.length, letters, positions, maxMatch: maxMatchOf(ast) };
  } catch (e) {
    if (e instanceof Refused) return [e.issue];
    throw e;
  }
}

/** Letter positions once repeats are counted out, stopping past `cap` (a count of 10^6 copies is never multiplied out). */
function positionsOf(a: Ast, cap: number): number {
  switch (a.t) {
    case 'set':
      return 1;
    case 'eps':
      return 0;
    case 'cat':
    case 'alt': {
      let n = 0;
      for (const x of a.xs) if ((n += positionsOf(x, cap)) > cap) return n;
      return n;
    }
    case 'rep': {
      const one = positionsOf(a.x, cap);
      const copies = a.max === -1 ? Math.max(a.min, 1) : a.max;
      return Math.min(one * copies, cap + 1);
    }
  }
}

/** The longest match, −1 when unbounded. */
function maxMatchOf(a: Ast): number {
  switch (a.t) {
    case 'set':
      return 1;
    case 'eps':
      return 0;
    case 'cat': {
      let n = 0;
      for (const x of a.xs) {
        const m = maxMatchOf(x);
        if (m < 0) return -1;
        n += m;
      }
      return n;
    }
    case 'alt': {
      let n = 0;
      for (const x of a.xs) {
        const m = maxMatchOf(x);
        if (m < 0) return -1;
        n = Math.max(n, m);
      }
      return n;
    }
    case 'rep': {
      const m = maxMatchOf(a.x);
      if (m === 0) return 0;
      return m < 0 || a.max === -1 ? -1 : m * a.max;
    }
  }
}

/** A Glushkov automaton: each position's set, `first`, `last` and `follow`, as bitsets of `words` words. */
export interface Glushkov {
  readonly positions: number;
  readonly words: number;
  readonly sets: readonly Uint8Array[];
  readonly first: Int32Array;
  readonly last: Int32Array;
  readonly follow: readonly Int32Array[];
}

/** The automaton of a parsed expression, counted repeats copied out (within the checked positions). */
export function glushkov(parsed: RegexParsed): Glushkov {
  const words = Math.max(1, Math.ceil(parsed.positions / 32));
  const sets: Uint8Array[] = [];
  const follow: Int32Array[] = [];
  const bits = () => new Int32Array(words);
  const or = (a: Int32Array, b: Int32Array) => {
    for (let k = 0; k < words; k++) a[k] |= b[k];
  };
  const each = (s: Int32Array, f: (p: number) => void) => {
    for (let k = 0; k < words; k++) for (let b = s[k]; b !== 0; b &= b - 1) f(k * 32 + (31 - Math.clz32(b & -b)));
  };
  interface Part {
    nullable: boolean;
    first: Int32Array;
    last: Int32Array;
  }
  const seq = (parts: Part[]): Part => {
    let acc: Part = { nullable: true, first: bits(), last: bits() };
    for (const p of parts) {
      each(acc.last, (q) => or(follow[q], p.first));
      const first = acc.first.slice();
      if (acc.nullable) or(first, p.first);
      const last = p.last.slice();
      if (p.nullable) or(last, acc.last);
      acc = { nullable: acc.nullable && p.nullable, first, last };
    }
    return acc;
  };
  const loop = (p: Part, nullable: boolean): Part => {
    each(p.last, (q) => or(follow[q], p.first));
    return { nullable: nullable || p.nullable, first: p.first, last: p.last };
  };
  const build = (a: Ast): Part => {
    switch (a.t) {
      case 'set': {
        const p = sets.length;
        sets.push(a.set);
        follow.push(bits());
        const s = bits();
        s[p >> 5] |= 1 << (p & 31);
        return { nullable: false, first: s, last: s.slice() };
      }
      case 'eps':
        return { nullable: true, first: bits(), last: bits() };
      case 'cat':
        return seq(a.xs.map(build));
      case 'alt': {
        const out: Part = { nullable: false, first: bits(), last: bits() };
        for (const x of a.xs) {
          const p = build(x);
          out.nullable ||= p.nullable;
          or(out.first, p.first);
          or(out.last, p.last);
        }
        return out;
      }
      case 'rep': {
        const parts: Part[] = [];
        if (a.max === -1) {
          for (let k = 1; k < a.min; k++) parts.push(build(a.x));
          parts.push(loop(build(a.x), a.min === 0));
        } else {
          for (let k = 0; k < a.min; k++) parts.push(build(a.x));
          for (let k = a.min; k < a.max; k++) {
            const p = build(a.x);
            parts.push({ ...p, nullable: true });
          }
        }
        return seq(parts);
      }
    }
  };
  const top = build(parsed.ast);
  return { positions: sets.length, words, sets, first: top.first, last: top.last, follow };
}

/** The tables a step reads: `T(D)` as a table a byte of the state, of `words` words each. */
function followTable(g: Glushkov, follow: readonly Int32Array[]): Int32Array {
  const w = g.words;
  const chunks = Math.ceil(g.positions / 8);
  const t = new Int32Array(chunks * 256 * w);
  for (let k = 0; k < chunks; k++)
    for (let byte = 1; byte < 256; byte++) {
      const low = byte & -byte;
      const p = k * 8 + (31 - Math.clz32(low));
      const at = (k * 256 + byte) * w;
      const rest = (k * 256 + (byte & (byte - 1))) * w;
      for (let x = 0; x < w; x++) t[at + x] = t[rest + x] | (p < follow.length ? follow[p][x] : 0);
    }
  return t;
}

/** `out = T(d)`: the union of the follows of `d`'s positions, a lookup a byte. */
function step(t: Int32Array, w: number, chunks: number, d: Int32Array, out: Int32Array): void {
  out.fill(0);
  for (let k = 0; k < chunks; k++) {
    const byte = (d[k >> 2] >>> ((k & 3) * 8)) & 255;
    if (byte === 0) continue;
    const at = (k * 256 + byte) * w;
    for (let x = 0; x < w; x++) out[x] |= t[at + x];
  }
}

/** The automaton's tables, shared by the kernels of one compile (read only). */
export interface RegexTables {
  readonly g: Glushkov;
  readonly forward: Int32Array;
  readonly backward: Int32Array;
}

export function regexTables(g: Glushkov): RegexTables {
  // the reversed automaton's follow: q follows p backward when p follows q forward
  const back: Int32Array[] = g.follow.map(() => new Int32Array(g.words));
  for (let p = 0; p < g.positions; p++)
    for (let q = 0; q < g.positions; q++) if (g.follow[p][q >> 5] & (1 << (q & 31))) back[q][p >> 5] |= 1 << (p & 31);
  return { g, forward: followTable(g, g.follow), backward: followTable(g, back) };
}

/**
 * The kernel of an automaton over one alphabet's codes (`accepts(set, code)`: whether a position's set takes a folded
 * sequence code). A hit a start, the longest match from it, at most `maxLength` long: a backward pass of the reversed
 * automaton marks where a match starts, then a bounded forward pass from each start finds its end, counted against the
 * budget (a pattern long at every start ends at the budget, never a hang).
 */
export function regexKernel(
  tables: RegexTables,
  parsed: RegexParsed,
  maxLength: number,
  accepts: (set: Uint8Array, code: number) => boolean,
): SearchKernel {
  const { g, forward, backward } = tables;
  const w = g.words;
  const chunks = Math.ceil(g.positions / 8);
  const masks = new Int32Array(256 * w);
  for (let c = 0; c < 256; c++)
    for (let p = 0; p < g.positions; p++) if (accepts(g.sets[p], c)) masks[c * w + (p >> 5)] |= 1 << (p & 31);
  const d = new Int32Array(w);
  const e = new Int32Array(w);
  let marks = new Uint8Array(0);
  const meets = (a: Int32Array, b: Int32Array) => {
    for (let x = 0; x < w; x++) if ((a[x] & b[x]) !== 0) return true;
    return false;
  };

  /** The longest match's end from `s`, or −1; its letters charged. */
  const longest = (seq: Uint8Array, n: number, s: number, atEnd: boolean, budget: StepBudget): number => {
    const lim = Math.min(n, s + maxLength);
    let c = seq[s] * w;
    let any = 0;
    for (let x = 0; x < w; x++) {
      d[x] = g.first[x] & masks[c + x];
      any |= d[x];
    }
    let out = -1;
    let u = s + 1;
    if (any !== 0) {
      if (meets(d, g.last) && (!parsed.anchorEnd || (u === n && atEnd))) out = u;
      for (; u < lim; u++) {
        step(forward, w, chunks, d, e);
        c = seq[u] * w;
        any = 0;
        for (let x = 0; x < w; x++) {
          d[x] = e[x] & masks[c + x];
          any |= d[x];
        }
        if (any === 0) break;
        if (meets(d, g.last) && (!parsed.anchorEnd || (u + 1 === n && atEnd))) out = u + 1;
      }
    }
    budget.steps += (u - s) * w;
    return out;
  };

  return {
    maxLength,
    words: w,
    find(seq: Uint8Array, n: number, from: number, to: number, edges: number, budget: StepBudget, visit: HitVisitor) {
      const atStart = (edges & EdgeStart) !== 0;
      const atEnd = (edges & EdgeEnd) !== 0;
      if (parsed.anchorStart && (!atStart || from > 0)) return true;
      if (parsed.anchorEnd && !atEnd) return true;
      if (marks.length < n) marks = new Uint8Array(Math.max(n, marks.length * 2));
      marks.fill(0, from, to);
      // the reversed automaton over the letters backward, a match started again at each letter (`$`: only at the end)
      d.fill(0);
      for (let u = n - 1; u >= from; u--) {
        step(backward, w, chunks, d, e);
        if (!parsed.anchorEnd || u === n - 1) for (let x = 0; x < w; x++) e[x] |= g.last[x];
        const c = seq[u] * w;
        for (let x = 0; x < w; x++) d[x] = e[x] & masks[c + x];
        if (u < to && meets(d, g.first)) marks[u] = 1;
      }
      const last = parsed.anchorStart ? Math.min(to, from + 1) : to;
      for (let s = from; s < last; s++) {
        if (!marks[s]) continue;
        const end = longest(seq, n, s, atEnd, budget);
        if (end > s) visit(s, end);
        if (budget.steps > budget.limit) return false;
      }
      return true;
    },
  };
}
