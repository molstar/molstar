/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Addresses (_spec/README.md "Addresses"; addressing.md): an id is local, unique among its siblings; a node's address
// is the path of ids from the root, `/data/p53/features`, as an IPv6 address is a path of groups. The root is `/`. A
// node's address is a prefix of every address below it, so `/data/*` names a subtree. There is no elision: a search at
// any depth is the query language's.
//
// A segment may also be a query form, `(kind :key value …)`: by default it selects the addressable children of the
// node(s) before it by kind and params (meso*'s child selection); a caller's form resolver table gives a form kind
// another meaning (a range, an item, a rect), whose result may be no nodes at all (locations-and-queries.md).
//
//   address = "/" | ("/" step)+ ["/*"]
//   step    = id | form
//   form    = "(" kind (":" key value)* ")"
//   value   = atom | list | form | address
//   list    = "[" value* "]"
//
// Written from meso*'s tree core. Changed: forms are read here by their own small reader (the tree knows no codec;
// the S-expression codec reads an address word through this one), values may be lists, strings and addresses, atoms
// are atoms.ts's (a word that parses as a JavaScript number is a number), relative addresses wait for definitions, and
// selection takes a form-resolver table.

import { Delimiter, KeywordPattern, printAtom, printKey, readQuoted, readWord } from './atoms.js';
import { valuesEqual } from './hash.js';
import { type Issue, issueOf, TreeIssues } from './issues.js';
import {
  type Atom,
  IdHint,
  IdPattern,
  isNode,
  isReservedKey,
  KindHint,
  KindPattern,
  MaxDepth,
  type Node,
  own,
  type Value,
} from './node.js';

/** The ids from the root, the root itself excluded (`[]` is the root). */
export type Address = readonly string[];

/** One step of an address pattern: an id, or a query form (a node with a kind and params, no id, no children). */
export type Segment = string | Node;

/** An address, or with `subtree` every address it is a prefix of (`/a/*`); with forms among its segments, a set. */
export interface AddressPattern {
  readonly segments: readonly Segment[];
  readonly subtree: boolean;
}

type Code = keyof typeof TreeIssues.codes;
const issue = (code: Code, message: string, offset?: number, hint?: string): Issue =>
  issueOf(TreeIssues, code, message, { offset, hint });

const FormHint = '(kind :key value …): values are atoms, lists, addresses or nested forms; no id, no children';

/** Whether a pattern is only ids: it then names at most one node (or one subtree), found in O(depth × fan-out). */
export function isExactAddress(pattern: AddressPattern): boolean {
  return pattern.segments.every((s) => typeof s === 'string');
}

// ---- printing ----

/** The canonical text, forms on one line. A segment that is not a query form is a programmer error: this throws. */
export function printAddress(segments: readonly Segment[], subtree = false): string {
  let path = '';
  for (const s of segments) {
    if (typeof s === 'string') {
      if (!IdPattern.test(s)) throw new Error('An address segment is not an id');
      path += `/${s}`;
    } else {
      const why = formIssue(s);
      if (why) throw new Error(`An address segment is not a query form: ${why}`);
      path += `/${printForm(s)}`;
    }
  }
  return subtree ? `${path}/*` : path || '/';
}

/** A query form on one line: `(range :from 94 :to 292)`. */
function printForm(form: Node): string {
  let out = `(${form.kind}`;
  const params = form.params ?? {};
  for (const key of Object.keys(params)) {
    if (params[key] === undefined) continue;
    out += ` ${printKey(key)} ${printFormValue(params[key])}`;
  }
  return `${out})`;
}

function printFormValue(v: Value): string {
  if (isNode(v)) return printForm(v);
  if (Array.isArray(v)) return `[${v.map(printFormValue).join(' ')}]`;
  return printAddressAtom(v as Atom);
}

/** An atom as an address or the S-expression codec spells it: a canonical address string bare, others by atoms.ts. */
export function printAddressAtom(v: Atom): string {
  if (typeof v === 'string' && v.startsWith('/') && isCanonicalAddress(v)) return v;
  return printAtom(v);
}

/**
 * Why a value is not a query form, or undefined: a node whose kind is a word, no id, no children, no args, no custom
 * data, no reserved key, values atoms, lists of form values, or nested query forms (no records, no columns).
 */
function formIssue(form: unknown): string | undefined {
  if (!isNode(form) || !KindPattern.test(form.kind)) return 'A query form is a node whose kind is a word';
  if (form.id !== undefined) return 'A query form has no id: an id is a segment of its own';
  if (form.children?.length || form.args?.length || form.custom) return 'A query form has params only';
  const params = form.params ?? {};
  for (const key of Object.keys(params)) {
    if (isReservedKey(key)) return 'A query form names a reserved key';
    const why = formValueIssue(params[key]);
    if (why) return why;
  }
  return undefined;
}

function formValueIssue(v: Value): string | undefined {
  if (isNode(v)) return formIssue(v);
  if (Array.isArray(v)) {
    for (const x of v) {
      const why = formValueIssue(x);
      if (why) return why;
    }
    return undefined;
  }
  return isAtom(v) || v === undefined ? undefined : 'A query form compares atoms, lists and nested forms only';
}

const isAtom = (v: unknown) =>
  v === null || typeof v === 'string' || typeof v === 'number' || typeof v === 'bigint' || typeof v === 'boolean';

/** The issue of a pattern built in code, checked as `parseAddress` checks text, or undefined. */
function patternIssue(pattern: AddressPattern): Issue | undefined {
  const { segments } = pattern;
  if (!Array.isArray(segments) || typeof pattern.subtree !== 'boolean')
    return issue('repr.tree.address.syntax', 'A pattern is {segments, subtree}');
  if (segments.length > MaxDepth) return issue('repr.tree.address.depth_exceeded', `More than ${MaxDepth} segments`);
  for (let i = 0; i < segments.length; i++) {
    const s: unknown = segments[i];
    if (typeof s === 'string') {
      if (!IdPattern.test(s))
        return issue('repr.tree.address.syntax', `Segment ${i + 1} is not an id`, undefined, IdHint);
      continue;
    }
    const why = formIssue(s);
    if (why) return issue('repr.tree.address.form', `${why} (segment ${i + 1})`, undefined, FormHint);
  }
  return undefined;
}

// ---- reading ----

/**
 * `/`, `/a/b`, `/a/b/*`, `/*`, `/data/p53/(range :from 94 :to 292)`; strictly. Issues carry an offset into the text,
 * never the text.
 */
export function parseAddress(text: string): { pattern?: AddressPattern; issues: readonly Issue[] } {
  return parseAddressIn(text, 0, text.length);
}

/** `parseAddress` of `text` between `start` and `end`, offsets into the whole text (an address word in a spec). */
export function parseAddressIn(
  text: string,
  start: number,
  end: number,
  depth = 0,
): { pattern?: AddressPattern; issues: readonly Issue[] } {
  const fail = (code: Code, message: string, offset: number, hint?: string) => ({
    issues: [issue(code, message, offset, hint)],
  });
  if (start >= end || text[start] !== '/') return fail('repr.tree.address.syntax', 'An address opens with /', start);
  if (end - start === 1) return { pattern: { segments: [], subtree: false }, issues: [] };
  const segments: Segment[] = [];
  let subtree = false;
  let i = start;
  while (i < end) {
    const from = ++i; // past a `/`
    if (text[from] === '(') {
      const r = readForm(text, from, end, depth + 1);
      if ('issue' in r) return { issues: [r.issue] };
      if (segments.length === MaxDepth)
        return fail('repr.tree.address.depth_exceeded', `More than ${MaxDepth} segments`, from);
      segments.push(r.form);
      i = r.next;
    } else {
      let next = text.indexOf('/', from);
      if (next < 0 || next > end) next = end;
      const part = text.slice(from, next);
      if (part === '*') {
        if (next !== end) return fail('repr.tree.address.syntax', '* ends an address', from);
        subtree = true;
      } else if (!IdPattern.test(part)) return fail('repr.tree.address.syntax', 'A segment is not an id', from, IdHint);
      else if (segments.length === MaxDepth)
        return fail('repr.tree.address.depth_exceeded', `More than ${MaxDepth} segments`, from);
      else segments.push(part);
      i = next;
    }
    if (i < end && text[i] !== '/') return fail('repr.tree.address.syntax', 'A segment ends at / or at the end', i);
    if (i === end - 1) return fail('repr.tree.address.syntax', 'An address does not end with /', i);
  }
  return { pattern: { segments, subtree }, issues: [] };
}

type Read<T> = { readonly issue: Issue } | ({ readonly next: number } & T);

const isSpace = (c: string) => c === ' ' || c === '\t' || c === '\r' || c === '\n';

function skipSpace(text: string, i: number, end: number): number {
  while (i < end && isSpace(text[i])) i++;
  return i;
}

/** The end of the word at `i`: up to a delimiter. */
function wordEnd(text: string, i: number, end: number): number {
  while (i < end && !Delimiter.test(text[i])) i++;
  return i;
}

/**
 * The end of an address word at `i` (at its `/`): through balanced query forms (strings inside honoured), up to a
 * delimiter outside them; -1 when a form does not close.
 */
export function addressWordEnd(text: string, i: number, end: number): number {
  while (i < end && (text[i] === '(' || !Delimiter.test(text[i]))) {
    if (text[i] !== '(') i++;
    else {
      const close = formEnd(text, i, end);
      if (close < 0) return -1;
      i = close;
    }
  }
  return i;
}

/** The offset just past the form opening at `start`, or -1: brackets balanced outside strings, no comments. */
export function formEnd(text: string, start: number, end: number): number {
  let depth = 0;
  for (let i = start; i < end; i++) {
    const c = text[i];
    if (c === '"') {
      i++;
      while (i < end && text[i] !== '"') i += text[i] === '\\' ? 2 : 1;
      if (i >= end) return -1;
    } else if (c === ';') return -1;
    else if (c === '(') depth++;
    else if (c === ')' && --depth === 0) return i + 1;
  }
  return -1;
}

/** A query form at `i` (at its `(`). */
function readForm(text: string, i: number, end: number, depth: number): Read<{ form: Node }> {
  const syntax = (message: string, at: number, hint?: string) => ({
    issue: issue('repr.tree.address.syntax', message, at, hint),
  });
  if (depth > MaxDepth)
    return {
      issue: issueOf(TreeIssues, 'repr.tree.depth_exceeded', `Nested deeper than ${MaxDepth} levels`, { offset: i }),
    };
  const open = i;
  i = skipSpace(text, i + 1, end);
  const headEnd = wordEnd(text, i, end);
  const head = text.slice(i, headEnd);
  if (!head || head.startsWith(':'))
    return { issue: issue('repr.tree.address.form', 'A query form opens with its kind', i, FormHint) };
  if (!KindPattern.test(head))
    return { issue: issueOf(TreeIssues, 'repr.tree.node.kind', 'Not a kind', { offset: i, hint: KindHint }) };
  i = headEnd;
  let params: Record<string, Value> | undefined;
  for (;;) {
    i = skipSpace(text, i, end);
    if (i >= end) return syntax('An unclosed query form', open);
    if (text[i] === ')') {
      i++;
      break;
    }
    if (text[i] !== ':')
      return { issue: issue('repr.tree.address.form', 'A query form holds keyword params only', i, FormHint) };
    const keyAt = i;
    let key: string;
    if (text[i + 1] === '"') {
      const q = readQuoted(text, i + 1, end);
      if (!q) return syntax('An unterminated or malformed keyword', i);
      [key, i] = q;
    } else {
      const e = wordEnd(text, i + 1, end);
      key = text.slice(i + 1, e);
      if (!KeywordPattern.test(key))
        return syntax('Not a keyword', i, 'a letter or _, then letters, digits, _ or -; or quote it, :"…"');
      i = e;
    }
    if (isReservedKey(key))
      return {
        issue: issueOf(TreeIssues, 'repr.tree.value.reserved_key', 'A reserved key', {
          offset: keyAt,
          hint: 'kind, __proto__, $…',
        }),
      };
    params ??= {};
    if (Object.hasOwn(params, key))
      return { issue: issue('repr.tree.address.form', 'A query form names a key twice', keyAt) };
    i = skipSpace(text, i, end);
    const v = readValue(text, i, end, depth);
    if ('issue' in v) return v;
    params[key] = v.value;
    i = v.next;
  }
  const form: Node = params ? { kind: head, params } : { kind: head };
  return { form, next: i };
}

/** A form value at `i`: an atom, a string, a list, a nested form or an address. */
function readValue(text: string, i: number, end: number, depth: number): Read<{ value: Value }> {
  if (i >= end || text[i] === ')' || text[i] === ']')
    return { issue: issue('repr.tree.address.form', 'A keyword has no value', i, FormHint) };
  const c = text[i];
  if (c === '(') {
    const r = readForm(text, i, end, depth + 1);
    return 'issue' in r ? r : { value: r.form, next: r.next };
  }
  if (c === '[') {
    if (depth + 1 > MaxDepth)
      return {
        issue: issueOf(TreeIssues, 'repr.tree.depth_exceeded', `Nested deeper than ${MaxDepth} levels`, { offset: i }),
      };
    const items: Value[] = [];
    let j = i + 1;
    for (;;) {
      j = skipSpace(text, j, end);
      if (j >= end) return { issue: issue('repr.tree.address.syntax', 'An unclosed list', i) };
      if (text[j] === ']') return { value: items, next: j + 1 };
      const v = readValue(text, j, end, depth + 1);
      if ('issue' in v) return v;
      items.push(v.value);
      j = v.next;
    }
  }
  if (c === '"') {
    const q = readQuoted(text, i, end);
    if (!q) return { issue: issue('repr.tree.address.syntax', 'An unterminated or malformed string', i) };
    return { value: q[0], next: q[1] };
  }
  if (c === '/') {
    const e = addressWordEnd(text, i, end);
    if (e < 0) return { issue: issue('repr.tree.address.syntax', 'An unclosed query form in an address', i) };
    const r = parseAddressIn(text, i, e, depth + 1);
    if (!r.pattern) return { issue: r.issues[0] };
    return { value: printAddress(r.pattern.segments, r.pattern.subtree), next: e };
  }
  if (c === ':') return { issue: issue('repr.tree.address.form', 'A keyword has no value', i, FormHint) };
  const e = wordEnd(text, i, end);
  const atom = e > i ? readWord(text.slice(i, e)) : undefined;
  if (atom === undefined) return { issue: issue('repr.tree.address.syntax', 'Not a word', i, 'quote it for a string') };
  return { value: atom, next: e };
}

/** Whether `text` is an address in its canonical spelling (the printers write such a string as a bare word). */
export function isCanonicalAddress(text: string): boolean {
  const { pattern } = parseAddress(text);
  return !!pattern && printAddress(pattern.segments, pattern.subtree) === text;
}

/**
 * The address an exact one-node address names, printed canonically; undefined when it is not one (not a string, not
 * an address, a query form, `/*`, or the root).
 */
export function exactAddress(text: unknown): string | undefined {
  if (typeof text !== 'string') return undefined;
  const p = parseAddress(text);
  if (p.issues.length || !p.pattern || p.pattern.subtree || !isExactAddress(p.pattern) || !p.pattern.segments.length)
    return undefined;
  return printAddress(p.pattern.segments as Address);
}

/** The node at an address, or undefined. */
export function resolveAddress(root: Node, address: Address): Node | undefined {
  let node: Node | undefined = root;
  for (const id of address) {
    node = node.children?.find((c) => c.id === id);
    if (!node) return undefined;
  }
  return node;
}

/** The address of `target` (by identity) in the tree at `root`; undefined if it is not there or has no address. */
export function addressOf(root: Node, target: Node): Address | undefined {
  if (root === target) return [];
  const path: string[] = [];
  const go = (node: Node): boolean => {
    for (const child of node.children ?? []) {
      if (child.id === undefined) continue;
      path.push(child.id);
      if (child === target || go(child)) return true;
      path.pop();
    }
    return false;
  };
  return go(root) ? path : undefined;
}

// ---- matching and selection ----

/**
 * Whether `node` fits the query form (meso*'s child selection): its kind, and for each of the form's params a value of
 * the node's own that is the same atom (`Object.is`: `0` is not `-0`, `nan` is `nan`, `1` is not `1n`, absent is not
 * `nil`), a list param holding that atom, an equal list, or a node fitting the nested form. Params the form does not
 * name, the node's id and its children are not looked at; defaults are not applied.
 */
export function matchesForm(form: Node, node: Node): boolean {
  if (node.kind !== form.kind) return false;
  const params = form.params;
  if (!params) return true;
  for (const key in params) {
    if (!Object.hasOwn(params, key)) continue;
    if (!matchesValue(params[key], node.params ? own(node.params, key) : undefined)) return false;
  }
  return true;
}

function matchesValue(want: Value, have: Value | undefined): boolean {
  if (isNode(want)) return isNode(have) && matchesForm(want, have);
  if (Array.isArray(want)) return valuesEqual(want, have);
  // an atom also matches a list param holding it (one scan of the node's own list, which its schema bounds)
  if (Array.isArray(have)) {
    for (const h of have) if (Object.is(want, h)) return true;
    return false;
  }
  return Object.is(want, have);
}

/**
 * Whether `address` is named by `pattern`: equal to it, or below it when the pattern names a subtree. A pattern test,
 * not existence: id segments compare ids, a form is tested (by child selection) on the node the address names at
 * that depth in `root`. A pattern with forms needs `root`; without it, that is a programmer error (it throws). A
 * pattern given as an object is checked as `parseAddress` checks text; one that is not a pattern matches nothing.
 */
export function matchesAddress(pattern: AddressPattern, address: Address, root?: Node): boolean {
  if (patternIssue(pattern)) return false;
  if (!root && !isExactAddress(pattern)) throw new Error('matchesAddress: a pattern with query forms needs the tree');
  const p = pattern.segments;
  if (address.length < p.length || (!pattern.subtree && address.length !== p.length)) return false;
  let node: Node | undefined = root;
  for (let i = 0; i < p.length; i++) {
    const s = p[i];
    if (typeof s === 'string') {
      if (s !== address[i]) return false;
      if (node) node = node.children?.find((c) => c.id === s);
    } else {
      node = node?.children?.find((c) => c.id === address[i]);
      if (!node || !matchesForm(s, node)) return false;
    }
  }
  return true;
}

/** A bound: a non-negative integer, or `Infinity` said explicitly (never NaN, a fraction or a negative). */
const isCount = (n: unknown) => n === Number.POSITIVE_INFINITY || (Number.isInteger(n) && (n as number) >= 0);
const refuse = (code: Code, message: string, hint?: string): Selection => ({
  matches: [],
  issues: [issue(code, message, undefined, hint)],
});

/**
 * The first key of a query form, or of a form nested in its values, that its kind's params lack; kinds not in `params`
 * pass.
 */
export function unknownFormKey(form: Node, params: ReadonlyMap<string, ReadonlySet<string>>): string | undefined {
  const values = form.params;
  if (!values) return undefined;
  const known = params.get(form.kind);
  for (const key of Object.keys(values)) {
    if (known && !known.has(key)) return key;
    for (const nested of nestedForms(values[key])) {
      const k = unknownFormKey(nested, params);
      if (k !== undefined) return k;
    }
  }
  return undefined;
}

function nestedForms(v: Value): Node[] {
  if (isNode(v)) return [v];
  if (Array.isArray(v)) return v.flatMap(nestedForms);
  return [];
}

/** Whether the form's kind and every nested form's kind are among `kinds`. */
function knownKinds(form: Node, kinds: ReadonlySet<string>): boolean {
  if (!kinds.has(form.kind)) return false;
  const params = form.params;
  if (params)
    for (const key of Object.keys(params))
      for (const n of nestedForms(params[key])) if (!knownKinds(n, kinds)) return false;
  return true;
}

/** A selected node and its address (ids only: every selected node is addressable). */
export interface Selected {
  readonly node: Node;
  readonly address: Address;
}

/**
 * What an address, or its prefix, names: nodes (exact addresses), or what a form resolver made of them (`type` is
 * the resolver's: a location, elements; `value` is opaque to the tree).
 */
export type Named =
  | { readonly type: 'nodes'; readonly matches: readonly Selected[] }
  | { readonly type: string; readonly value: unknown };

export interface FormResolverInput {
  /** The form, as read. */
  readonly form: Node;
  /** What the address before the form names (the root, as nodes, for a first segment). */
  readonly prefix: Named;
  readonly root: Node;
  /** The caller's limit: a result of more nodes or elements is the resolver's issue, never a truncation. */
  readonly limit: number;
}

/**
 * What a form of one kind means (locations-and-queries.md): it receives what the prefix names and returns one of its
 * result types, or issues. A result of nodes lets ids, child selection and `/*` follow it.
 */
export type FormResolver = (input: FormResolverInput) => {
  readonly result?: Named;
  readonly issues?: readonly Issue[];
};

export interface SelectOptions {
  /**
   * At most this many matches, a non-negative integer or an explicit `Infinity` (else `repr.tree.address.option`);
   * more is `repr.tree.address.limit_exceeded` and no matches.
   */
  readonly limit: number;
  /**
   * How many levels below the root the walk may go: a count (else `repr.tree.address.option`), at most and by default
   * MaxDepth, so the recursive walk never overflows the stack; deeper is `repr.tree.address.depth_exceeded`.
   */
  readonly maxDepth?: number;
  /** The kinds a form, or a form nested in its values, may name; another is `repr.tree.address.unknown_kind`. */
  readonly kinds?: ReadonlySet<string>;
  /**
   * The params each kind has: a form keyed by a param its kind lacks is `repr.tree.address.unknown_param` (it could
   * match nothing). A kind absent is not checked. Forms with a resolver are the resolver's to check.
   */
  readonly params?: ReadonlyMap<string, ReadonlySet<string>>;
  /**
   * The form-resolver table, by form kind; a form whose kind has none is meso*'s child selection (`matchesForm`), which
   * applies to nodes only.
   */
  readonly resolvers?: Readonly<Record<string, FormResolver>>;
  /** Cancels the walk: read before it and every `CancelEvery` nodes it visits (and between segments). */
  readonly signal?: AbortSignal;
  /** A `performance.now()` time the walk stops at, read as the signal is. */
  readonly deadline?: number;
}

/** How many nodes a walk visits between two reads of its signal. */
export const CancelEvery = 1024;

export interface Selection {
  /** The nodes named, in document order, each once; empty when there is an issue or the result is not nodes. */
  readonly matches: readonly Selected[];
  /** What a form resolver made, when the address names something other than nodes. */
  readonly result?: { readonly type: string; readonly value: unknown };
  readonly issues: readonly Issue[];
}

/**
 * What a pattern names, with the addresses of the nodes: an id segment steps to the child with that id, a form to
 * every addressable child fitting it (or to its resolver's result), and `subtree` adds every addressable node below
 * each. All ids is a walk of O(depth × fan-out). The walk keeps one path stack; an address is copied per match only.
 */
export function selectAddress(root: Node, address: AddressPattern | string, options: SelectOptions): Selection {
  let pattern: AddressPattern;
  if (typeof address === 'string') {
    const parsed = parseAddress(address);
    if (!parsed.pattern) return { matches: [], issues: parsed.issues };
    pattern = parsed.pattern;
  } else {
    pattern = address;
    const wrong = patternIssue(pattern);
    if (wrong) return { matches: [], issues: [wrong] };
  }
  const { limit, kinds, resolvers } = options;
  if (!isCount(limit)) return refuse('repr.tree.address.option', 'The limit is not a non-negative integer or Infinity');
  if (options.maxDepth !== undefined && !isCount(options.maxDepth))
    return refuse('repr.tree.address.option', 'The maxDepth is not a non-negative integer or Infinity');
  const maxDepth = Math.min(options.maxDepth ?? MaxDepth, MaxDepth);
  const segments = pattern.segments;
  const resolved = (s: Segment) => typeof s !== 'string' && !!resolvers && own(resolvers, s.kind) !== undefined;
  if (kinds) {
    for (let i = 0; i < segments.length; i++) {
      const s = segments[i];
      if (typeof s === 'string' || resolved(s) || knownKinds(s, kinds)) continue;
      const known = [...kinds].slice(0, 20).join(', ') + (kinds.size > 20 ? ', …' : '');
      const message = `A query form (segment ${i + 1}) names a kind the language does not have`;
      return { matches: [], issues: [issue('repr.tree.address.unknown_kind', message, undefined, `kinds: ${known}`)] };
    }
  }
  if (options.params) {
    for (let i = 0; i < segments.length; i++) {
      const s = segments[i];
      if (typeof s === 'string' || resolved(s)) continue;
      if (unknownFormKey(s, options.params) === undefined) continue;
      // the key is the spec's, not quoted: the segment says where
      const message = `A query form (segment ${i + 1}) names a param its kind does not have`;
      return { matches: [], issues: [issue('repr.tree.address.unknown_param', message)] };
    }
  }
  const { signal, deadline } = options;
  const timed = deadline !== undefined;
  /** The issue when the signal aborted or the deadline passed (each read once), else undefined. */
  const over = (): Issue | undefined =>
    signal?.aborted
      ? issue('repr.tree.address.cancelled', 'The selection was cancelled')
      : timed && performance.now() > deadline
        ? issue('repr.tree.address.deadline_exceeded', 'The selection ran past its deadline', undefined, 'narrow it')
        : undefined;
  const before = over();
  if (before) return { matches: [], issues: [before] };
  if (segments.some(resolved)) return staged(root, pattern, { limit, maxDepth, resolvers: resolvers!, over });

  const matches: Selected[] = [];
  const path: string[] = [];
  let stop: Issue | undefined;
  let visits = 0;
  /** Counts a visited node; every CancelEvery reads the signal and the clock (no allocation). */
  const visit = () => {
    if ((!signal && !timed) || ++visits < CancelEvery) return true;
    visits = 0;
    stop = over();
    return !stop;
  };
  const emit = (node: Node) => {
    if (matches.length >= limit) {
      stop = limitIssue(limit);
      return;
    }
    matches.push({ node, address: path.slice() });
  };
  const deeper = () => {
    if (path.length < maxDepth) return true;
    stop = issue('repr.tree.address.depth_exceeded', `The walk goes deeper than ${maxDepth} levels`);
    return false;
  };
  const below = (node: Node) => {
    emit(node);
    const children = node.children;
    if (stop || !children) return;
    for (let k = 0; k < children.length; k++) {
      const id = children[k].id;
      if (id === undefined) continue;
      if (!visit() || !deeper()) return;
      path.push(id);
      below(children[k]);
      path.pop();
      if (stop) return;
    }
  };
  const step = (node: Node, i: number) => {
    if (i === segments.length) {
      if (pattern.subtree) below(node);
      else emit(node);
      return;
    }
    const children = node.children;
    if (!children) return;
    const s = segments[i];
    for (let k = 0; k < children.length; k++) {
      const child = children[k];
      const id = child.id;
      if (id === undefined) continue;
      if (!visit()) return;
      if (typeof s === 'string' ? id !== s : !matchesForm(s, child)) continue;
      if (!deeper()) return;
      path.push(id);
      step(child, i + 1);
      path.pop();
      if (stop || typeof s === 'string') return; // ids are unique among siblings
    }
  };
  step(root, 0);
  return stop ? { matches: [], issues: [stop] } : { matches, issues: [] };
}

const limitIssue = (limit: number) =>
  issue(
    'repr.tree.address.limit_exceeded',
    `More than ${limit} matches`,
    undefined,
    'narrow the address, or raise the limit',
  );

/**
 * A pattern with resolved forms, segment by segment: what each prefix names is handed to the next step, nodes or a
 * resolver's result. Bounded by the limit per step and by the signal and the deadline between steps.
 */
function staged(
  root: Node,
  pattern: AddressPattern,
  o: {
    limit: number;
    maxDepth: number;
    resolvers: Readonly<Record<string, FormResolver>>;
    over: () => Issue | undefined;
  },
): Selection {
  const fail = (message: string): Selection => ({
    matches: [],
    issues: [issue('repr.tree.address.resolver', message)],
  });
  let named: Named = { type: 'nodes', matches: [{ node: root, address: [] }] };
  const nodesOf = (n: Named) => (n.type === 'nodes' ? (n as { matches: readonly Selected[] }).matches : undefined);
  for (let i = 0; i < pattern.segments.length; i++) {
    const s = pattern.segments[i];
    const late = o.over();
    if (late) return { matches: [], issues: [late] };
    const resolver = typeof s === 'string' ? undefined : own(o.resolvers, s.kind);
    if (resolver) {
      const r = resolver({ form: s as Node, prefix: named, root, limit: o.limit });
      if (r.issues?.length) return { matches: [], issues: r.issues };
      if (!r.result) return fail(`The resolver of segment ${i + 1} gave no result`);
      named = r.result;
      continue;
    }
    const prefix = nodesOf(named);
    if (!prefix) return fail(`Segment ${i + 1} applies to nodes, and its prefix names a ${named.type}`);
    const next: Selected[] = [];
    for (const m of prefix) {
      if (m.address.length >= o.maxDepth)
        return {
          matches: [],
          issues: [issue('repr.tree.address.depth_exceeded', `The walk goes deeper than ${o.maxDepth} levels`)],
        };
      for (const child of m.node.children ?? []) {
        if (child.id === undefined) continue;
        if (typeof s === 'string' ? child.id !== s : !matchesForm(s, child)) continue;
        if (next.length >= o.limit) return { matches: [], issues: [limitIssue(o.limit)] };
        next.push({ node: child, address: [...m.address, child.id] });
        if (typeof s === 'string') break;
      }
    }
    named = { type: 'nodes', matches: next };
  }
  const nodes = nodesOf(named);
  if (!nodes) {
    if (pattern.subtree) return fail(`/* applies to nodes, and the address names a ${named.type}`);
    return { matches: [], result: named as { type: string; value: unknown }, issues: [] };
  }
  if (!pattern.subtree) {
    if (nodes.length > o.limit) return { matches: [], issues: [limitIssue(o.limit)] };
    return { matches: nodes, issues: [] };
  }
  // every addressable node below each, each once (a resolver's nodes may hold one another)
  const out: Selected[] = [];
  const seen = new Set<Node>();
  const below = (m: Selected): Issue | undefined => {
    if (seen.has(m.node)) return undefined;
    seen.add(m.node);
    if (out.length >= o.limit) return limitIssue(o.limit);
    out.push(m);
    for (const child of m.node.children ?? []) {
      if (child.id === undefined) continue;
      if (m.address.length >= o.maxDepth)
        return issue('repr.tree.address.depth_exceeded', `The walk goes deeper than ${o.maxDepth} levels`);
      const stop = below({ node: child, address: [...m.address, child.id] });
      if (stop) return stop;
    }
    return undefined;
  };
  for (const m of nodes) {
    const stop = below(m);
    if (stop) return { matches: [], issues: [stop] };
  }
  return { matches: out, issues: [] };
}
