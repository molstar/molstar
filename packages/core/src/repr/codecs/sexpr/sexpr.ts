/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The S-expression printing of the tree and its reader (_spec/README.md "Codecs"; format.md): for people, docs, tests
// and agents; never a file type, never stored. The reader lowers the same spelling to the same tree: a column in a
// declared column param (the printer is given the language) prints as a list, any other tagged,
// `(:$column int32 :values [1 2])`; past its head a column is elided in a brace form the reader refuses,
// `{float32 × 3,000,000: 0.12 0.15 …}`.
//
//   node    (kind [id] :param value … :$custom (…) child …)   a word head; an optional id; keyword params; custom
//                                                 data under the reserved keyword `:$custom`; children are forms
//   value   atom | [value …] | (:key value …) | node      a list, a record (a form opening with a keyword or empty), a node
//   atoms   tree/atoms.ts: "string" or a bare word (94A, #ff0000), numbers as JavaScript reads them, nan, +inf, -inf,
//           123n, true, false, nil
//   address /data/p53/(range :from 94 :to 292)   a word opening with `/` is an address, read through its query forms
//                                                 (tree/address.ts): a string value, in its canonical spelling
//   ;       a comment, to the end of the line
//
// A printing, not the canonical text (that is the hash's). Atoms are tree/atoms.ts's (a word that is no JavaScript
// number is a symbol, so no number issue and no units table in the reader), every word opening with `/` is read as an
// address, `:$custom` holds custom data, columns are elided, and definitions' declared params wait for definitions.

import { CodecIssues, type Codec, type DecodeOptions, type Decoded } from '../../codec/codec.js';
import { type ColumnPlaces, declaredColumnPlaces, type Language } from '../../schema/index.js';
import {
  addressWordEnd,
  assertKeys,
  assertPrintable,
  type Atom,
  CategoryColumn,
  type Column,
  type ColumnType,
  ColumnTypes,
  categoryFromCodes,
  columnType,
  Delimiter,
  IdHint,
  IdPattern,
  type Issue,
  isColumn,
  isNode,
  isReservedKey,
  issueOf,
  KeywordPattern,
  KindHint,
  KindPattern,
  MaxDepth,
  type Node,
  parseAddressIn,
  printAddress,
  printAddressAtom,
  printKey,
  readQuoted,
  readWord,
  TreeIssues,
  type TreePath,
  toColumn,
  type Value,
  type ValueRecord,
} from '../../tree/index.js';

/** Reserved: never in a word, so a later spelling may use them. */
const RESERVED = /[{},]/;
/** The keyword custom data is spelled under: `$` keeps it apart from every param. */
const CustomKey = '$custom';
/** The keyword a tagged column is spelled with: `(:$column int32 :values [1 2])`. */
const ColumnKey = '$column';
/** An elided column's opening, a brace form no value reads as (braces are reserved): `{float32 × `. */
const ELIDED = new RegExp(`^\\{(${ColumnTypes.join('|')}) ×`);

// ---- reading ----

type SNode =
  | { readonly t: 'list' | 'vector'; readonly items: SNode[]; readonly start: number }
  | { readonly t: 'atom'; readonly value: Atom; readonly word: boolean; readonly start: number }
  | { readonly t: 'keyword'; readonly name: string; readonly start: number };

type Code = keyof typeof CodecIssues.codes | keyof typeof TreeIssues.codes;
const issue = (code: Code, message: string, offset: number, path?: TreePath, hint?: string): Issue =>
  code.startsWith('repr.codec.')
    ? issueOf(CodecIssues, code as keyof typeof CodecIssues.codes, message, { offset, path, hint })
    : issueOf(TreeIssues, code as keyof typeof TreeIssues.codes, message, { offset, path, hint });

/**
 * The forms of `text` between `start` and `end`, or its first syntax issue. Messages never quote the text: the offset
 * says where. A leading byte order mark is skipped; offsets count it.
 */
function read(text: string, start = 0, end = text.length): { forms: SNode[]; error?: Issue } {
  const stack: { t: 'list' | 'vector'; items: SNode[]; start: number }[] = [];
  const top: SNode[] = [];
  const push = (s: SNode) => (stack.length ? stack[stack.length - 1].items : top).push(s);
  const fail = (message: string, offset: number, hint?: string, code: Code = 'repr.codec.sexpr.syntax') => ({
    forms: top,
    error: issue(code, message, offset, undefined, hint),
  });
  let i = start === 0 && text.charCodeAt(0) === 0xfeff ? 1 : start;
  while (i < end) {
    const c = text[i];
    if (c === ' ' || c === '\t' || c === '\r' || c === '\n') i++;
    else if (c === ';') while (i < end && text[i] !== '\n') i++;
    else if (c === '(' || c === '[') {
      if (stack.length === MaxDepth)
        return fail(`Nested deeper than ${MaxDepth} levels`, i, undefined, 'repr.tree.depth_exceeded');
      stack.push({ t: c === '(' ? 'list' : 'vector', items: [], start: i++ });
    } else if (c === ')' || c === ']') {
      const open = stack.pop();
      if (!open) return fail(`Unexpected ${c}`, i);
      if ((open.t === 'list') !== (c === ')')) return fail(`${c} closes the ${open.t} opened at ${open.start}`, i);
      push(open);
      i++;
    } else if (c === '{' && ELIDED.test(text.slice(i, i + 16)))
      return fail('An elided column', i, CodecIssues.codes['repr.codec.sexpr.elided'].hint, 'repr.codec.sexpr.elided');
    else if (RESERVED.test(c)) return fail(`${c} is reserved`, i, 'quote a string that holds it');
    else if (c === '"') {
      const s = readQuoted(text, i, end);
      if (!s) return fail('An unterminated or malformed string', i);
      push({ t: 'atom', value: s[0], word: false, start: i });
      i = s[1];
    } else if (c === ':' && text[i + 1] === '"') {
      const s = readQuoted(text, i + 1, end);
      if (!s) return fail('An unterminated or malformed keyword', i);
      push({ t: 'keyword', name: s[0], start: i });
      i = s[1];
    } else if (c === '/') {
      // an address (tree/address.ts): a word through balanced query forms, `/data/p53/(range :from 94 :to 292)`
      const close = addressWordEnd(text, i, end);
      if (close < 0) return fail('An unclosed or malformed query form in an address', i);
      const { pattern, issues } = parseAddressIn(text, i, close);
      if (!pattern) return { forms: top, error: issues[0] };
      push({ t: 'atom', value: printAddress(pattern.segments, pattern.subtree), word: true, start: i });
      i = close;
    } else {
      const from = i;
      while (i < end && !Delimiter.test(text[i])) i++;
      const raw = text.slice(from, i);
      if (c === ':') {
        const name = raw.slice(1);
        if (name !== CustomKey && name !== ColumnKey && !KeywordPattern.test(name))
          return fail('Not a keyword', from, 'a letter or _, then letters, digits, _ or -; or quote it, :"…"');
        push({ t: 'keyword', name, start: from });
        continue;
      }
      const atom = readWord(raw);
      if (atom === undefined) return fail('Not a word', from, 'quote it for a string');
      push({ t: 'atom', value: atom, word: typeof atom === 'string', start: from });
    }
  }
  if (stack.length) return fail(`An unclosed ${stack[stack.length - 1].t}`, stack[stack.length - 1].start);
  return { forms: top };
}

/**
 * A form as the reader reads it, before it is a tree: a list `( … )`, a vector `[ … ]`, an atom (`word`: written bare)
 * or a keyword, each with its offset. For spellings over the tree's own (a patch's positional arguments): they read
 * their forms here and their values through `formToValue`, so atoms, nodes and issues are the tree's.
 */
export type Form = SNode;

/** How a reader reads forms: which kinds are data forms. */
export interface ReadOptions {
  /** Kinds whose positional items are values (`args`), not an id and children; none by default. */
  readonly dataForms?: ReadonlySet<string>;
}

const NoDataForms: ReadonlySet<string> = new Set();

/** Every top-level form of `text`, or its first syntax issue (an offset into the text). */
export function readForms(text: string): { readonly forms: readonly Form[]; readonly issues: readonly Issue[] } {
  const { forms, error } = read(text);
  return error ? { forms: [], issues: [error] } : { forms, issues: [] };
}

/** A form read as a tree value (an atom, a list, a record or a node), its issues at `path` (and their offsets). */
export function formToValue(
  form: Form,
  path: TreePath = [],
  options: ReadOptions = {},
): { readonly value: Value; readonly issues: Issue[] } {
  const issues: Issue[] = [];
  const value = toValue(form, path, issues, options.dataForms ?? NoDataForms);
  return { value, issues };
}

/** The tree of a text holding exactly one node, strictly; every shape issue with its offset and path. */
export function parseSexpr(text: string, options: ReadOptions = {}): Decoded {
  const { forms, error } = read(text);
  if (error) return { issues: [error] };
  const issues: Issue[] = [];
  if (forms.length !== 1 || forms[0].t !== 'list') {
    const at = forms.length > 1 ? forms[1].start : (forms[0]?.start ?? 0);
    return { issues: [issue('repr.codec.sexpr.shape', 'The text must hold exactly one node, a form', at)] };
  }
  const node = toNode(forms[0], [], issues, options.dataForms ?? NoDataForms);
  return issues.length ? { issues } : { node, issues };
}

type DataForms = ReadonlySet<string>;

function toNode(s: Extract<SNode, { items: SNode[] }>, path: TreePath, issues: Issue[], data: DataForms): Node {
  const [head, ...rest] = s.items;
  if (head?.t !== 'atom' || !head.word) {
    issues.push(issue('repr.codec.sexpr.shape', 'A node opens with its kind, a word', head?.start ?? s.start, path));
    return { kind: '' };
  }
  if (!KindPattern.test(head.value as string))
    issues.push(issue('repr.tree.node.kind', 'Not a kind', head.start, [...path, 'kind'], KindHint));
  const node: {
    kind: string;
    id?: string;
    params?: Record<string, Value>;
    custom?: Record<string, Value>;
    children?: Node[];
    args?: Value[];
  } = { kind: head.value as string };
  const custom = (value: SNode, at: number) => {
    const v = toValue(value, [...path, 'custom'], issues, data);
    if (typeof v === 'object' && v !== null && !Array.isArray(v) && !isNode(v) && !isColumn(v))
      node.custom = v as Record<string, Value>;
    else
      issues.push(issue('repr.codec.sexpr.shape', 'Custom data is a record, (:key value …)', at, [...path, 'custom']));
  };
  if (data.has(node.kind)) {
    const params = pairs(
      rest,
      [...path, 'params'],
      issues,
      data,
      (item) => {
        node.args ??= [];
        node.args.push(toValue(item, [...path, 'args', node.args.length], issues, data));
      },
      custom,
    );
    if (params) node.params = params;
    return node;
  }
  let k = 0;
  if (rest[0]?.t === 'atom') {
    const id = rest[0];
    if (typeof id.value !== 'string' || !IdPattern.test(id.value))
      issues.push(issue('repr.tree.node.id', 'Not an id', id.start, [...path, 'id'], IdHint));
    else node.id = id.value;
    k = 1;
  }
  const params = pairs(
    rest.slice(k),
    [...path, 'params'],
    issues,
    data,
    (child, i) => {
      if (child.t === 'list' && child.items[0]?.t === 'atom' && child.items[0].word) {
        node.children ??= [];
        node.children.push(toNode(child, [...path, 'children', node.children.length], issues, data));
      } else
        issues.push(
          issue('repr.codec.sexpr.shape', 'After the id, a node holds keyword params and child nodes', i, path),
        );
    },
    custom,
  );
  if (params) node.params = params;
  return node;
}

/** Keyword pairs; anything else goes to `other`; `:$custom` to `custom` where a node takes it. */
function pairs(
  items: readonly SNode[],
  path: TreePath,
  issues: Issue[],
  data: DataForms,
  other: (s: SNode, offset: number) => void,
  custom?: (value: SNode, offset: number) => void,
): Record<string, Value> | undefined {
  let out: Record<string, Value> | undefined;
  let customSeen = false;
  for (let i = 0; i < items.length; i++) {
    const s = items[i];
    if (s.t !== 'keyword') {
      other(s, s.start);
      continue;
    }
    const value = items[i + 1];
    if (!value || value.t === 'keyword') {
      issues.push(issue('repr.codec.sexpr.shape', 'A keyword has no value', s.start, [...path, s.name]));
      continue;
    }
    i++;
    if (s.name === CustomKey && custom) {
      if (customSeen) issues.push(issue('repr.codec.sexpr.duplicate_key', 'A keyword given twice', s.start, path));
      else custom(value, s.start);
      customSeen = true;
      continue;
    }
    out ??= {};
    if (isReservedKey(s.name))
      issues.push(
        issue('repr.tree.value.reserved_key', 'A reserved key', s.start, [...path, s.name], 'kind, __proto__, $…'),
      );
    else if (Object.hasOwn(out, s.name))
      issues.push(issue('repr.codec.sexpr.duplicate_key', 'A keyword given twice', s.start, [...path, s.name]));
    else out[s.name] = toValue(value, [...path, s.name], issues, data);
  }
  return out;
}

function toValue(s: SNode, path: TreePath, issues: Issue[], data: DataForms): Value {
  if (s.t === 'atom') return s.value;
  if (s.t === 'vector') return s.items.map((x, i) => toValue(x, [...path, i], issues, data));
  if (s.t === 'list') {
    const head = s.items[0];
    if (head?.t === 'atom' && head.word) return toNode(s, path, issues, data);
    if (head?.t === 'keyword' && head.name === ColumnKey) return taggedColumn(s.items, s.start, path, issues, data);
    if (!head || head.t === 'keyword')
      return (
        pairs(s.items, path, issues, data, (_, at) =>
          issues.push(issue('repr.codec.sexpr.shape', 'A record holds keyword pairs only', at, path)),
        ) ?? {}
      );
  }
  issues.push(
    issue(
      'repr.codec.sexpr.shape',
      'Not a value: a form opens with a kind word, a record with a keyword',
      s.start,
      path,
    ),
  );
  return null;
}

/**
 * A tagged column, `(:$column int32 :values [1 2])` or `(:$column category :dict [a b] :codes [0 1 -1])`, typed here;
 * one that does not read is `repr.tree.value.tag`.
 */
function taggedColumn(items: readonly SNode[], at: number, path: TreePath, issues: Issue[], data: DataForms): Value {
  const fields: Record<string, Value> = {};
  for (let i = 0; i + 1 < items.length; i += 2) {
    const k = items[i];
    if (k.t !== 'keyword') break;
    fields[k.name] = toValue(items[i + 1], path, issues, data);
  }
  const type = fields[ColumnKey];
  const keys = Object.keys(fields).sort().join(',');
  let column: Column | undefined;
  if (items.length % 2 === 0 && type === 'category' && keys === '$column,codes,dict') {
    const { dict, codes } = fields;
    if (
      Array.isArray(dict) &&
      dict.every((d) => typeof d === 'string') &&
      Array.isArray(codes) &&
      codes.every((c) => Number.isInteger(c) && (c as number) >= -1 && (c as number) < dict.length)
    )
      column = categoryFromCodes(dict as string[], codes as number[]);
  } else if (items.length % 2 === 0 && typeof type === 'string' && keys === '$column,values')
    column = ColumnTypes.includes(type as ColumnType) ? toColumn(type as ColumnType, fields.values) : undefined;
  if (column) return column;
  issues.push(issue('repr.tree.value.tag', 'A malformed tagged column', at, path, '(:$column int32 :values [1 2])'));
  return null;
}

// ---- printing ----

const WIDTH = 80;

/** A layout spec: an atom's text, or a form with its head and its parts (keyword or positional). */
type Doc =
  | string
  | { readonly open: '(' | '['; readonly head?: string; readonly parts: readonly [string | null, Doc][] };

export interface PrintOptions {
  /** A column longer than this prints elided, its first `head` values shown (default 8). */
  readonly head?: number;
  /**
   * The language whose declared column params print as plain lists; every other column is tagged
   * (`(:$column int32 :values [1 2])`). Without one, every column is tagged.
   */
  readonly language?: Language;
  /** The languages its slots name. */
  readonly languages?: (id: string) => Language | undefined;
}

const close = (o: '(' | '[') => (o === '(' ? ')' : ']');

/** A record's or params' entries, `undefined` ones absent (as in JSON and the hash). */
const entries = (r: ValueRecord) => Object.keys(r).filter((k) => r[k] !== undefined);

/** `3000000` as `3,000,000`. */
const grouped = (n: number) => String(n).replace(/\B(?=(\d{3})+(?!\d))/g, ',');

function printer(head: number, places?: ColumnPlaces) {
  const declared = (container: object, key: string | number) => places?.get(container)?.has(key) ?? false;
  function nodeDoc(node: Node): Doc {
    assertPrintable(node);
    const parts: [string | null, Doc][] = [];
    if (node.id !== undefined) parts.push([null, printAddressAtom(node.id)]);
    // a data form's leading values that are not forms go first, as an id does; its forms after the params, as children
    // do; the args keep their order either way
    const args = node.args ?? [];
    let lead = 0;
    while (lead < args.length && !isNode(args[lead])) parts.push([null, valueDoc(args[lead], declared(args, lead++))]);
    for (const key of entries(node.params ?? {}))
      parts.push([key, valueDoc(node.params![key], declared(node.params!, key))]);
    if (node.custom && entries(node.custom).length) parts.push([CustomKey, valueDoc(node.custom)]);
    for (const child of node.children ?? []) parts.push([null, nodeDoc(child)]);
    for (let i = lead; i < args.length; i++) parts.push([null, valueDoc(args[i], declared(args, i))]);
    return { open: '(', head: node.kind, parts };
  }

  const list = (items: readonly string[]): Doc => ({ open: '[', parts: items.map((x) => [null, x] as [null, Doc]) });

  /** A declared column a plain list, any other tagged; past the head, elided in a brace form the reader refuses. */
  function columnDoc(column: Column, plain: boolean): Doc {
    const type = columnType(column)!;
    const n = Math.min(column.length, head);
    const items: string[] = new Array(n);
    for (let i = 0; i < n; i++)
      items[i] = printAddressAtom(column instanceof CategoryColumn ? column.get(i) : column[i]);
    if (column.length > head) return `{${type} × ${grouped(column.length)}:${items.map((x) => ` ${x}`).join('')} …}`;
    if (plain) return list(items);
    if (column instanceof CategoryColumn)
      return {
        open: '(',
        parts: [
          [ColumnKey, 'category'],
          ['dict', list(column.dict.map(printAddressAtom))],
          ['codes', list(Array.from(column.codes, String))],
        ],
      };
    return {
      open: '(',
      parts: [
        [ColumnKey, type],
        ['values', list(items)],
      ],
    };
  }

  function valueDoc(v: Value, plain = false): Doc {
    if (v === null || typeof v !== 'object') return printAddressAtom(v);
    if (isColumn(v)) return columnDoc(v, plain);
    if (Array.isArray(v))
      return { open: '[', parts: v.map((x, i) => [null, valueDoc(x, declared(v, i))] as [null, Doc]) };
    if (isNode(v)) return nodeDoc(v);
    const r = v as ValueRecord;
    assertKeys(r);
    return { open: '(', parts: entries(r).map((k) => [k, valueDoc(r[k], declared(r, k))] as [string, Doc]) };
  }
  return { nodeDoc, valueDoc };
}

const keyText = (k: string) => (k === CustomKey || k === ColumnKey ? `:${k}` : printKey(k));

function flat(d: Doc): string {
  if (typeof d === 'string') return d;
  const parts = d.parts.map(([k, v]) => (k === null ? flat(v) : `${keyText(k)} ${flat(v)}`));
  return `${d.open}${[...(d.head !== undefined ? [d.head] : []), ...parts].join(' ')}${close(d.open)}`;
}

/**
 * One line when it fits. Else the head, the id and, when they all fit there, the keyword params on the first line;
 * then one argument a line, two spaces in from the form's column; a keyword whose value does not fit puts the value on
 * the next line, two spaces in from the keyword; a list of atoms fills its lines.
 */
function layout(d: Doc, column: number): string {
  const line = flat(d);
  if (typeof d === 'string' || column + line.length <= WIDTH) return line;
  if (d.open === '[' && d.parts.every(([, v]) => typeof v === 'string')) return fill(d.parts, column);
  let out = d.open + (d.head ?? '');
  let k = 0;
  if (d.head !== undefined) {
    if (d.parts[0]?.[0] === null && typeof d.parts[0][1] === 'string') out += ` ${d.parts[k++][1]}`;
    let end = k;
    while (end < d.parts.length && d.parts[end][0] !== null) end++;
    const params = d.parts.slice(k, end).map(([key, v]) => `${keyText(key!)} ${flat(v)}`);
    if (params.length && column + out.length + 1 + params.join(' ').length <= WIDTH) {
      out += ` ${params.join(' ')}`;
      k = end;
    }
  }
  const inner = column + (d.head === undefined ? 1 : 2);
  for (; k < d.parts.length; k++) {
    const [key, v] = d.parts[k];
    const first = out.length === 1;
    const lead = key === null ? '' : `${keyText(key)} `;
    const text = flat(v);
    let body: string;
    if (typeof v === 'string' || inner + lead.length + text.length <= WIDTH) body = lead + text;
    else if (key === null) body = layout(v, inner);
    else body = `${lead.trimEnd()}\n${' '.repeat(inner + 2)}${layout(v, inner + 2)}`;
    out += first ? body : `\n${' '.repeat(inner)}${body}`;
  }
  return out + close(d.open);
}

function fill(parts: readonly [string | null, Doc][], column: number): string {
  const lines: string[] = [];
  let cur = '[';
  let at = column; // the column `cur` starts at: the list's own on its first line, 0 after (`cur` holds its indent)
  for (const [, v] of parts) {
    const t = v as string;
    if (cur !== '[' && at + cur.length + 1 + t.length > WIDTH) {
      lines.push(cur);
      cur = ' '.repeat(column + 1) + t;
      at = 0;
    } else cur += (cur === '[' ? '' : ' ') + t;
  }
  return [...lines, `${cur}]`].join('\n');
}

function printerFor(node: Node, options: PrintOptions) {
  const places = options.language
    ? declaredColumnPlaces(options.language, node, { languages: options.languages })
    : undefined;
  return printer(options.head ?? 8, places);
}

/** The tree's printing: `parseSexpr(printSexpr(x))` is `x` once normalized, unless a column was elided. */
export function printSexpr(node: Node, options: PrintOptions = {}): string {
  return layout(printerFor(node, options).nodeDoc(node), 0);
}

/**
 * A form of a spelling over the tree's (a patch's), laid out as the tree's forms are, from `column`: its head (a word),
 * then its arguments, each a positional value (`key` null) or a keyword's. The head is the caller's to check.
 */
export function printForm(
  head: string,
  args: readonly (readonly [string | null, Value])[],
  column = 0,
  options: PrintOptions = {},
): string {
  const { valueDoc } = printer(options.head ?? 8);
  return layout({ open: '(', head, parts: args.map(([k, v]) => [k, valueDoc(v)] as [string | null, Doc]) }, column);
}

/** The tree's printing on one line, however long. */
export function printSexprLine(node: Node, options: PrintOptions = {}): string {
  return flat(printerFor(node, options).nodeDoc(node));
}

/** The S-expression codec: a printing for people that elides long columns (`head`), whose reader refuses them. */
export function sexprCodec(options: PrintOptions = {}): Codec<string> {
  return {
    id: 'sexpr',
    output: 'text',
    typedColumns: false,
    elides: true,
    encode: (node) => printSexpr(node, options),
    decode: (text: string, o: DecodeOptions = {}) => parseSexpr(text, o),
  };
}

export const SexprCodec = sexprCodec();
