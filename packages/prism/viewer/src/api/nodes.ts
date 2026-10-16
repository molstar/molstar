/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type State, type StateObjectCell, StateSelection, StateTransform } from '@molstar/core/state';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { type ApiContext, type ApiFailure, ApiResult } from './area.js';
import { mergeParams } from './params.js';

// Node selectors (api.md §1 "Nodes", nodes.md): a node or a subtree of the plugin's state tree, said the way CSS
// says elements: `structure[label="1CBS"] component[label=Ligand] > representation[type=cartoon]`, `#ref`, `:hidden`.
// A selector compiles to a `StateSelection.Builder` (Mol*'s own monadic selection), so an operation that takes a `node`
// acts on every cell it matches, and a script can keep composing the builder. In scripts, `node(...)` writes the same
// text (`node('structure', { label: '1CBS' }).find('representation')`); the text is what the history records.

/** The kinds a selector names; `*` is any. */
export const NodeKinds = ['data', 'trajectory', 'model', 'structure', 'component', 'representation', 'volume'] as const;
export type NodeKind = (typeof NodeKinds)[number];

/** The attributes a selector tests (`[name op value]`); text compares without case, except `ref`. */
export const NodeAttributes = ['label', 'type', 'ref', 'tag', 'transformer'] as const;
export type NodeAttribute = (typeof NodeAttributes)[number];

type AttrOp = '=' | '*=' | '^=' | '$=' | 'has';
type Pseudo = 'hidden' | 'visible' | 'error';

interface Compound {
  readonly kind?: NodeKind;
  readonly attrs: readonly { readonly name: NodeAttribute; readonly op: AttrOp; readonly value: string }[];
  readonly pseudos: readonly Pseudo[];
  readonly not: readonly Compound[];
}

interface Complex {
  /** The compounds left to right; `combinators[i]` joins `parts[i]` and `parts[i + 1]`. */
  readonly parts: readonly Compound[];
  readonly combinators: readonly ('descendant' | 'child')[];
}

/** A parsed selector: one or more complex selectors, the cells any of them matches. */
export interface NodeSelector {
  readonly text: string;
  readonly alternatives: readonly Complex[];
}

const NodeQueryBrand = Symbol.for('molstar.prism.node-query');

/** What `node(...)` builds in scripts: the selector's text, chainable; an operation's `node` takes it as the text. */
export interface NodeQuery {
  readonly [NodeQueryBrand]: true;
  readonly text: string;
  /** Its children (`>`) of a kind; a decorator between (a transform) does not count. */
  child(kind?: NodeKind | '*', attrs?: NodeAttrs): NodeQuery;
  /** Anything below it (` `) of a kind. */
  find(kind?: NodeKind | '*', attrs?: NodeAttrs): NodeQuery;
  /** Only the hidden ones (`:hidden`). */
  hidden(): NodeQuery;
  /** Only the shown ones (`:visible`). */
  visible(): NodeQuery;
  /** These or those (`,`). */
  or(other: NodeQuery | string): NodeQuery;
  /** Mol*'s selection builder for it, to keep composing. */
  compile(): StateSelection.Builder;
  toString(): string;
  toJSON(): string;
}

/** Attribute tests for `node(...)`: equal text; `{ label: { contains: 'lig' } }` and the like for the other operators. */
export type NodeAttrs = Partial<
  Record<NodeAttribute, string | { contains: string } | { startsWith: string } | { endsWith: string }>
>;

const quote = (value: string) => (/^[\w-]+$/.test(value) ? value : JSON.stringify(value));

function compoundText(kind: NodeKind | '*' | undefined, attrs?: NodeAttrs): string {
  let text = kind ?? '*';
  for (const [name, test] of Object.entries(attrs ?? {})) {
    if (test === undefined) continue;
    if (typeof test === 'string')
      text += name === 'ref' && kind === undefined ? `#${test}` : `[${name}=${quote(test)}]`;
    else if ('contains' in test) text += `[${name}*=${quote(test.contains)}]`;
    else if ('startsWith' in test) text += `[${name}^=${quote(test.startsWith)}]`;
    else text += `[${name}$=${quote(test.endsWith)}]`;
  }
  return text.length > 1 && text.startsWith('*') ? text.slice(1) : text;
}

function nodeQuery(text: string): NodeQuery {
  // Each alternative of a list goes on with the step (`a, b` then `> c` is `a > c, b > c`).
  const each = (step: string) =>
    nodeQuery(
      splitTop(text)
        .map((t) => `${t}${step}`)
        .join(', '),
    );
  return {
    [NodeQueryBrand]: true,
    text,
    child: (kind, attrs) => each(` > ${compoundText(kind, attrs)}`),
    find: (kind, attrs) => each(` ${compoundText(kind, attrs)}`),
    hidden: () => each(':hidden'),
    visible: () => each(':visible'),
    or: (other) => nodeQuery(`${text}, ${typeof other === 'string' ? other : other.text}`),
    compile: () => {
      const parsed = parseNodeSelector(text);
      if ('error' in parsed) throw new Error(parsed.error);
      return compileNodeSelector(parsed);
    },
    toString: () => text,
    toJSON: () => text,
  };
}

/**
 * A node selector in scripts: `node('structure', { label: '1CBS' }).child('component').find('representation')`.
 * `node.ref(ref)` is one cell by its ref; `node.parse(text)` takes selector text as is.
 */
export function node(kind?: NodeKind | '*', attrs?: NodeAttrs): NodeQuery {
  return nodeQuery(compoundText(kind, attrs));
}
node.ref = (ref: string): NodeQuery => nodeQuery(`#${ref}`);
node.parse = (text: string): NodeQuery => nodeQuery(text.trim());

export function isNodeQuery(value: unknown): value is NodeQuery {
  return !!value && typeof value === 'object' && (value as NodeQuery)[NodeQueryBrand] === true;
}

/** The text of what an operation was given as a node: selector text, or a `node(...)` query. */
export function nodeText(value: unknown): string {
  if (isNodeQuery(value)) return value.text;
  return typeof value === 'string' ? value.trim() : '';
}

// ---- Parsing ----

/** The top-level alternatives of a list (commas outside brackets, parentheses and quotes). */
function splitTop(text: string): string[] {
  const out: string[] = [];
  let depth = 0;
  let quoteChar = '';
  let start = 0;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoteChar) {
      if (c === '\\') i++;
      else if (c === quoteChar) quoteChar = '';
    } else if (c === '"' || c === "'") quoteChar = c;
    else if (c === '[' || c === '(') depth++;
    else if (c === ']' || c === ')') depth--;
    else if (c === ',' && depth === 0) {
      out.push(text.slice(start, i).trim());
      start = i + 1;
    }
  }
  out.push(text.slice(start).trim());
  return out;
}

const NameChar = /[\w-]/;
// A ref may hold '=' (the root's is `-=root=-`).
const RefChar = /[\w\-=]/;
const AttrOps: readonly AttrOp[] = ['*=', '^=', '$=', '='];

class SelectorParser {
  private i = 0;
  constructor(private readonly text: string) {}

  private fail(what: string): never {
    throw new SyntaxError(`${what} at ${this.i + 1} in '${this.text}'.`);
  }

  private get c() {
    return this.text[this.i];
  }

  private spaces(): boolean {
    const start = this.i;
    while (this.i < this.text.length && /\s/.test(this.c)) this.i++;
    return this.i > start;
  }

  private ident(what: string, chars = NameChar): string {
    const start = this.i;
    while (this.i < this.text.length && chars.test(this.c)) this.i++;
    if (this.i === start) this.fail(`Expected ${what}`);
    return this.text.slice(start, this.i);
  }

  private value(): string {
    const q = this.c;
    if (q !== '"' && q !== "'") return this.ident('a value');
    let out = '';
    this.i++;
    while (this.i < this.text.length && this.c !== q) {
      if (this.c === '\\') this.i++;
      out += this.c;
      this.i++;
    }
    if (this.c !== q) this.fail('Unclosed quote');
    this.i++;
    return out;
  }

  list(): Complex[] {
    const out: Complex[] = [];
    while (true) {
      this.spaces();
      out.push(this.complex());
      this.spaces();
      if (this.i >= this.text.length || this.c === ')') return out;
      if (this.c !== ',') this.fail(`Unexpected '${this.c}'`);
      this.i++;
    }
  }

  private complex(): Complex {
    const parts: Compound[] = [this.compound()];
    const combinators: ('descendant' | 'child')[] = [];
    while (true) {
      const spaced = this.spaces();
      if (this.i >= this.text.length || this.c === ',' || this.c === ')') return { parts, combinators };
      if (this.c === '>') {
        this.i++;
        this.spaces();
        combinators.push('child');
      } else if (spaced) combinators.push('descendant');
      else this.fail(`Unexpected '${this.c}'`);
      parts.push(this.compound());
    }
  }

  compound(): Compound {
    let kind: NodeKind | undefined;
    const attrs: { name: NodeAttribute; op: AttrOp; value: string }[] = [];
    const pseudos: Pseudo[] = [];
    const not: Compound[] = [];
    const start = this.i;
    if (this.c === '*') this.i++;
    else if (this.c && /[a-z]/i.test(this.c)) {
      const name = this.ident('a kind').toLowerCase();
      if (!NodeKinds.includes(name as NodeKind)) {
        this.i = start;
        this.fail(`Unknown kind '${name}' (kinds: ${NodeKinds.join(', ')}, *)`);
      }
      kind = name as NodeKind;
    }
    while (this.i < this.text.length) {
      const c = this.c;
      if (c === '#') {
        this.i++;
        attrs.push({ name: 'ref', op: '=', value: this.ident('a ref', RefChar) });
      } else if (c === '.') {
        this.i++;
        attrs.push({ name: 'tag', op: 'has', value: this.ident('a tag') });
      } else if (c === '[') {
        this.i++;
        this.spaces();
        const name = this.ident('an attribute').toLowerCase() as NodeAttribute;
        if (!NodeAttributes.includes(name)) this.fail(`Unknown attribute '${name}' (${NodeAttributes.join(', ')})`);
        this.spaces();
        const op = AttrOps.find((o) => this.text.startsWith(o, this.i));
        if (!op) this.fail(`Expected =, *=, ^= or $=`);
        this.i += op.length;
        this.spaces();
        const value = this.value();
        this.spaces();
        if (this.c !== ']') this.fail(`Expected ']'`);
        this.i++;
        attrs.push({ name, op: name === 'tag' && op === '=' ? 'has' : op, value });
      } else if (c === ':') {
        this.i++;
        const name = this.ident('a state').toLowerCase();
        if (name === 'not') {
          if (this.c !== '(') this.fail(`Expected '(' after :not`);
          this.i++;
          this.spaces();
          not.push(this.compound());
          this.spaces();
          if ((this.c as string) !== ')') this.fail(`:not takes one kind with its tests; expected ')'`);
          this.i++;
        } else if (name === 'hidden' || name === 'visible' || name === 'error') pseudos.push(name);
        else this.fail(`Unknown state ':${name}' (:hidden, :visible, :error, :not(…))`);
      } else break;
    }
    if (this.i === start) this.fail(this.i >= this.text.length ? 'Expected a selector' : `Unexpected '${this.c}'`);
    return { kind, attrs, pseudos, not };
  }

  done() {
    if (this.i < this.text.length) this.fail(`Unexpected '${this.c}'`);
  }
}

/** Parses selector text; a failure says where and what it expected. */
export function parseNodeSelector(text: string): NodeSelector | { error: string } {
  const trimmed = text.trim();
  if (!trimmed) return { error: 'An empty selector.' };
  try {
    const parser = new SelectorParser(trimmed);
    const alternatives = parser.list();
    parser.done();
    return { text: trimmed, alternatives };
  } catch (e) {
    return { error: e instanceof Error ? e.message : String(e) };
  }
}

// ---- Matching ----

const isDecorator = (cell: StateObjectCell) => !!cell.transform.transformer.definition.isDecorator;

/** The cell's parent, past decorators (a transform folded into the structure it moves). */
function baseParent(state: State, cell: StateObjectCell): StateObjectCell | undefined {
  let ref = cell.transform.parent;
  while (ref !== cell.transform.ref) {
    const parent = state.cells.get(ref);
    if (!parent || ref === StateTransform.RootRef) return parent;
    if (!isDecorator(parent)) return parent;
    if (parent.transform.parent === ref) return undefined;
    ref = parent.transform.parent;
  }
  return undefined;
}

function hasStructureAbove(state: State, cell: StateObjectCell): boolean {
  for (let p = baseParent(state, cell); p && p.transform.ref !== StateTransform.RootRef; p = baseParent(state, p)) {
    if (SO.Molecule.Structure.is(p.obj) && !isDecorator(p)) return true;
  }
  return false;
}

/** The cell's kind, as the Structure panel's hierarchy reads it (a structure inside a structure is a component). */
export function nodeKind(state: State, cell: StateObjectCell): NodeKind | undefined {
  const obj = cell.obj;
  if (!obj || isDecorator(cell)) return undefined;
  if (SO.isRepresentation3D(obj)) return 'representation';
  if (SO.Molecule.Structure.is(obj)) return hasStructureAbove(state, cell) ? 'component' : 'structure';
  if (SO.Molecule.Model.is(obj)) return 'model';
  if (SO.Molecule.Trajectory.is(obj)) return 'trajectory';
  if (SO.Volume.Data.is(obj)) return 'volume';
  if (obj.type.typeClass === 'Data') return 'data';
  return undefined;
}

/** A representation's type, a component's static kind (ligand, water…). */
function typeOf(cell: StateObjectCell): string | undefined {
  const t = (cell.transform.params as { type?: unknown } | undefined)?.type as
    | { name?: unknown; params?: unknown }
    | undefined;
  if (!t || typeof t !== 'object' || typeof t.name !== 'string') return undefined;
  if (t.name === 'static' && typeof t.params === 'string') return t.params;
  return t.name;
}

function attrValues(cell: StateObjectCell, name: NodeAttribute): readonly string[] {
  switch (name) {
    case 'label':
      return [cell.obj?.label ?? cell.transform.transformer.definition.display.name];
    case 'type': {
      const t = typeOf(cell);
      return t ? [t] : [];
    }
    case 'ref':
      return [cell.transform.ref];
    case 'tag':
      return cell.transform.tags ?? [];
    case 'transformer':
      return [cell.transform.transformer.id];
  }
}

function testAttr(cell: StateObjectCell, a: Compound['attrs'][number]): boolean {
  const exact = a.name === 'ref';
  const want = exact ? a.value : a.value.toLowerCase();
  return attrValues(cell, a.name).some((v) => {
    const have = exact ? v : v.toLowerCase();
    switch (a.op) {
      case '=':
      case 'has':
        return have === want;
      case '*=':
        return have.includes(want);
      case '^=':
        return have.startsWith(want);
      case '$=':
        return have.endsWith(want);
    }
  });
}

function matches(state: State, cell: StateObjectCell, c: Compound): boolean {
  if (cell.transform.ref === StateTransform.RootRef)
    return c.attrs.some((a) => a.name === 'ref') && c.attrs.every((a) => testAttr(cell, a));
  if (cell.state.isGhost) return false;
  if (c.kind && nodeKind(state, cell) !== c.kind) return false;
  for (const a of c.attrs) if (!testAttr(cell, a)) return false;
  for (const p of c.pseudos) {
    if (p === 'hidden' && !cell.state.isHidden) return false;
    if (p === 'visible' && cell.state.isHidden) return false;
    if (p === 'error' && cell.status !== 'error') return false;
  }
  for (const n of c.not) if (matches(state, cell, n)) return false;
  return true;
}

/** Every cell below `ref`, in tree order. */
function below(state: State, ref: string, out: StateObjectCell[]): void {
  for (const c of state.tree.children.get(ref).toArray()) {
    const cell = state.cells.get(c);
    if (!cell) continue;
    out.push(cell);
    below(state, c, out);
  }
}

/** Its children, a decorator between counted as part of its base (`structure > component` through a transform). */
function childrenOf(state: State, ref: string, out: StateObjectCell[]): void {
  for (const c of state.tree.children.get(ref).toArray()) {
    const cell = state.cells.get(c);
    if (!cell) continue;
    out.push(cell);
    if (isDecorator(cell)) childrenOf(state, c, out);
  }
}

function evaluate(state: State, complex: Complex): StateObjectCell[] {
  const all = [state.cells.get(StateTransform.RootRef)!];
  below(state, StateTransform.RootRef, all);
  let current = all.filter((cell) => matches(state, cell, complex.parts[0]));
  for (let i = 1; i < complex.parts.length; i++) {
    const part = complex.parts[i];
    const seen = new Set<string>();
    const next: StateObjectCell[] = [];
    for (const cell of current) {
      const candidates: StateObjectCell[] = [];
      if (complex.combinators[i - 1] === 'child') childrenOf(state, cell.transform.ref, candidates);
      else below(state, cell.transform.ref, candidates);
      for (const c of candidates) {
        if (seen.has(c.transform.ref) || !matches(state, c, part)) continue;
        seen.add(c.transform.ref);
        next.push(c);
      }
    }
    current = next;
  }
  return current;
}

/** Mol*'s selection builder for a parsed selector: the cells any alternative matches, each once, in tree order. */
export function compileNodeSelector(selector: NodeSelector): StateSelection.Builder {
  return StateSelection.flatMap(StateSelection.Generators.root, (_, state) => {
    if (selector.alternatives.length === 1) return evaluate(state, selector.alternatives[0]);
    const found = new Set<string>();
    for (const alt of selector.alternatives) for (const cell of evaluate(state, alt)) found.add(cell.transform.ref);
    const ordered = [state.cells.get(StateTransform.RootRef)!];
    below(state, StateTransform.RootRef, ordered);
    return ordered.filter((cell) => found.has(cell.transform.ref));
  });
}

/** The cells a selector (text or a `node(...)` query) matches in a state; a failure for text that does not parse. */
export function selectNodes(state: State, selector: string | NodeQuery): StateObjectCell[] | ApiFailure {
  const text = nodeText(selector);
  // A ref alone is that cell (what state.get and the results list).
  const byRef = state.cells.get(text);
  if (byRef) return [byRef];
  const parsed = parseNodeSelector(text);
  // A plain word that is no kind: a ref that is not (or no longer) there.
  if ('error' in parsed && /^[\w\-=]+$/.test(text)) {
    return ApiResult.fail(
      'not-found',
      `There is no node '${text}'.`,
      `A ref (state.get lists them), or a kind: ${NodeKinds.join(', ')}.`,
    );
  }
  if ('error' in parsed) {
    return ApiResult.fail(
      'invalid',
      `The node selector cannot be read: ${parsed.error}`,
      'For example structure[label="1CBS"] > representation[type=cartoon], or #ref (state.get lists refs).',
    );
  }
  return compileNodeSelector(parsed).select(state);
}

/** A cell for a line: "Cartoon (abc123)". */
export function nodeLine(cell: StateObjectCell): string {
  return `${cell.obj?.label ?? cell.transform.transformer.definition.display.name} (${cell.transform.ref})`;
}

/**
 * The one cell an operation acts on: none or several fail, the several listed so the caller can pick (`#ref`).
 * `kind` narrows what counts and names it in the message.
 */
export function selectOneNode(
  state: State,
  selector: string | NodeQuery,
  kind?: NodeKind,
): StateObjectCell | ApiFailure {
  const cells = selectNodes(state, selector);
  if (ApiResult.isFailure(cells)) return cells;
  const of = kind ? cells.filter((c) => nodeKind(state, c) === kind) : cells;
  const what = kind ?? 'node';
  const text = nodeText(selector);
  if (of.length === 0) {
    return ApiResult.fail(
      'not-found',
      `No ${what} matches '${text}'${cells.length ? ` (${cells.length} other node${cells.length > 1 ? 's' : ''} do)` : ''}.`,
      'state.select lists what a selector matches; state.get lists the tree.',
    );
  }
  if (of.length > 1) {
    const shown = of.slice(0, 8).map(nodeLine).join(', ');
    return ApiResult.fail(
      'ambiguous',
      `${of.length} ${what}s match '${text}': ${shown}${of.length > 8 ? ', …' : ''}.`,
      'Narrow it, or pick one by its ref (#ref).',
    );
  }
  return of[0];
}

/**
 * The cells an operation acts on (api.md §1 "Nodes"): every match of the selector that `accept` takes (a
 * representation, a volume…), mapped to what the operation works with. None fails, saying what matched instead.
 */
export function nodesFor<T>(
  state: State,
  selector: string | NodeQuery,
  what: string,
  accept: (cell: StateObjectCell) => T | undefined,
): T[] | ApiFailure {
  const text = nodeText(selector);
  if (!text) {
    return ApiResult.fail('invalid', `Which ${what}?`, 'Give a node selector (representation[type=cartoon]) or a ref.');
  }
  const cells = selectNodes(state, text);
  if (ApiResult.isFailure(cells)) return cells;
  const out: T[] = [];
  for (const cell of cells) {
    const t = accept(cell);
    if (t !== undefined) out.push(t);
  }
  if (!out.length) {
    return ApiResult.fail(
      'not-found',
      `No ${what} matches '${text}'${cells.length ? ` (it matches ${cells.slice(0, 3).map(nodeLine).join(', ')}${cells.length > 3 ? ', …' : ''})` : ''}.`,
      'state.select lists what a selector matches; state.get lists the tree.',
    );
  }
  return out;
}

/**
 * The inspectors' forms on what a selector matched: each cell's params merged with the patch and updated where they
 * change; none changing is already so.
 */
export async function updateEach(
  ctx: ApiContext,
  refs: readonly { readonly cell: StateObjectCell }[],
  patch: Record<string, unknown>,
  what: string,
): Promise<ApiResult> {
  const changed: StateObjectCell[] = [];
  for (const { cell } of refs) {
    const current = cell.transform.params as Record<string, unknown>;
    const next = mergeParams(current, patch);
    if (JSON.stringify(next) === JSON.stringify(current)) continue;
    await ctx.plugin.state.updateTransform(ctx.plugin.state.data, cell.transform.ref, next, what);
    changed.push(cell);
  }
  const name = (list: readonly StateObjectCell[]) =>
    list.length === 1 ? nodeLine(list[0]) : `${list.length} ${what.toLowerCase()}s`;
  if (!changed.length)
    return ApiResult.already(`${name(refs.map((r) => r.cell))} ${refs.length > 1 ? 'are' : 'is'} already so.`);
  return ApiResult.ok(`Updated ${name(changed)}.`, { nodes: changed.map((c) => c.transform.ref) });
}
