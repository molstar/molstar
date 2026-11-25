/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The builder's engine (_spec/README.md "Builder"): a spec written in TypeScript, typed from a language's own
// description. Nothing is generated: a language is `defineLanguage<const L>`, so its kinds, their parents, ids and
// param types are in its type, and the types below read them: a kind's method exists on exactly the bodies of the
// kinds it may sit under (`parents`), takes an id exactly when the kind says so (`ids`), and its params are the
// declared ones, required or not, of the declared type. It builds the normalized tree the read pipeline gives: a column
// param holds its typed array, whether the builder was handed one or a plain array. The engine knows no language
// (a language's own builder, with handles and its vocabulary, is written on it). It normalizes as it builds; `print` is
// the codecs' (`printSexpr`).

import {
  type KindSpec,
  type Language,
  normalize,
  type ParamValues,
  type UnitsHook,
  validate,
} from '../schema/index.js';
import { assertPrintable, IdHint, IdPattern, type Issue, type Node, printAddress, type Value } from '../tree/index.js';

type Kinds = Readonly<Record<string, KindSpec>>;

// ---- the types: read from the description ----

/** A number with a unit, authored: a number, or a quantity (`4.2nm`) the units hook reads. */
export type Quantity = number | `${number}${string}`;

/** The params of a node of kind `C`: the schema's `ParamValues`, quantities allowed (no params: no keys). */
export type ParamsFor<K extends Kinds, C extends keyof K> = K[C] extends { params: infer F }
  ? ParamValues<F, Quantity>
  : Readonly<Record<string, never>>;

type HasRequired<P> = {} extends P ? false : true;

/** The kinds that may sit under a node of kind `P` (their `parents`). */
export type ChildKinds<K extends Kinds, P extends string> = {
  [C in keyof K & string]: P extends K[C]['parents'][number] ? C : never;
}[keyof K & string];

/**
 * What a kind's callback builds into: the kinds that sit under it. `A` says whether the node it builds into has an
 * address (the root has; a node with an id under one has; a node without an id, or in a slot, has not).
 */
export type NodeBuilder<K extends Kinds, P extends string, A extends boolean = true> = {
  readonly [C in ChildKinds<K, P>]: KindFn<K, C, A>;
};

type Body<K extends Kinds, C extends string, A extends boolean> = (node: NodeBuilder<K, C, A>) => void;

/** A handle: the new node's address when its parent has one, which a language's builder can wrap. */
type Handle<A extends boolean> = A extends true ? string : undefined;

type WithBody<K extends Kinds, C extends string, P> =
  HasRequired<P> extends true
    ? (params: P, body?: Body<K, C, false>) => undefined
    : (params?: P, body?: Body<K, C, false>) => undefined;

type WithIdAndBody<K extends Kinds, C extends string, P, A extends boolean> =
  HasRequired<P> extends true
    ? (id: string, params: P, body?: Body<K, C, A>) => Handle<A>
    : (id: string, params?: P, body?: Body<K, C, A>) => Handle<A>;

/**
 * A kind's method: `id` (required where the kind says `ids: 'required'`, absent where `'forbidden'`, either way
 * otherwise), its params, a body that adds its children. It returns the new node's address (`/data/p53`) when the
 * node has one: an id, under a parent with an address. A data form takes its positional values and its params.
 */
export type KindFn<K extends Kinds, C extends keyof K & string, A extends boolean = true> = K[C] extends {
  data: unknown;
}
  ? (args: readonly Value[], params?: ParamsFor<K, C>) => undefined
  : K[C] extends { ids: 'required' }
    ? WithIdAndBody<K, C, ParamsFor<K, C>, A>
    : K[C] extends { ids: 'forbidden' }
      ? WithBody<K, C, ParamsFor<K, C>>
      : WithBody<K, C, ParamsFor<K, C>> & WithIdAndBody<K, C, ParamsFor<K, C>, A>;

/**
 * A kind as a value: the node itself, built at once. `A`: whether its body builds addressable nodes (the root's
 * does; a slot's node has no address).
 */
export type FormFn<K extends Kinds, C extends keyof K & string, A extends boolean = false> = K[C] extends {
  data: unknown;
}
  ? (args: readonly Value[], params?: ParamsFor<K, C>) => Node
  : K[C] extends { ids: 'required' }
    ? FormWithId<K, C, A>
    : K[C] extends { ids: 'forbidden' }
      ? FormBare<K, C, A>
      : FormBare<K, C, A> & FormWithId<K, C, A>;

type FormBare<K extends Kinds, C extends keyof K & string, A extends boolean> =
  HasRequired<ParamsFor<K, C>> extends true
    ? (params: ParamsFor<K, C>, body?: Body<K, C, A>) => Node
    : (params?: ParamsFor<K, C>, body?: Body<K, C, A>) => Node;
type FormWithId<K extends Kinds, C extends keyof K & string, A extends boolean> =
  HasRequired<ParamsFor<K, C>> extends true
    ? (id: string, params: ParamsFor<K, C>, body?: Body<K, C, A>) => Node
    : (id: string, params?: ParamsFor<K, C>, body?: Body<K, C, A>) => Node;

export interface Builder<K extends Kinds, R extends keyof K & string> {
  readonly language: Language;
  /** The tree of a spec: the root kind, its params and a body that adds its children; normalized. */
  readonly root: FormFn<K, R, true>;
  /** Every kind of the language as a node built at once (what a slot holds); normalized. */
  readonly forms: { readonly [C in keyof K & string]: FormFn<K, C> };
  /** The description's issues for a tree (none: valid), the slots' languages resolved by `languages`. */
  check(node: Node): Issue[];
}

// ---- the runtime ----

interface Draft {
  readonly kind: string;
  id?: string;
  params?: Record<string, Value>;
  args?: readonly Value[];
  readonly children: (Draft | Node)[];
  /** The ids of its children, so a duplicate is refused where it is written (a Set, never a scan). */
  readonly ids: Set<string>;
  /** Its address (ids from the root), when it has one: the root, or an id under a node with one. */
  readonly address: readonly string[] | undefined;
}

const isDraft = (x: Draft | Node): x is Draft => 'ids' in x;

function cleanParams(params: Readonly<Record<string, unknown>> | undefined): Record<string, Value> | undefined {
  if (!params) return undefined;
  const out: Record<string, Value> = {};
  let any = false;
  for (const key of Object.keys(params)) {
    const value = params[key];
    if (value === undefined) continue;
    out[key] = value as Value;
    any = true;
  }
  return any ? out : undefined;
}

function finish(draft: Draft | Node): Node {
  if (!isDraft(draft)) return draft;
  const node: { kind: string; id?: string; params?: Record<string, Value>; children?: Node[]; args?: Value[] } = {
    kind: draft.kind,
  };
  if (draft.id !== undefined) node.id = draft.id;
  if (draft.params) node.params = draft.params;
  if (draft.args) node.args = [...draft.args];
  if (draft.children.length) node.children = draft.children.map(finish);
  assertPrintable(node);
  return node;
}

export interface BuilderOptions {
  /** The languages the slots name, by id (`validate`'s and `normalize`'s). */
  readonly languages?: (id: string) => Language | undefined;
  /** The tie to a units table, for `check`. */
  readonly units?: UnitsHook;
}

/**
 * The builder of a language. Its type says what the language's kinds allow, read from the language's own type, so a
 * `defineLanguage` result needs no annotation.
 */
export function createBuilder<const L extends Language>(
  language: L,
  options?: BuilderOptions,
): Builder<L['kinds'], L['root'] & keyof L['kinds'] & string>;
export function createBuilder(language: Language, options: BuilderOptions = {}): unknown {
  const kinds: Kinds = language.kinds;
  /** The draft of a node of `kind` from a method's arguments: `(id?, params?, body?)` or a data form's `(args, params?)`. */
  function draftOf(kind: string, argv: readonly unknown[], parent: readonly string[] | undefined, root = false): Draft {
    const spec = kinds[kind];
    const id = typeof argv[0] === 'string' && !spec.data ? argv[0] : undefined;
    // the root is `/` (its own id is no part of addresses); a node is addressable with an id under an addressable one
    const address = root ? [] : parent && id !== undefined ? [...parent, id] : undefined;
    const draft: Draft = { kind, children: [], ids: new Set(), address };
    if (spec.data) {
      if (!Array.isArray(argv[0])) throw new Error(`The data form ${kind} takes its values as a list`);
      draft.args = argv[0] as Value[];
      draft.params = cleanParams(argv[1] as Record<string, unknown> | undefined);
      return draft;
    }
    let i = 0;
    const ids = spec.ids ?? 'optional';
    if (typeof argv[0] === 'string') {
      if (ids === 'forbidden') throw new Error(`The kind ${kind} takes no id`);
      if (!IdPattern.test(argv[0])) throw new Error(`The id of ${kind} is not an id: ${IdHint}`);
      draft.id = argv[0];
      i = 1;
    } else if (ids === 'required') throw new Error(`The kind ${kind} needs an id`);
    const params = argv[i];
    if (params !== undefined && (typeof params !== 'object' || params === null || Array.isArray(params)))
      throw new Error(`The params of ${kind} are an object`);
    draft.params = cleanParams(params as Record<string, unknown> | undefined);
    const body = argv[i + 1];
    if (typeof body === 'function') (body as (n: object) => void)(bodyOf(draft));
    else if (body !== undefined) throw new Error(`The body of ${kind} is a function`);
    return draft;
  }

  /**
   * What a body offers: a method for each kind that sits under `parent`; any other kind's name says where it does sit
   * (a call on a kind that is no child is the builder's error, not a bare TypeError).
   */
  function bodyOf(parent: Draft): object {
    return new Proxy(
      {},
      {
        get(_, name) {
          if (typeof name !== 'string' || !Object.hasOwn(kinds, name)) return undefined;
          const where = kinds[name].parents;
          if (!where.includes(parent.kind)) {
            return () => {
              throw new Error(
                `The kind ${name} does not sit under ${parent.kind}: ${where.length ? `it sits under ${where.join(', ')}` : 'it is a root or a slot only'}`,
              );
            };
          }
          return (...argv: unknown[]) => {
            const draft = draftOf(name, argv, parent.address);
            if (draft.id !== undefined) {
              if (parent.ids.has(draft.id))
                throw new Error(`The id ${draft.id} is used twice under one ${parent.kind}`);
              parent.ids.add(draft.id);
            }
            parent.children.push(draft);
            return draft.address ? printAddress(draft.address) : undefined;
          };
        },
      },
    );
  }

  const normalizeOptions = options.languages ? { languages: options.languages } : {};
  /** A node built at once: the root (its body builds addressable nodes, the root being `/`), or a slot's node. */
  const built =
    (kind: string, root = false) =>
    (...argv: unknown[]) =>
      normalize(language, finish(draftOf(kind, argv, undefined, root)), normalizeOptions);
  const forms: Record<string, unknown> = {};
  for (const kind of Object.keys(kinds)) forms[kind] = built(kind);

  return {
    language,
    root: built(language.root, true) as never,
    forms: forms as never,
    check: (node: Node) => validate(language, node, { languages: options.languages, units: options.units }),
  };
}
