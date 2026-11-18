/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The language builder (_spec/README.md, Language builder): the way most specs are written. It reads as the spec does
// (`b.data.sequence(...)`, `b.views.view(...)`), every call returns the address of what it adds, branded by the kind
// it names, and a reference param takes a handle of a kind it allows (or an address written out). `build` gives the
// normalized tree the generic engine of `@molstar/core/repr/builder` builds. Calls are recorded as they are written
// and replayed on the engine, so misuse that types do not reach (an id twice, a malformed id) is thrown where it is
// written.

import { createBuilder, type ParamsFor } from '@molstar/core/repr/builder';
import { IdHint, IdPattern, type Issue, type Node, printAddress } from '@molstar/core/repr/tree';
import { checkSeqSpec } from './check.js';
import { type SeqKinds, SeqLanguage, SeqLanguageVersion } from './language.js';

declare const KindBrand: unique symbol;

/** The address of a node of kind `K`, as the builder returns it: a string that says what it names. */
export type Address<K extends string = string> = string & { readonly [KindBrand]: K };

/** A builder that has an address: it stands where the address does. */
export interface Addressed<K extends string = string> {
  readonly address: Address<K>;
  toString(): string;
}

/** What a reference param takes: a handle of an allowed kind, or an address written out (any string that is no handle). */
export type RefValue<K extends string> = Address<K> | Addressed<K> | (string & { readonly [KindBrand]?: undefined });

type Fields<C extends keyof SeqKinds> = SeqKinds[C] extends { params: infer F } ? F : Record<never, never>;

/** The params of a kind as the builder takes them: the description's, its reference params typed by their `kinds`. */
export type Params<C extends keyof SeqKinds & string> = {
  [P in keyof ParamsFor<SeqKinds, C>]: P extends keyof Fields<C>
    ? Fields<C>[P] extends { type: 'ref'; kinds: readonly (infer K extends string)[] }
      ? RefValue<K>
      : ParamsFor<SeqKinds, C>[P]
    : ParamsFor<SeqKinds, C>[P];
};

export interface DataBuilder {
  source(id: string, params: Params<'source'>): Address<'source'>;
  sequence(id: string, params: Params<'sequence'>): Address<'sequence'>;
  ranges(id: string, params: Params<'ranges'>): Address<'ranges'>;
  values(id: string, params: Params<'values'>): Address<'values'>;
  map(id: string, params: Params<'map'>): Address<'map'>;
}

export interface ViewsBuilder {
  view(id: string, body?: (view: ViewBuilder) => void): Address<'view'>;
  view(id: string, params: Params<'view'>, body?: (view: ViewBuilder) => void): Address<'view'>;
}

export interface ViewBuilder extends Addressed<'view'> {
  section(id: string, params: Params<'section'>, body?: (section: SectionBuilder) => void): SectionBuilder;
}

export interface SectionBuilder extends Addressed<'section'> {
  track(id: string, params?: Params<'track'>, body?: (track: TrackBuilder) => void): TrackBuilder;
}

export interface TrackBuilder extends Addressed<'track'> {
  /** A layer and its address: its id is required, so hiding it is a one-key patch. */
  layer(id: string, params: Params<'layer'>): Address<'layer'>;
}

/** What a ref names, as the builder takes it: an address (a handle or written out), or a set of them. */
export type RefTo = RefValue<string> | readonly RefValue<string>[];

export interface RefsBuilder {
  /** A user-defined ref (`/refs/<id>`): a name for what an address, or a set of them, names. */
  ref(id: string, params: { readonly to: RefTo; readonly label?: string; readonly color?: string }): Address<'ref'>;
}

export interface SeqSpecBuilder {
  readonly data: DataBuilder;
  readonly views: ViewsBuilder;
  readonly refs: RefsBuilder;
  /** The spec as a normalized tree (print it with `printSexpr` from `@molstar/core/repr/codecs/sexpr`). */
  build(): Node;
  /** The issues of the spec so far: the description's, then the language's. */
  check(): Issue[];
  toString(): string;
}

export interface SeqSpecOptions {
  readonly title?: string;
  /** The language version the spec says; this one by default. */
  readonly version?: string;
}

interface Rec {
  readonly kind: string;
  readonly id?: string;
  readonly params?: Record<string, unknown>;
  readonly children: Rec[];
  readonly ids: Set<string>;
  readonly address: string;
}

type Engine = Record<string, (...args: unknown[]) => unknown>;

type Described = Readonly<
  Record<string, { readonly params?: Readonly<Record<string, { readonly type: string; readonly of?: unknown }>> }>
>;

/** Whether a param type holds references: a ref, or a union or list with one in it. */
function holdsRefs(t: { readonly type: string; readonly of?: unknown } | undefined): boolean {
  if (!t) return false;
  if (t.type === 'ref') return true;
  if (t.type === 'list') return holdsRefs(t.of as { type: string });
  if (t.type === 'union') return (t.of as { type: string }[]).some(holdsRefs);
  return false;
}

/** A handle read as its address, in a list too. */
function addressOfHandle(v: unknown): unknown {
  if (Array.isArray(v)) return v.map(addressOfHandle);
  return v && typeof v === 'object' && 'address' in v ? (v as Addressed).address : v;
}

/**
 * The params with each reference param's handle read as its address (in a set of them too); which params are
 * references is the description's (`type: 'ref'`, or a union or list holding one), never a list here.
 */
function addresses(kind: string, params: unknown): Record<string, unknown> | undefined {
  if (!params || typeof params !== 'object') return params as undefined;
  const fields = (SeqLanguage.kinds as Described)[kind]?.params;
  const out = { ...(params as Record<string, unknown>) };
  for (const key of Object.keys(out)) if (holdsRefs(fields?.[key])) out[key] = addressOfHandle(out[key]);
  return out;
}

/** A spec being written: `build()` is its normalized tree. */
export function createSeqSpec(options: SeqSpecOptions = {}): SeqSpecBuilder {
  const engine = createBuilder(SeqLanguage);
  const data: Rec = { kind: 'data', id: 'data', children: [], ids: new Set(), address: '/data' };
  const views: Rec = { kind: 'views', id: 'views', children: [], ids: new Set(), address: '/views' };
  const refs: Rec = { kind: 'refs', id: 'refs', children: [], ids: new Set(), address: '/refs' };

  /** A node added under `parent`: its id checked where it is written, its address the parent's and the id. */
  function add(parent: Rec, kind: string, id: string | undefined, params: unknown): Rec {
    if (id !== undefined) {
      if (!IdPattern.test(id)) throw new Error(`The id of ${kind} is not an id: ${IdHint}`);
      if (parent.ids.has(id)) throw new Error(`The id ${id} is used twice under one ${parent.kind}`);
      parent.ids.add(id);
    }
    const address = printAddress([...parent.address.split('/').filter(Boolean), ...(id !== undefined ? [id] : [])]);
    const rec: Rec = { kind, id, params: addresses(kind, params), children: [], ids: new Set(), address };
    parent.children.push(rec);
    return rec;
  }

  const dataKind = (kind: string) => (id: string, params: unknown) => add(data, kind, id, params).address;
  const dataBuilder = {
    source: dataKind('source'),
    sequence: dataKind('sequence'),
    ranges: dataKind('ranges'),
    values: dataKind('values'),
    map: dataKind('map'),
  } as unknown as DataBuilder;

  function trackBuilder(rec: Rec): TrackBuilder {
    return {
      address: rec.address as Address<'track'>,
      layer(id: string, params: unknown) {
        if (typeof id !== 'string') throw new Error('A layer needs an id');
        return add(rec, 'layer', id, params).address;
      },
      toString: () => rec.address,
    } as TrackBuilder;
  }
  function sectionBuilder(rec: Rec): SectionBuilder {
    return {
      address: rec.address as Address<'section'>,
      track(id, params, body) {
        const t = trackBuilder(add(rec, 'track', id, params));
        body?.(t);
        return t;
      },
      toString: () => rec.address,
    };
  }
  function viewBuilder(rec: Rec): ViewBuilder {
    return {
      address: rec.address as Address<'view'>,
      section(id, params, body) {
        const s = sectionBuilder(add(rec, 'section', id, params));
        body?.(s);
        return s;
      },
      toString: () => rec.address,
    };
  }
  const viewsBuilder = {
    view(id: string, a?: unknown, b?: unknown) {
      const [params, body] = typeof a === 'function' ? [undefined, a] : [a, b];
      const rec = add(views, 'view', id, params);
      (body as ((v: ViewBuilder) => void) | undefined)?.(viewBuilder(rec));
      return rec.address;
    },
  } as unknown as ViewsBuilder;

  const refsBuilder: RefsBuilder = {
    ref: (id, params) => add(refs, 'ref', id, params).address as Address<'ref'>,
  };

  /** Replays a recorded node on the engine's body. */
  function emit(body: Engine, rec: Rec): void {
    const inner = (b: Engine) => {
      for (const c of rec.children) emit(b, c);
    };
    if (rec.id !== undefined) body[rec.kind](rec.id, rec.params, inner);
    else body[rec.kind](rec.params, inner);
  }

  function build(): Node {
    return engine.root({ version: options.version ?? SeqLanguageVersion, title: options.title }, (root) => {
      const r = root as unknown as Engine;
      for (const container of [data, views, refs]) if (container.children.length) emit(r, container);
    });
  }

  return {
    data: dataBuilder,
    views: viewsBuilder,
    refs: refsBuilder,
    build,
    check: () => checkSeqSpec(build()),
    toString: () => '/',
  };
}
