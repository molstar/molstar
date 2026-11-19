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
import { type Patch, patchToNode } from '@molstar/core/repr/patch';
import { IdHint, IdPattern, type Issue, type Node, printAddress, type Value } from '@molstar/core/repr/tree';
import { checkSeqSpec } from './check.js';
import { type SeqKinds, SeqLanguage, SeqLanguageVersion, seqLanguages } from './language.js';
import { writeLayerTheme } from './theme-builder.js';

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
  alignment(id: string, params: Params<'alignment'>): Address<'alignment'>;
  tree(id: string, params: Params<'tree'>): Address<'tree'>;
  pairs(id: string, params: Params<'pairs'>): Address<'pairs'>;
  grid(id: string, params: Params<'grid'>): Address<'grid'>;
}

export interface ViewsBuilder {
  view(id: string, body?: (view: ViewBuilder) => void): Address<'view'>;
  view(id: string, params: Params<'view'>, body?: (view: ViewBuilder) => void): Address<'view'>;
  /** A shared scale (`/views/<id>`): layers whose themes name it share its domain and palette. */
  scale(id: string, params?: Params<'scale'>): Address<'scale'>;
}

export interface ViewBuilder extends Addressed<'view'> {
  section(id: string, params: Params<'section'>, body?: (section: SectionBuilder) => void): SectionBuilder;
  /** A navigator (`<view>/navigators/<id>`): its params, or a form's (`stripNavigator`, `ideogramNavigator`, …). */
  navigator(id: string, params: NavigatorParams): Address<'navigator'>;
}

/** The sections a navigator binds: one, or several (a handle or an address each). */
export type NavigatorOf = RefValue<'section'> | readonly RefValue<'section'>[];

/** A navigator's params as the builder takes them: its references as handles or addresses. */
export type NavigatorParams = Omit<Params<'navigator'>, 'of' | 'tracks'> & {
  readonly of: NavigatorOf;
  readonly tracks?: readonly RefValue<'track'>[];
};

type NavigatorExtra = Omit<NavigatorParams, 'of' | 'axes' | 'content'>;

/**
 * The navigator's forms as params (_spec/README.md, Language, Navigators): a strip of the axis, or of the section's
 * tracks small (`content: 'tracks'`), across the view.
 */
export function stripNavigator(
  of: NavigatorOf,
  extra: NavigatorExtra & { readonly content?: 'axis' | 'tracks' } = {},
): NavigatorParams {
  return { of, axes: 'positions', content: 'axis', ...extra };
}

/** An ideogram: a chromosome's stained bands (`bands`), on `system` when the section shows a region of it. */
export function ideogramNavigator(
  of: NavigatorOf,
  bands: NavigatorParams['bands'],
  extra: Omit<NavigatorExtra, 'bands'> = {},
): NavigatorParams {
  return { of, axes: 'positions', content: 'bands', bands, ...extra };
}

/** A vertical overview of an alignment's rows: the rows summarized down a column, the window full width. */
export function overviewNavigator(of: NavigatorOf, extra: NavigatorExtra = {}): NavigatorParams {
  return { of, axes: 'rows', content: 'rows', place: 'right', ...extra };
}

/** A miniature: the section's tracks and rows, both axes whole. */
export function miniatureNavigator(of: NavigatorOf, extra: NavigatorExtra = {}): NavigatorParams {
  return { of, axes: 'both', content: 'tracks', ...extra };
}

export interface SectionBuilder extends Addressed<'section'> {
  track(id: string, params?: Params<'track'>, body?: (track: TrackBuilder) => void): TrackBuilder;
  /** A block: a piece of the section on a system of its own (a chain, a ligand), after the blocks before it. */
  block(id: string, params: Params<'block'>): Address<'block'>;
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

/** A story step's params: a patch (no base) or a full state, and its start, view and timing. */
export interface SnapshotParams {
  readonly patch?: Patch | Node;
  readonly state?: Node;
  readonly label?: string;
  /** The step the patch applies over (`/snapshots/<id>`), instead of the spec. */
  readonly from?: string;
  readonly keep?: 'full' | 'spec';
  readonly view?: Value;
  readonly duration?: number;
  readonly easing?: 'linear' | 'ease-in' | 'ease-out' | 'ease-in-out';
  readonly linger?: number;
}

export interface SnapshotsBuilder {
  /**
   * A story step (`/snapshots/<id>`): a stored patch with no base, shown over whatever the spec is. Its addresses are
   * refs of the spec (a rename or a move rewrites them).
   */
  snapshot(id: string, params: SnapshotParams): Address<'snapshot'>;
}

export interface SeqSpecBuilder {
  readonly data: DataBuilder;
  readonly views: ViewsBuilder;
  readonly refs: RefsBuilder;
  readonly snapshots: SnapshotsBuilder;
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
  const engine = createBuilder(SeqLanguage, { languages: seqLanguages });
  const data: Rec = { kind: 'data', id: 'data', children: [], ids: new Set(), address: '/data' };
  const views: Rec = { kind: 'views', id: 'views', children: [], ids: new Set(), address: '/views' };
  const refs: Rec = { kind: 'refs', id: 'refs', children: [], ids: new Set(), address: '/refs' };
  const snapshots: Rec = { kind: 'snapshots', id: 'snapshots', children: [], ids: new Set(), address: '/snapshots' };

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
    alignment: dataKind('alignment'),
    tree: dataKind('tree'),
    pairs: dataKind('pairs'),
    grid: dataKind('grid'),
  } as unknown as DataBuilder;

  /** The inline rows of the data node an address names (a ranges node's), for the layer's theme. */
  const rowsOf = (handle: unknown): Readonly<Record<string, unknown>> | undefined => {
    if (handle === undefined) return undefined;
    const a = addressOfHandle(handle);
    const rows = data.children.find((r) => r.address === a)?.params?.rows;
    return rows && typeof rows === 'object' ? (rows as Record<string, unknown>) : undefined;
  };

  function trackBuilder(rec: Rec): TrackBuilder {
    return {
      address: rec.address as Address<'track'>,
      layer(id: string, params: unknown) {
        if (typeof id !== 'string') throw new Error('A layer needs an id');
        const p = params as Record<string, unknown>;
        return add(rec, 'layer', id, writeLayerTheme(p, rowsOf(p?.data))).address;
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
      block: (id, params) => add(rec, 'block', id, params).address as Address<'block'>,
      toString: () => rec.address,
    };
  }
  function viewBuilder(rec: Rec): ViewBuilder {
    let navigators: Rec | undefined;
    return {
      address: rec.address as Address<'view'>,
      section(id, params, body) {
        const s = sectionBuilder(add(rec, 'section', id, params));
        body?.(s);
        return s;
      },
      navigator(id, params) {
        navigators ??= add(rec, 'navigators', 'navigators', undefined);
        return add(navigators, 'navigator', id, params).address as Address<'navigator'>;
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
    scale(id: string, params?: unknown) {
      return add(views, 'scale', id, params ?? {}).address;
    },
  } as unknown as ViewsBuilder;

  const refsBuilder: RefsBuilder = {
    ref: (id, params) => add(refs, 'ref', id, params).address as Address<'ref'>,
  };

  const snapshotsBuilder: SnapshotsBuilder = {
    snapshot: (id, { patch, ...rest }) => {
      if (patch && 'base' in patch && patch.base !== undefined) throw new Error('A story step is a patch with no base');
      const params: Record<string, Value> = {};
      for (const [k, v] of Object.entries(rest)) if (v !== undefined) params[k] = v as Value;
      if (patch) params.patch = 'kind' in patch ? patch : patchToNode(patch);
      return add(snapshots, 'snapshot', id, params).address as Address<'snapshot'>;
    },
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
      for (const container of [data, views, refs, snapshots]) if (container.children.length) emit(r, container);
    });
  }

  return {
    data: dataBuilder,
    views: viewsBuilder,
    refs: refsBuilder,
    snapshots: snapshotsBuilder,
    build,
    check: () => checkSeqSpec(build()),
    toString: () => '/',
  };
}
