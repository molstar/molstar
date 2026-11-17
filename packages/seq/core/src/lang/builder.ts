/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The language builder (_spec/README.md, Language builder): the way most specs are written. It reads as the spec does
// (`b.data.sequence(...)`, `b.views.view(...)`), every call returns the address of what it adds, so a handle stands
// wherever an address does, and `build` gives the normalized tree the generic engine of `@molstar/core/repr/builder`
// builds. Calls are recorded as they are written and replayed on the engine, so misuse that types do not reach (an id
// twice, a malformed id) is thrown where it is written.

import { createBuilder, type ParamsFor } from '@molstar/core/repr/builder';
import { IdHint, IdPattern, type Issue, type Node, printAddress } from '@molstar/core/repr/tree';
import { checkSeqSpec } from './check.js';
import { type SeqKinds, SeqLanguage, SeqLanguageVersion } from './language.js';

type Params<C extends keyof SeqKinds & string> = ParamsFor<SeqKinds, C>;

/** What stands where an address does: the address, or a builder that has one. */
export interface Addressed {
  readonly address: string;
}

type Followed = Omit<Params<'section'>, 'follow'> & { readonly follow?: string | Addressed };

export interface DataBuilder {
  source(id: string, params: Params<'source'>): string;
  sequence(id: string, params: Params<'sequence'>): string;
  ranges(id: string, params: Params<'ranges'>): string;
  values(id: string, params: Params<'values'>): string;
  map(id: string, params: Params<'map'>): string;
}

export interface ViewsBuilder {
  view(id: string, body?: (view: ViewBuilder) => void): string;
  view(id: string, params: Params<'view'>, body?: (view: ViewBuilder) => void): string;
}

export interface ViewBuilder extends Addressed {
  section(id: string, params: Followed, body?: (section: SectionBuilder) => void): SectionBuilder;
}

export interface SectionBuilder extends Addressed {
  track(id: string, params?: Params<'track'>, body?: (track: TrackBuilder) => void): TrackBuilder;
  toString(): string;
}

export interface TrackBuilder extends Addressed {
  /** A layer; with an id it has an address and is returned as a builder, without one it is anonymous. */
  layer(params: Params<'layer'>): TrackBuilder;
  layer(id: string, params: Params<'layer'>): TrackBuilder;
  toString(): string;
}

export interface SeqSpecBuilder {
  readonly data: DataBuilder;
  readonly views: ViewsBuilder;
  /** The spec as a normalized tree (print it with `printSexpr` from `@molstar/core/repr/codecs/sexpr`). */
  build(): Node;
  /** The issues of the spec so far: the description's, then the language's. */
  check(): Issue[];
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

function refs(
  params: Record<string, unknown> | undefined,
  keys: readonly string[],
): Record<string, unknown> | undefined {
  if (!params) return params;
  const out = { ...params };
  for (const k of keys) {
    const v = out[k];
    if (v && typeof v === 'object' && 'address' in v) out[k] = (v as Addressed).address;
  }
  return out;
}

/** A spec being written: `build()` is its normalized tree. */
export function createSeqSpec(options: SeqSpecOptions = {}): SeqSpecBuilder {
  const engine = createBuilder(SeqLanguage);
  const data: Rec = { kind: 'data', id: 'data', children: [], ids: new Set(), address: '/data' };
  const views: Rec = { kind: 'views', id: 'views', children: [], ids: new Set(), address: '/views' };

  /** A node added under `parent`: its id checked where it is written, its address the parent's and the id. */
  function add(parent: Rec, kind: string, id: string | undefined, params: unknown): Rec {
    if (id !== undefined) {
      if (!IdPattern.test(id)) throw new Error(`The id of ${kind} is not an id: ${IdHint}`);
      if (parent.ids.has(id)) throw new Error(`The id ${id} is used twice under one ${parent.kind}`);
      parent.ids.add(id);
    }
    const address = printAddress([...parent.address.split('/').filter(Boolean), ...(id !== undefined ? [id] : [])]);
    const rec: Rec = {
      kind,
      id,
      params: params as Record<string, unknown> | undefined,
      children: [],
      ids: new Set(),
      address,
    };
    parent.children.push(rec);
    return rec;
  }

  const dataKind = (kind: string) => (id: string, params: unknown) => add(data, kind, id, params).address;
  const dataBuilder: DataBuilder = {
    source: dataKind('source'),
    sequence: dataKind('sequence'),
    ranges: dataKind('ranges'),
    values: dataKind('values'),
    map: dataKind('map'),
  };

  function trackBuilder(rec: Rec): TrackBuilder {
    const self: TrackBuilder = {
      address: rec.address,
      layer(...argv: [Params<'layer'>] | [string, Params<'layer'>]) {
        const [id, params] = typeof argv[0] === 'string' ? [argv[0], argv[1]] : [undefined, argv[0]];
        add(rec, 'layer', id, refs(params, ['data']));
        return self;
      },
      toString: () => rec.address,
    };
    return self;
  }
  function sectionBuilder(rec: Rec): SectionBuilder {
    return {
      address: rec.address,
      track(id, params, body) {
        const t = trackBuilder(add(rec, 'track', id, params));
        body?.(t);
        return t as TrackBuilder;
      },
      toString: () => rec.address,
    };
  }
  function viewBuilder(rec: Rec): ViewBuilder {
    return {
      address: rec.address,
      section(id, params, body) {
        const s = sectionBuilder(add(rec, 'section', id, refs(params as Record<string, unknown>, ['follow'])));
        body?.(s);
        return s;
      },
    };
  }
  const viewsBuilder: ViewsBuilder = {
    view(id: string, a?: unknown, b?: unknown) {
      const [params, body] = typeof a === 'function' ? [undefined, a] : [a, b];
      const rec = add(views, 'view', id, params);
      (body as ((v: ViewBuilder) => void) | undefined)?.(viewBuilder(rec));
      return rec.address;
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
      for (const container of [data, views]) if (container.children.length) emit(r, container);
    });
  }

  return {
    data: dataBuilder,
    views: viewsBuilder,
    build,
    check: () => checkSeqSpec(build()),
  };
}
