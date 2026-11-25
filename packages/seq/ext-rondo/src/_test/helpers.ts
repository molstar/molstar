/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Test helpers: a toy genome (chromosomes of a few Mb), seeded links between them (held inline or served by a `seq`
// source by region, one picture), a spec through `rondoRing`, and a frame settled on a composition with the extension
// (or without it: the fallback), its ribbons' lay at hand.

import type { SeqConfigItem, SeqExtension } from '@molstar/seqx-core/context';
import { createSeqContext } from '@molstar/seqx-core/context/default';
import { DataStore } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import {
  bundledTextMetrics,
  DarkTokens,
  type Frame,
  type FrameInput,
  LightTokens,
  openModel,
  type RenderConfigOverrides,
  Renderer,
  renderConfig,
} from '@molstar/seqx-core/render';
import { rondo } from '../index.js';
import type { RondoLay } from '../lay.js';
import { rondoRing, type RondoRingOptions } from '../ring.js';

/** A spec (the builder's tree). */
export type Node = ReturnType<ReturnType<typeof createSeqSpec>['build']>;

/** A small generator (Park–Miller): the tests' data is deterministic without the core's hashes. */
export function rng(seed: number): () => number {
  let x = seed % 2147483647 || 1;
  return () => (x = (x * 48271) % 2147483647) / 2147483647;
}

/** `n` chromosomes, `chr1` the longest, each a few Mb. */
export function toyChromosomes(n = 6): { name: string; length: number }[] {
  return Array.from({ length: n }, (_, i) => ({ name: `chr${i + 1}`, length: 8_000_000 - i * 700_000 }));
}

/** `n` chromosomes of a human-sized genome (250 Mb down to about 45 Mb): the performance cases' 23. */
export function genomeChromosomes(n = 23): { name: string; length: number }[] {
  return Array.from({ length: n }, (_, i) => ({ name: `chr${i + 1}`, length: Math.round(250e6 / (1 + i * 0.2)) }));
}

export interface Link {
  readonly a: number;
  readonly b: number;
  readonly as: number;
  readonly bs: number;
  readonly type: string;
}

/** Seeded links: `within` of them within a chromosome near each other, the rest between hub pairs and anywhere. */
export function toyLinks(chroms: readonly { length: number }[], n: number, seed = 7, within = 0.6): Link[] {
  const r = rng(seed);
  const out: Link[] = [];
  const hubs = [
    [0, 2],
    [1, 4],
    [3, 5 % chroms.length],
  ];
  for (let i = 0; i < n; i++) {
    let a = Math.floor(r() * chroms.length);
    let b = a;
    if (r() >= within) {
      if (r() < 0.6) [a, b] = hubs[Math.floor(r() * hubs.length)];
      else b = Math.floor(r() * chroms.length);
    }
    const la = chroms[a].length;
    const lb = chroms[b].length;
    const as = Math.floor(r() * (la - 1000));
    const bs = a === b ? Math.min(lb - 1000, as + Math.floor(r() ** 3 * 2_000_000)) : Math.floor(r() * (lb - 1000));
    out.push({ a, b, as, bs, type: r() < 0.7 ? 'contact' : 'other' });
  }
  return out;
}

/**
 * The links as BEDPE (`chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2 type`): read whole, by either
 * end; its `type` column mapped by the source (`fields: { type: '10' }`).
 */
export function linksBedpe(chroms: readonly { name: string }[], links: readonly Link[]): string {
  return `${links
    .map((l, i) =>
      [chroms[l.a].name, l.as, l.as + 500, chroms[l.b].name, l.bs, l.bs + 500, `l${i}`, 0, '+', '+', l.type].join('\t'),
    )
    .join('\n')}\n`;
}

export interface RingSpecOptions extends Partial<Omit<RondoRingOptions, 'assembly' | 'links'>> {
  readonly chroms?: readonly { name: string; length: number }[];
  readonly links?: readonly Link[];
  /** Held inline (`a_name`, `b_name`) rather than read by region. */
  readonly held?: boolean;
  readonly include?: readonly string[];
  readonly layout?: 'rings' | 'tracks';
  readonly extensions?: readonly { readonly id: string; readonly version: string }[];
  /** The layer's data filter, a query after the links' address (`/(except :of (item :type other))`). */
  readonly query?: string;
  /** Params set on every ribbons layer after the build (a colour, an opacity, a width), as a patch would. */
  readonly layer?: Readonly<Record<string, unknown>>;
}

/** A genome ring through `rondoRing` over the toy genome and its links. */
export function ringSpec(o: RingSpecOptions = {}): Node {
  const chroms = o.chroms ?? toyChromosomes();
  const links = o.links ?? toyLinks(chroms, 2000);
  const b = createSeqSpec({ title: 'genome ring', extensions: o.extensions ?? [rondo] });
  const g = b.data.assembly('toy', assembly('Toy', chroms));
  const pairs = o.held
    ? b.data.pairs('links', {
        a: g,
        rows: {
          a_name: links.map((l) => chroms[l.a].name),
          b_name: links.map((l) => chroms[l.b].name),
          a_start: Int32Array.from(links, (l) => l.as),
          a_end: Int32Array.from(links, (l) => l.as + 500),
          b_start: Int32Array.from(links, (l) => l.bs),
          b_end: Int32Array.from(links, (l) => l.bs + 500),
          type: links.map((l) => l.type),
          // a numeric field a width sums (deterministic, never 0)
          q: Float64Array.from(links, (_, i) => 0.001 + ((i * 37) % 100) / 1000),
        },
      } as never)
    : b.data.pairs('links', {
        a: g,
        source: b.data.source('src', {
          format: 'bedpe',
          value: linksBedpe(chroms, links) as never,
          // no lines, no columns to map
          ...(links.length ? { fields: { type: '10' } } : {}),
        } as never),
      } as never);
  b.views.view('v', (v) => {
    const { ring } = rondoRing(v, {
      assembly: g,
      links: (o.query ? `${pairs}${o.query}` : pairs) as never,
      ...(o.square !== undefined ? { square: o.square } : { square: false }),
      ...(o.hues !== undefined ? { hues: o.hues } : {}),
      ...(o.opens !== undefined ? { opens: o.opens } : {}),
      ...(o.kinds !== undefined ? { kinds: o.kinds } : {}),
      ...(o.scale !== undefined ? { scale: o.scale } : {}),
    });
    void ring;
  });
  let root = b.build();
  if (o.layer) root = withRibbons(root, o.layer);
  // the section's `include`, `layout` as patches would set them
  if (o.include || o.layout) {
    const views = root.children!.find((c) => c.kind === 'views')!;
    const view = views.children!.find((c) => c.id === 'v')!;
    const section = view.children!.find((c) => c.id === 'ring')!;
    const params = {
      ...section.params,
      ...(o.include ? { include: [...o.include] } : {}),
      ...(o.layout ? { layout: o.layout } : {}),
    };
    if (o.layout === 'tracks') delete (params as Record<string, unknown>).center;
    const s2 = { ...section, params };
    const v2 = { ...view, children: view.children!.map((c) => (c === section ? s2 : c)) };
    const views2 = { ...views, children: views.children!.map((c) => (c === view ? v2 : c)) };
    root = { ...root, children: root.children!.map((c) => (c === views ? views2 : c)) };
  }
  return root;
}

/** The spec with params set on every `rondo:ribbons` layer (as a patch would set them). */
export function withRibbons(root: Node, params: Readonly<Record<string, unknown>>): Node {
  const map = (n: Node): Node => {
    if (n.kind === 'layer' && n.params?.type === 'rondo:ribbons')
      return { ...n, params: { ...n.params, ...params } as never };
    return n.children ? { ...n, children: n.children.map(map) } : n;
  };
  return map(root);
}

/** Every `rondo:ribbons` lay of a frame, in draw order (the first drawn first). */
export function lays(f: Frame): RondoLay[] {
  const out: RondoLay[] = [];
  for (const s of f.layout.sections)
    for (const l of s.sectionLays) if (l.type.name === 'rondo:ribbons') out.push(l.lay as RondoLay);
  return out;
}

export interface Settled {
  readonly store: DataStore;
  readonly renderer: Renderer;
  frame(input: FrameInput): Promise<Frame>;
}

/** A frame settled on a composition with these extensions (`rondo` by default; none: the fallback). */
export async function settle(
  spec: Node,
  extensions: readonly SeqExtension[] = [rondo],
  config: RenderConfigOverrides = {},
  dark = false,
  /** The composition's settings (`[item, value]` pairs), read by the layer types as composed. */
  settings: readonly (readonly [SeqConfigItem<any>, unknown])[] = [],
): Promise<Settled> {
  const ctx = createSeqContext({ extensions, config: settings });
  const store = new DataStore();
  const { model, view } = await openModel(spec, store, { runtimes: [BedRuntime], kinds: ctx.kinds });
  const renderer = new Renderer({
    types: ctx.types,
    layouts: ctx.layouts,
    model,
    view: view!,
    data: store,
    tokens: dark ? DarkTokens : LightTokens,
    text: bundledTextMetrics(),
    config: ctx.renderConfig.settings
      ? { ...renderConfig(config), settings: ctx.renderConfig.settings }
      : renderConfig(config),
    ...(ctx.schemes ? { schemes: ctx.schemes } : {}),
  } as never);
  return {
    store,
    renderer,
    async frame(input) {
      let f = renderer.frame(input);
      for (let i = 0; i < 6; i++) {
        store.want(f.wants);
        await store.whenSettled();
        const next = renderer.frame(input);
        if (!next.rebuilt) break;
        f = next;
      }
      // the view settled: nothing pending (a host says so once its store has)
      return renderer.frame({ ...input, settled: input.settled ?? true });
    },
  };
}

/** The ribbons' lay of a frame (the first section's section layer of the type). */
export function lay(f: Frame): RondoLay | undefined {
  for (const s of f.layout.sections)
    for (const l of s.sectionLays) if (l.type.name === 'rondo:ribbons') return l.lay as RondoLay;
  return undefined;
}

export const Size = { width: 720, height: 760 };
