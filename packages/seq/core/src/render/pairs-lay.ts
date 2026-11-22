/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Key, type Location, location } from '@molstar/seqx-core/query';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import type { PairsData, StaticModel } from '@molstar/seqx-core/resolve';
import type { Channel, ThemeDatum } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer, ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, Extents, type LayerChannels } from './channels.js';
import type { RenderConfig } from './config.js';
import type { PassWriter } from './draw-lists.js';
import { compileRowFilter } from './filter.js';
import { layerSource } from './layer-data.js';
import { type PairEnds, pairEnds } from './pairs-ends.js';
import type { PairsRegion } from './pairs-region.js';
import { splitCrossing } from './pairs-split.js';
import { type Box, type EmphasisHighlight, type Hit } from './layer-type.js';
import type { SeqTokens } from './tokens.js';

// The pairs layer's prepared data and lay (_spec/README.md "Render", The pairs layer): a pair as the themes and the
// filter read it, the layer's data and params read once, what a layout keeps, a pair's hit.

/** An arc's least height, px: never flat. */
export const MinArc = 3;

/** A dot's least side, CSS px. */
export const DotMin = 4;

/** The feet's inset from the band's edge, px: the stroke's half and a hair, so an arc is drawn within its band. */
export const Inset = 1;

/** The points a curve is sampled at for a hit. */
export const HitSamples = 16;

/** The caret's half width at a placed end whose partner is not placed, px. */
export const Caret = 3;

export type Side = 'above' | 'below';

/** The default color: by kind (the `type` column), a categorical palette whose colors are fixed by value. */
export const DefaultColor = themeForm.field('type', { list: false });

/** A pair as the themes and the filter read it: its kind, label, key, strands and ends. */
export class PairDatum implements ThemeDatum {
  value = Number.NaN;
  letter = '';
  start = 0;
  end = 0;
  k = 0;
  key: string | number = 0;
  constructor(readonly data: PairsData) {}
  set(k: number): this {
    const d = this.data;
    this.k = k;
    this.key = d.keys[k];
    this.start = d.aStart[k];
    this.end = d.aEnd[k];
    return this;
  }
  field(name: string): string | number | undefined {
    const d = this.data;
    const k = this.k;
    switch (name) {
      case 'type':
        return d.kind[k] < 0 ? undefined : d.kinds[d.kind[k]];
      case 'label':
        return d.labels?.[k] ?? undefined;
      case 'id':
        return d.keys[k];
      case 'a_strand':
        return d.aStrand[k] < 0 ? '-' : '+';
      case 'b_strand':
        return d.bStrand[k] < 0 ? '-' : '+';
      case 'a_start':
        return d.aStart[k];
      case 'a_end':
        return d.aEnd[k];
      case 'b_start':
        return d.bStart[k];
      case 'b_end':
        return d.bEnd[k];
      default:
        return d.extra?.[name]?.[k] ?? undefined;
    }
  }
}

const PairEnds = new Set(['a_start', 'a_end', 'b_start', 'b_end']);
const PairEnvs = new WeakMap<PairsData, { numeric(n: string): boolean }>();
/** The color channel's environment: the ends and the numeric extra columns are numbers, the rest text. */
export function pairEnv(d: PairsData | undefined): { numeric(n: string): boolean } {
  let env = d && PairEnvs.get(d);
  if (env) return env;
  const numbers = new Set(
    Object.entries(d?.extra ?? {}).flatMap(([k, v]) => (v.some((x) => typeof x === 'number') ? [k] : [])),
  );
  env = { numeric: (n: string) => PairEnds.has(n) || numbers.has(n) };
  if (d) PairEnvs.set(d, env);
  return env;
}

export interface PairsPrep {
  readonly owner: string;
  readonly data?: PairsData;
  readonly aSystem: string;
  readonly bSystem: string;
  /** A complex's links: the units their ends are on (`PairsData.aUnit`, `bUnit` index them), each row's own. */
  readonly endUnits?: readonly string[];
  readonly side: Side;
  readonly reach: number;
  readonly acrossLines: boolean;
  readonly mirror: boolean;
  readonly triangle: 'both' | 'upper' | 'lower';
  readonly channels: LayerChannels;
  /** 1 where the layer's data filter keeps a pair (all when none). */
  readonly shown: Uint8Array;
  /** The pairs the filter left out (the track's "n filtered"); whether its form is not read row by row (said). */
  readonly filtered: number;
  readonly unreadFilter: boolean;
  readonly datum?: PairDatum;
  /** With `split`: 1 where a pair is drawn on the other side. */
  readonly other?: Uint8Array;
  /** Each pair's index by its key. */
  readonly index: ReadonlyMap<string | number, number>;
  /** Each pair's element and extent (both ends), made when first hit. */
  readonly elements: (Key | undefined)[];
  readonly extents: (Location | undefined)[];
  /** The resolver's steps: a block of a unit shows what it instantiates (an end on a molecule reaches its units). */
  readonly model: StaticModel;
  /** The ends drawn in a section (`pairEnds`), by section. */
  readonly units: WeakMap<ResolvedSection, PairEnds>;
  /** Read by region (a source's pairs): where, and what the frame read is (`pairs-region.ts`). */
  readonly region?: PairsRegion;
  /** Bins drawn as arcs: each one's count (its stroke and its cull's order) and the level read. */
  readonly weight?: Float32Array;
  readonly level?: number;
}

export type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

export class PairsLay {
  key = '';
  section = '';
  owner = '';
  tokens: SeqTokens | undefined;
  config: RenderConfig | undefined;
  prep: PairsPrep | undefined;
  dpr = 1;
  /** `arcs` (a section's) or `dots` (a grid's strip). */
  mode: 'arcs' | 'dots' | 'none' = 'none';
  /** The color channel as last laid out, its open domains' extents. */
  colorChannel: Channel | undefined;
  readonly extents = new Extents();
  /** Arcs: the pair's index, and its curve on the page (start, two controls, end: 8 numbers), its color. */
  n = 0;
  pair = new Int32Array(64);
  curve = new Float64Array(512);
  color = new Uint32Array(64);
  /** Each arc's stroke, CSS px (a bin pair's by its count, else 1). */
  width = new Float32Array(64);
  /** Carets at ends whose partner is not placed: x, y, 1 above or −1 below, the pair. */
  carets: number[] = [];
  /** The density strip of ends past `maxArcs`: counts a device column from `densityX`, at `densityY`. */
  density = new Uint16Array(0);
  densityN = 0;
  densityX = 0;
  densityY = 0;
  densityH = 0;
  densityMax = 0;
  /** The page's box the arcs are clipped to (the data's, for translated strips). */
  clip: Box | undefined;
  /** Dots: x, y, w, h each, the pair, its color; segments: x0, y0, x1, y1, width each, the pair, its color. */
  dots = new Float64Array(64);
  dotPair = new Int32Array(16);
  dotColor = new Uint32Array(16);
  nd = 0;
  segs = new Float64Array(80);
  segPair = new Int32Array(16);
  segColor = new Uint32Array(16);
  ns = 0;
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
  /** A path's points for a segment (two). */
  readonly seg = new Float64Array(4);
}

/** A lay's pairs lit by the highlight channels (module-private): what they were made for, and each pair's channel. */
export interface Lit {
  /** The highlight channels and the prep `lit` was made for (the viewer keeps the list until a channel changes). */
  for: readonly EmphasisHighlight[] | undefined;
  prep: PairsPrep | undefined;
  /** Each pair's highlight channel (its index + 1; 0 none): the last channel listing its key. */
  lit: Int32Array;
}

export const LitOf = new WeakMap<PairsLay, Lit>();

/** The annotation a layer's data names: its address, else (a query) the annotation it starts at. */
export function annotationOf(model: StaticModel, address: string) {
  const a = model.annotations.get(address);
  if (a) return { address, entry: a };
  for (const [at, entry] of model.annotations) if (address.startsWith(`${at}/`)) return { address: at, entry };
  return undefined;
}

export function prepare(layer: ResolvedLayer, model: StaticModel): PairsPrep {
  const named = annotationOf(model, layer.data?.address ?? '');
  const owner = named?.address ?? '';
  const a = named?.entry;
  const data = a?.shape === 'pairs' ? a.pairs : undefined;
  const p = layer.params;
  const index = new Map<string | number, number>();
  if (data) for (let k = 0; k < data.n; k++) if (!index.has(data.keys[k])) index.set(data.keys[k], k);
  const t = p.triangle;
  // the layer's data filter, read pair by pair once (_spec/README.md "Themes", Filters)
  const form = layer.data?.filter;
  const filter = form ? compileRowFilter(form) : undefined;
  const baseForm = layer.data?.base;
  const base = baseForm ? compileRowFilter(baseForm) : undefined;
  const datum = data ? new PairDatum(data) : undefined;
  const shown = new Uint8Array(data?.n ?? 0).fill(1);
  let filtered = 0;
  if ((filter || base) && datum)
    for (let k = 0; k < shown.length; k++)
      if (base && !base(datum.set(k))) shown[k] = 0;
      else if (filter && !filter(datum.set(k))) {
        shown[k] = 0;
        filtered++;
      }
  // a source's pairs, read by region on the layer's first system (`pairs-region.ts`)
  const src = !data && a?.shape === 'pairs' && a.access?.kind === 'source' ? layerSource(model, layer) : undefined;
  const bName = a?.ySystem ? model.systems.get(a.ySystem)?.name : undefined;
  const region = src ? { src, aName: src.system, bName: bName ?? src.system } : undefined;
  return {
    owner,
    ...(data ? { data } : {}),
    ...(region ? { region } : {}),
    ...(datum ? { datum } : {}),
    channels: channelForms(layer.theme),
    shown,
    filtered,
    unreadFilter: (form !== undefined && !filter) || (baseForm !== undefined && !base) || layer.data?.unread === true,
    aSystem: a?.target ?? '',
    bSystem: a?.ySystem ?? '',
    ...(a?.units ? { endUnits: a.units } : {}),
    side: p.side === 'below' ? 'below' : 'above',
    reach: typeof p.reach === 'number' && p.reach >= 0 ? p.reach : 0.25,
    acrossLines: p.across_lines === true,
    mirror: p.mirror === true,
    triangle: t === 'upper' || t === 'lower' ? t : 'both',
    ...(data && p.crossing === 'split' && a?.target === a?.ySystem ? { other: splitCrossing(data, shown) } : {}),
    index,
    elements: [],
    extents: [],
    model,
    units: new WeakMap(),
  };
}

/** Each pair's ends in a section, kept by the section and the pieces they were placed on. */
export function endsIn(prep: PairsPrep, section: ResolvedSection, layout: AxisTransform): PairEnds {
  const kept = prep.units.get(section);
  if (kept && kept.pieces === (section.axis.ranges ?? layout.pieces)) return kept;
  const ends = pairEnds({ ...prep, data: prep.data! }, section, layout);
  prep.units.set(section, ends);
  return ends;
}

/** A pair's color from the layer's channel (its kind's by default); the ranges' ink where it gives none. */
export function pairColor(lay: PairsLay, prep: PairsPrep, k: number): number {
  const c = lay.colorChannel?.at(prep.datum!.set(k));
  return typeof c === 'number' ? c : lay.tokens!.colors['seq-range'];
}

/** A pair's `a` end's system: its row's unit under a complex, else the layer's. */
export function aOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.aUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.aSystem) : prep.aSystem;
}

/** A pair's `b` end's system: its row's unit under a complex, else the layer's. */
export function bOf(prep: PairsPrep, k: number): string {
  const u = prep.endUnits && prep.data!.bUnit;
  return u ? (prep.endUnits![u[k]] ?? prep.bSystem) : prep.bSystem;
}

/** A pair's hit: its element and its extent, both ends (made once a pair). */
export function pairHit(lay: PairsLay, k: number): Hit {
  const prep = lay.prep!;
  const d = prep.data!;
  let element = prep.elements[k];
  if (!element) element = prep.elements[k] = { owner: prep.owner, key: d.keys[k] };
  let extent = prep.extents[k];
  const as = aOf(prep, k);
  const bs = bOf(prep, k);
  if (!extent) {
    const parts =
      as === bs
        ? [{ system: as, ranges: [d.aStart[k], d.aEnd[k], d.bStart[k], d.bEnd[k]] }]
        : [
            { system: as, ranges: [d.aStart[k], d.aEnd[k]] },
            { system: bs, ranges: [d.bStart[k], d.bEnd[k]] },
          ];
    extent = prep.extents[k] = location({ parts });
  }
  const h = lay.scratch;
  h.layer = lay.key;
  h.section = lay.section;
  h.system = as;
  h.element = element;
  h.extent = extent;
  h.position = undefined;
  h.ySystem = undefined;
  h.yPosition = undefined;
  h.target = undefined;
  h.bin = prep.level;
  return h;
}

export function clipTo(out: PassWriter, b: Box) {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

/** The squared distance from (x, y) to the segment (x0, y0)–(x1, y1) (small: inlined in the hit loops). */
export function segmentDistance2(x0: number, y0: number, x1: number, y1: number, x: number, y: number): number {
  const dx = x1 - x0;
  const dy = y1 - y0;
  const l2 = dx * dx + dy * dy;
  let f = l2 > 0 ? ((x - x0) * dx + (y - y0) * dy) / l2 : 0;
  f = f < 0 ? 0 : f > 1 ? 1 : f;
  const ex = x0 + f * dx - x;
  const ey = y0 + f * dy - y;
  return ex * ex + ey * ey;
}
