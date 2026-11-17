/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view (_spec/README.md "Resolved view"): tier 1, over the static model only. Views, sections (an axis:
// a system and its range), tracks, layers with their params, themes as values, each layer's fallback chain and its
// data as an address or a query with the static path from the data's system to its section's. Follow links resolved
// to their paths and scales. No data, no dynamic node, no pixels, no target: the same for every renderer.

import {
  exactAddress,
  type Issue,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import { OpenLength, type Path, projectRanges, type StaticModel } from '@molstar/seqx-core/resolve/index';
import { ViewIssues } from './issues.js';

/** A layer type: what it draws and the simpler type a reader that lacks it draws instead. */
export interface LayerType {
  readonly shape: 'letters' | 'ranges' | 'values';
  readonly fallback?: string;
}

/** The first slice's layer types; none has a fallback yet, so each chain is the type alone. */
export const LayerTypes: Readonly<Record<string, LayerType>> = {
  letters: { shape: 'letters' },
  ranges: { shape: 'ranges' },
  values: { shape: 'values' },
};

/** An axis: a system and the positions shown, `[start, end, …]` pairs (two when a cyclic range wraps the origin). */
export interface Axis {
  readonly system: string;
  /** Undefined: open (its source tells). */
  readonly length?: number;
  readonly cyclic: boolean;
  readonly ranges: readonly number[];
}

/** How a section follows another: the path from the leader's system to its own, and the scale of a position. */
export interface Follow {
  readonly leader: string;
  readonly path: Path;
  /** Follower positions per leader position: 3 from residues to nucleotides, 1/3 back, else 1. */
  readonly scale: number;
}

/** What a layer draws: an address of a node, or a query; the system its data is on and the shape it reads. */
export interface LayerData {
  readonly kind: 'address' | 'query';
  readonly address: string;
  readonly system: string;
  readonly shape: 'letters' | 'ranges' | 'values';
}

export interface ResolvedLayer {
  /** A layer has an address only with an id; `index` is its place in its track either way. */
  readonly address?: string;
  readonly index: number;
  readonly type: string;
  /** Its type, then the fallbacks a reader without it draws, in order. */
  readonly fallbacks: readonly string[];
  /** Its params other than `type` and `data`, as written. */
  readonly params: Readonly<Record<string, Value>>;
  /** Themes as values (the theme language comes later): the channel's value as written. */
  readonly theme: { readonly color?: Value };
  readonly data?: LayerData;
  /** The static path from the data's system to the section's; undefined when the same system or none. */
  readonly path?: Path;
  /** Whether its data reaches its section: false with an issue. */
  readonly placed: boolean;
}

export interface ResolvedTrack {
  readonly address: string;
  readonly label: string;
  readonly height?: number;
  readonly layers: readonly ResolvedLayer[];
}

export interface ResolvedSection {
  readonly address: string;
  readonly axis: Axis;
  readonly follow?: Follow;
  readonly tracks: readonly ResolvedTrack[];
}

export interface ResolvedView {
  readonly address: string;
  readonly layout: string;
  readonly sections: readonly ResolvedSection[];
  readonly issues: readonly Issue[];
}

const cache = new WeakMap<StaticModel, Map<string, ResolvedView>>();

/** The resolved views of a model's spec, in spec order; kept per model and view. */
export function resolveViews(model: StaticModel): ResolvedView[] {
  const out: ResolvedView[] = [];
  for (const c of model.spec.children ?? [])
    if (c.kind === 'views' && c.id !== undefined)
      for (const v of c.children ?? [])
        if (v.kind === 'view' && v.id !== undefined) out.push(resolveView(model, printAddress([c.id, v.id]))!);
  return out;
}

/** A view resolved over the static model, or undefined when the spec has no view at the address. */
export function resolveView(model: StaticModel, address: string): ResolvedView | undefined {
  let views = cache.get(model);
  if (!views) cache.set(model, (views = new Map()));
  const hit = views.get(address);
  if (hit) return hit;
  const node = nodeAt(model.spec, address);
  if (node?.kind !== 'view') return undefined;
  const v = build(model, node, address);
  views.set(address, v);
  return v;
}

function nodeAt(root: Node, address: string): Node | undefined {
  const p = parseAddress(address).pattern;
  if (!p) return undefined;
  let n: Node | undefined = root;
  for (const s of p.segments) {
    if (typeof s !== 'string') return undefined;
    n = n?.children?.find((c) => c.id === s);
  }
  return n;
}

const param = (n: Node, key: string): Value | undefined => n.params?.[key];

function build(model: StaticModel, view: Node, address: string): ResolvedView {
  const issues: Issue[] = [];
  const report = (code: keyof typeof ViewIssues.codes, message: string, at: string, cause?: Issue) =>
    issues.push(
      issueOf(ViewIssues, code, message, {
        address: at,
        hint: cause ? `${cause.message}: ${cause.hint ?? ''}`.trim() : ViewIssues.codes[code].hint,
      }),
    );
  const sections: ResolvedSection[] = [];
  const systemOf = new Map<string, string>();
  const sectionNodes = (view.children ?? []).filter((c) => c.kind === 'section' && c.id !== undefined);
  for (const s of sectionNodes) systemOf.set(`${address}/${s.id}`, exactAddress(param(s, 'of')) ?? '');

  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    const axis = axisOf(model, systemOf.get(at)!, param(s, 'range'));
    let follow: Follow | undefined;
    const leader = exactAddress(param(s, 'follow'));
    if (leader !== undefined) {
      const r = model.path(systemOf.get(leader) ?? '', axis.system);
      if (r.path) follow = { leader, path: r.path, scale: r.path.step === 1 ? 3 : r.path.step === -1 ? 1 / 3 : 1 };
      else report('seq.view.follow.unplaced', 'The section cannot follow its leader', at, r.issue);
    }
    const tracks: ResolvedTrack[] = [];
    for (const t of s.children ?? []) {
      if (t.kind !== 'track' || t.id === undefined) continue;
      const ta = `${at}/${t.id}`;
      const layers: ResolvedLayer[] = [];
      (t.children ?? []).forEach((l, index) => {
        if (l.kind !== 'layer') return;
        layers.push(layerOf(model, l, index, ta, axis.system, report));
      });
      const label = param(t, 'label');
      const height = param(t, 'height');
      tracks.push({
        address: ta,
        label: typeof label === 'string' ? label : t.id,
        ...(typeof height === 'number' ? { height } : {}),
        layers,
      });
    }
    sections.push({ address: at, axis, ...(follow ? { follow } : {}), tracks });
  }
  const layout = param(view, 'layout');
  return { address, layout: typeof layout === 'string' ? layout : 'tracks', sections, issues };
}

function axisOf(model: StaticModel, system: string, range: Value | undefined): Axis {
  const s = model.systems.get(system);
  const cyclic = s?.cyclic ?? false;
  const length = s?.length;
  const end = length ?? OpenLength;
  let ranges: number[] = [0, end];
  if (Array.isArray(range) && range.length === 2) {
    const [a, b] = range as number[];
    ranges = a > b && cyclic ? [0, b, a, end] : [a, b];
  }
  return { system, ...(length !== undefined ? { length } : {}), cyclic, ranges };
}

function layerOf(
  model: StaticModel,
  l: Node,
  index: number,
  track: string,
  system: string,
  report: (code: keyof typeof ViewIssues.codes, message: string, at: string, cause?: Issue) => void,
): ResolvedLayer {
  const type = String(param(l, 'type'));
  const address = l.id === undefined ? undefined : `${track}/${l.id}`;
  const at = address ?? track;
  const fallbacks: string[] = [];
  for (let t: string | undefined = type; t !== undefined && !fallbacks.includes(t); t = LayerTypes[t]?.fallback)
    fallbacks.push(t);
  const params: Record<string, Value> = {};
  for (const [k, v] of Object.entries(l.params ?? {})) if (k !== 'type' && k !== 'data') params[k] = v;
  const color = param(l, 'color');
  const base = {
    ...(address ? { address } : {}),
    index,
    type,
    fallbacks,
    params,
    theme: color === undefined ? {} : { color },
  };

  const data = dataOf(model, param(l, 'data'), type, system);
  if (!data) {
    report('seq.view.layer.data', 'The layer’s data names nothing the static model holds', at);
    return { ...base, placed: false };
  }
  if (data.system === system) return { ...base, data, placed: true };
  const r = model.path(data.system, system);
  if (!r.path) {
    report('seq.view.layer.unplaced', 'The layer’s data does not reach its section', at, r.issue);
    return { ...base, data, placed: false };
  }
  return { ...base, data, path: r.path, placed: true };
}

/** A layer's data: the section's letters when it names none; an address's node; a query's first node. */
function dataOf(model: StaticModel, value: Value | undefined, type: string, section: string): LayerData | undefined {
  const shape = LayerTypes[type]?.shape ?? 'ranges';
  if (value === undefined)
    return type === 'letters' ? { kind: 'address', address: section, system: section, shape } : undefined;
  if (typeof value !== 'string') return undefined;
  const exact = exactAddress(value);
  let prefix = exact;
  if (prefix === undefined) {
    const p = parseAddress(value).pattern;
    const ids: string[] = [];
    for (const s of p?.segments ?? []) {
      if (typeof s !== 'string') break;
      ids.push(s);
    }
    if (!ids.length) return undefined;
    prefix = printAddress(ids);
  }
  const kind = exact === undefined ? 'query' : 'address';
  const address = exact ?? printAddress(parseAddress(value).pattern!.segments);
  if (model.systems.has(prefix)) return { kind, address, system: prefix, shape };
  const a = model.annotations.get(prefix);
  return a ? { kind, address, system: a.target, shape } : undefined;
}

/**
 * A follower's range from its leader's `[start, end)`: the image along the follow's path, its hull; on a cyclic
 * system the shortest arc that holds the image (EC-12), `start > end` when it wraps. Undefined when nothing maps.
 */
export function followRange(
  section: ResolvedSection,
  start: number,
  end: number,
): { start: number; end: number } | undefined {
  if (!section.follow) return undefined;
  const image = projectRanges(section.follow.path, [Math.floor(start), Math.ceil(end)]);
  if (!image.length) return undefined;
  const pairs: [number, number][] = [];
  for (let i = 0; i < image.length; i += 2) pairs.push([image[i], image[i + 1]]);
  pairs.sort((a, b) => a[0] - b[0]);
  const { cyclic, length } = section.axis;
  if (!cyclic || length === undefined) return { start: pairs[0][0], end: Math.max(...pairs.map((p) => p[1])) };
  // the arc is the circle less its widest gap
  let best = -1;
  let from = 0;
  let to = 0;
  let reach = pairs[0][1];
  for (let i = 1; i < pairs.length; i++) {
    const gap = pairs[i][0] - reach;
    if (gap > best) [best, from, to] = [gap, pairs[i][0], reach];
    reach = Math.max(reach, pairs[i][1]);
  }
  const around = length - reach + pairs[0][0];
  if (around >= best) return { start: pairs[0][0], end: reach };
  return { start: from, end: to };
}
