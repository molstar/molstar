/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The resolved view (_spec/README.md "Resolved view"): tier 1, over the static model only. Views, sections (an axis:
// a system and its range), tracks, layers with their params, themes as values, each layer's fallback chain and its
// data as an address or a query with the static path from the data's system to its section's. Follow links resolved
// to their paths, directions and scales, and every path between the view's sections precomposed. No data, no dynamic
// node, no pixels, no target: the same for every renderer.

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import {
  exactAddress,
  type Issue,
  issueOf,
  type Node,
  parseAddress,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import {
  identityPath,
  type Path,
  pathDir,
  projectPoint,
  projectRanges,
  type StaticModel,
} from '@molstar/seqx-core/resolve/index';
import { ViewIssues } from './issues.js';

type Shape = 'letters' | 'ranges' | 'values';

/** A layer type: what it draws and the simpler type a reader that lacks it draws instead. */
export interface LayerType {
  readonly shape: Shape;
  readonly fallback?: string;
}

/** The first slice's layer types; none has a fallback yet, so each chain is the type alone. */
export const LayerTypes: Readonly<Record<string, LayerType>> = {
  letters: { shape: 'letters' },
  ranges: { shape: 'ranges' },
  values: { shape: 'values' },
};

/**
 * An axis: a system and the positions shown, `[start, end, …]` pairs in display order: one piece, or two when a cyclic
 * range wraps the origin (`[8, 10, 0, 2]`). `ranges` is undefined while the system's length is unknown and the
 * section shows all of it, or wraps: the layout waits for the length (pending, never a guess).
 */
export interface Axis {
  readonly system: string;
  readonly length?: number;
  readonly cyclic: boolean;
  readonly ranges?: readonly number[];
}

/** How a section follows another: the path from the leader's system to its own, its direction and scale. */
export interface Follow {
  readonly leader: string;
  /** The leader's axis: a leader range that wraps its origin is split on it. */
  readonly leaderAxis: Axis;
  readonly path: Path;
  /** Follower positions per leader position: 3 from residues to nucleotides, 1/3 back, else 1. */
  readonly scale: number;
  /** −1 where the follower runs against the leader (a minus-strand gene), else 1. */
  readonly dir: 1 | -1;
  /**
   * The continuous image of a leader position (a viewport's fractional edge): where element `floor(x)` starts plus
   * the fraction times the scale in the run's direction; undefined where it maps nowhere.
   */
  point(x: number): number | undefined;
}

/** What a layer draws: an address of a node, or a query; the system its data is on and the shape it reads. */
export interface LayerData {
  readonly kind: 'address' | 'query';
  readonly address: string;
  readonly system: string;
  /** The shape of the node it names (a query on a sequence: undefined, a location of any shape). */
  readonly shape?: Shape;
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
  /** Whether it can be drawn: its type known, its data fitting it and reaching its section; false with an issue. */
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

/** A path between two sections of a view, precomposed: hover and selection are projected along it. */
export interface ViewLink {
  readonly from: string;
  readonly to: string;
  readonly path: Path;
}

export interface ResolvedView {
  readonly address: string;
  readonly layout: string;
  readonly sections: readonly ResolvedSection[];
  /** Every ordered pair of its sections a static path joins (the same system: the identity). */
  readonly links: readonly ViewLink[];
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

type Report = (code: keyof typeof ViewIssues.codes, message: string, at: string, cause?: Issue) => void;

function build(model: StaticModel, view: Node, address: string): ResolvedView {
  const issues: Issue[] = [];
  const report: Report = (code, message, at, cause) =>
    issues.push(
      issueOf(ViewIssues, code, message, {
        address: at,
        hint: cause ? `${cause.message}: ${cause.hint ?? ''}`.trim() : ViewIssues.codes[code].hint,
      }),
    );
  const sectionNodes = (view.children ?? []).filter((c) => c.kind === 'section' && c.id !== undefined);
  const axes = new Map<string, Axis>();
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    const system = exactAddress(param(s, 'of'));
    if (system === undefined || !model.systems.has(system))
      report('seq.view.section.system', 'The section shows no system of the spec', at);
    axes.set(at, axisOf(model, system ?? '', param(s, 'range')));
  }

  const sections: ResolvedSection[] = [];
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    const axis = axes.get(at)!;
    const known = model.systems.has(axis.system);
    let follow: Follow | undefined;
    const leader = exactAddress(param(s, 'follow'));
    const leaderAxis = leader === undefined ? undefined : axes.get(leader);
    if (leader !== undefined && leaderAxis && known && model.systems.has(leaderAxis.system)) {
      const r = model.path(leaderAxis.system, axis.system);
      if (r.path) follow = followOf(leader, leaderAxis, r.path);
      else report('seq.view.follow.unplaced', 'The section cannot follow its leader', at, r.issue);
    }
    const tracks: ResolvedTrack[] = [];
    for (const t of s.children ?? []) {
      if (t.kind !== 'track' || t.id === undefined) continue;
      const ta = `${at}/${t.id}`;
      const layers: ResolvedLayer[] = [];
      (t.children ?? []).forEach((l, index) => {
        if (l.kind === 'layer') layers.push(layerOf(model, l, index, ta, known ? axis.system : undefined, report));
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

  // every path between the view's sections, composed now so hover never composes
  const links: ViewLink[] = [];
  for (const a of sections)
    for (const b of sections) {
      if (a === b || !model.systems.has(a.axis.system) || !model.systems.has(b.axis.system)) continue;
      const path =
        a.axis.system === b.axis.system ? identityPath(a.axis.system) : model.path(a.axis.system, b.axis.system).path;
      if (path) links.push({ from: a.address, to: b.address, path });
    }
  const layout = param(view, 'layout');
  return { address, layout: typeof layout === 'string' ? layout : 'tracks', sections, links, issues };
}

function followOf(leader: string, leaderAxis: Axis, path: Path): Follow {
  return {
    leader,
    leaderAxis,
    path,
    scale: path.step === 1 ? 3 : path.step === -1 ? 1 / 3 : 1,
    dir: pathDir(path),
    point: (x) => projectPoint(path, x),
  };
}

function axisOf(model: StaticModel, system: string, range: Value | undefined): Axis {
  const s = model.systems.get(system);
  const cyclic = s?.cyclic ?? false;
  const length = s?.length;
  const base = { system, ...(length !== undefined ? { length } : {}), cyclic };
  if (Array.isArray(range) && range.length === 2) {
    const [a, b] = range as number[];
    if (a <= b || !cyclic) return { ...base, ranges: [a, b] };
    // display order: from the start to the origin, then on
    return length === undefined ? base : { ...base, ranges: [a, length, 0, b] };
  }
  return length === undefined ? base : { ...base, ranges: [0, length] };
}

function layerOf(
  model: StaticModel,
  l: Node,
  index: number,
  track: string,
  system: string | undefined,
  report: Report,
): ResolvedLayer {
  const type = String(param(l, 'type'));
  const address = l.id === undefined ? undefined : `${track}/${l.id}`;
  const at = address ?? `${track}`;
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
  const layerType = LayerTypes[type];
  if (!layerType) {
    report('seq.view.layer.type', 'The layer’s type is not one this version draws', at);
    return { ...base, placed: false };
  }
  // a section with no system was reported; its layers are not placed
  if (system === undefined) return { ...base, placed: false };

  const data = dataOf(model, param(l, 'data'), type, system);
  if (!data) {
    report('seq.view.layer.data', 'The layer’s data names nothing the static model holds', at);
    return { ...base, placed: false };
  }
  if (data.shape !== undefined && data.shape !== layerType.shape) {
    report('seq.view.layer.shape', 'The layer’s type does not draw the shape of its data', at);
    return { ...base, data, placed: false };
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
  if (value === undefined)
    return type === 'letters' ? { kind: 'address', address: section, system: section, shape: 'letters' } : undefined;
  if (typeof value !== 'string') return undefined;
  const parsed = parseAddress(value).pattern;
  if (!parsed) return undefined;
  const exact = exactAddress(value);
  const ids: string[] = [];
  for (const s of parsed.segments) {
    if (typeof s !== 'string') break;
    ids.push(s);
  }
  if (!ids.length) return undefined;
  const prefix = printAddress(ids);
  const kind = exact === undefined ? 'query' : 'address';
  const address = exact ?? printAddress(parsed.segments);
  if (model.systems.has(prefix))
    return { kind, address, system: prefix, ...(kind === 'address' ? { shape: 'letters' as const } : {}) };
  const a = model.annotations.get(prefix);
  return a ? { kind, address, system: a.target, shape: a.shape } : undefined;
}

/**
 * A follower's range from its leader's, both as `[start, end, …]` pieces in display order (an axis's shape). A leader
 * piece that wraps its cyclic origin (`start > end`) is split on it first. The image's hull, or on a cyclic follower
 * the shortest arc that holds it (EC-12), as two pieces when it wraps; `arc: 'hull'` takes the hull there too.
 * Undefined when nothing maps.
 */
export function followRange(
  section: ResolvedSection,
  leader: readonly number[],
  options: { readonly arc?: 'shortest' | 'hull' } = {},
): number[] | undefined {
  const follow = section.follow;
  if (!follow) return undefined;
  const pieces: number[] = [];
  const ll = follow.leaderAxis.length;
  for (let i = 0; i + 1 < leader.length; i += 2) {
    const s = Math.floor(leader[i]);
    const e = Math.ceil(leader[i + 1]);
    if (s > e && follow.leaderAxis.cyclic && ll !== undefined) pieces.push(s, ll, 0, e);
    else pieces.push(s, e);
  }
  const image = HalfOpenRanges.normalize(projectRanges(follow.path, pieces));
  const n = image.length;
  if (!n) return undefined;
  const { cyclic, length } = section.axis;
  if (!cyclic || length === undefined || options.arc === 'hull') return [image[0], image[n - 1]];
  // the arc is the circle less its widest gap: around the origin, or between two pieces
  let best = length - image[n - 1] + image[0];
  let arc = [image[0], image[n - 1]];
  for (let i = 2; i < n; i += 2) {
    const gap = image[i] - image[i - 1];
    if (gap > best) {
      best = gap;
      arc = [image[i], length, 0, image[i - 1]];
    }
  }
  return arc;
}
