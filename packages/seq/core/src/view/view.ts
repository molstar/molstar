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
  type AccessPlan,
  codingRanges,
  identityPath,
  inversePath,
  type MapEntry,
  type Path,
  pathDir,
  pathDirection,
  projectEnd,
  projectPoint,
  projectRanges,
  type StaticModel,
} from '@molstar/seqx-core/resolve';
import { followRef } from '@molstar/seqx-core/lang';
import { type AxisShape, nearestOffsetOf, positionAt, wholePieces } from './axis.js';
import { ViewIssues } from './issues.js';

type Shape = 'letters' | 'ranges' | 'values';

/** A layer type: what it draws and the simpler type a reader that lacks it draws instead. */
export interface LayerType {
  readonly shape: Shape;
  readonly fallback?: string;
}

/**
 * The layer types: `secondary` and `variants` draw `ranges` where a reader lacks them; the `ruler` has none (omitted)
 * and, like `letters`, reads a sequence (the section's own when it names none).
 */
export const LayerTypes: Readonly<Record<string, LayerType>> = {
  letters: { shape: 'letters' },
  ruler: { shape: 'letters' },
  ranges: { shape: 'ranges' },
  values: { shape: 'values' },
  secondary: { shape: 'ranges', fallback: 'ranges' },
  variants: { shape: 'ranges', fallback: 'ranges' },
};

/**
 * An axis: a system, its direction and the positions shown, `[start, end, …]` pieces (start below end) in display
 * order: the whole system or its range (two pieces when a cyclic range wraps the origin, `[8, 10, 0, 2]`), or the
 * section's `pieces` (with `pad`, within `range`), the positions between two pieces hidden as a break (`axisBreaks`).
 * On a reversed axis the last piece is first, each shown from its end down (EC-32). `ranges` is undefined while the
 * system's length is unknown and the section shows all of it, or wraps: the layout waits for the length (pending,
 * never a guess).
 */
export interface Axis {
  readonly system: string;
  readonly length?: number;
  readonly cyclic: boolean;
  readonly ranges?: readonly number[];
  /** 1: positions count up left to right; −1: down (the section's `direction`, resolved). */
  readonly dir: 1 | -1;
  /** The form of the pieces shown: `whole`, `coding`, or the refs (comma-joined) as written; what an export says. */
  readonly pieces: string;
  /** What the direction reads along, when the section spells `{ along }`: a system's or a map's address. */
  readonly along?: string;
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
  /** A dynamic node's access plan and its system's length (a template's instance: the overlay's, not the model's). */
  readonly access?: AccessPlan;
  readonly length?: number;
}

export interface ResolvedLayer {
  /**
   * Its track's address and its id: every layer of a valid spec has one (_spec/README.md, Language), so undefined only
   * on a tree read without the language's checks; `index` is its place in its track.
   */
  readonly address?: string;
  readonly index: number;
  readonly type: string;
  /** Its type, then the fallbacks a reader without it draws, in order. */
  readonly fallbacks: readonly string[];
  /** Its params other than `type` and `data`, as written. */
  readonly params: Readonly<Record<string, Value>>;
  /**
   * Its channels' theme forms as written (_spec/README.md "Resolved view"; "Themes"): a slot's body is read by the
   * render's channels, renderer-independent (another reader evaluates it with `@molstar/seqx-core/theme`).
   */
  readonly theme: Readonly<Partial<Record<(typeof ThemeChannelNames)[number], Value>>>;
  readonly data?: LayerData;
  /** The static path from the data's system to the section's; undefined when the same system or none. */
  readonly path?: Path;
  /** Whether it can be drawn: its type known, its data fitting it and reaching its section; false with an issue. */
  readonly placed: boolean;
}

/**
 * A track template (`for_each`): a track for each dynamic node of `kind` that `source` (or a source nested `depth`
 * levels below it) holds, in `order`, filtered; expanded at layout, never here (the resolved view stays pure).
 */
export interface TrackTemplate {
  readonly source: string;
  readonly kind: string;
  readonly depth: number;
  /** `key`, or a param of the nodes; descending with a leading `-`. */
  readonly order: string;
  /** Params a node must have, each the value given (compared as written). */
  readonly filter: Readonly<Record<string, string>>;
}

export interface ResolvedTrack {
  readonly address: string;
  readonly label: string;
  readonly height?: number;
  /** Its layers; a template's layers without data are not placed here: each instance places them on its node. */
  readonly layers: readonly ResolvedLayer[];
  readonly template?: TrackTemplate;
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

/** The channel slots a layer has (`color`, `ink`, `size`, `opacity`, `line`, `arrow`, `shape`). */
export const ThemeChannelNames = ['color', 'ink', 'size', 'opacity', 'line', 'arrow', 'shape'] as const;

export interface ResolvedView {
  readonly address: string;
  readonly layout: string;
  /** The spec's scale nodes by address (`/views/<id>`): their params as written, a theme's `scale` names one. */
  readonly scales: Readonly<Record<string, Readonly<Record<string, Value>>>>;
  /** The view's export params (`export`: `legend`), as written. */
  readonly export?: Readonly<Record<string, Value>>;
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
  const byAddress = new Map<string, Node>();
  const systemOf = new Map<string, string>();
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    byAddress.set(at, s);
    const system = exactAddress(param(s, 'of'));
    if (system === undefined || !model.systems.has(system))
      report('seq.view.section.system', 'The section shows no system of the spec', at);
    systemOf.set(at, system ?? '');
  }
  // each follow's path, before the axes: a follower's direction and coding pieces read it
  const follows = new Map<string, { leader: string; path: Path }>();
  for (const [at, s] of byAddress) {
    const leader = exactAddress(param(s, 'follow'));
    const system = systemOf.get(at)!;
    const ls = leader === undefined ? undefined : systemOf.get(leader);
    if (leader === undefined || ls === undefined || !model.systems.has(system) || !model.systems.has(ls)) continue;
    const r = model.path(ls, system);
    if (r.path) follows.set(at, { leader, path: r.path });
    else report('seq.view.follow.unplaced', 'The section cannot follow its leader', at, r.issue);
  }
  const axes = new Map<string, Axis>();
  const dirs = new Map<string, 1 | -1>();
  const along = new Map<string, Path>();
  const dirOf = (at: string, depth = 0): 1 | -1 => {
    const known = dirs.get(at);
    if (known !== undefined) return known;
    const d = resolveDirection(model, at, byAddress.get(at)!, systemOf.get(at)!, follows.get(at), along, report, (l) =>
      depth < byAddress.size ? dirOf(l, depth + 1) : 1,
    );
    dirs.set(at, d);
    return d;
  };
  for (const [at, s] of byAddress) {
    const dir = dirOf(at);
    const axis = axisOf(model, s, at, systemOf.get(at)!, dir, follows.get(at)?.path, along.get(at), report);
    const d = param(s, 'direction');
    const ref = d && typeof d === 'object' && !Array.isArray(d) ? exactAddress((d as Record<string, Value>).along) : '';
    axes.set(at, ref ? { ...axis, along: ref } : axis);
  }

  const sections: ResolvedSection[] = [];
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    const axis = axes.get(at)!;
    const known = model.systems.has(axis.system);
    let follow: Follow | undefined;
    const f = follows.get(at);
    const leaderAxis = f && axes.get(f.leader);
    if (f && leaderAxis) follow = followOf(f.leader, leaderAxis, f.path);
    const tracks: ResolvedTrack[] = [];
    for (const t of s.children ?? []) {
      if (t.kind !== 'track' || t.id === undefined) continue;
      // a hidden track is kept in the spec and left out of the view: its row collapses
      if (param(t, 'hidden') === true) continue;
      const ta = `${at}/${t.id}`;
      const layers: ResolvedLayer[] = [];
      const template = templateOf(param(t, 'for_each'));
      (t.children ?? []).forEach((l, index) => {
        if (l.kind !== 'layer') return;
        const system = known ? axis.system : undefined;
        layers.push(layerOf(model, l, index, ta, system, report, template !== undefined));
      });
      const label = param(t, 'label');
      const height = param(t, 'height');
      tracks.push({
        address: ta,
        label: typeof label === 'string' ? label : t.id,
        ...(typeof height === 'number' ? { height } : {}),
        layers,
        ...(template ? { template } : {}),
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
  const scales: Record<string, Readonly<Record<string, Value>>> = {};
  for (const c of model.spec.children ?? [])
    if (c.kind === 'views' && c.id !== undefined)
      for (const n of c.children ?? [])
        if (n.kind === 'scale' && n.id !== undefined) scales[`/${c.id}/${n.id}`] = n.params ?? {};
  const ex = param(view, 'export');
  return {
    address,
    layout: typeof layout === 'string' ? layout : 'tracks',
    scales,
    ...(ex && typeof ex === 'object' && !Array.isArray(ex) ? { export: ex as Record<string, Value> } : {}),
    sections,
    links,
    issues,
  };
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

/**
 * A section's direction (_spec/README.md "Resolved view", Axis): its `direction` as written, `follow` for a follower
 * and `genomic` otherwise when absent, never read off its layers. `follow` is the leader's direction times the follow
 * path's; `along` the direction of the path from that system (or a map's own) to the section's. A path both ways takes
 * the direction of most of its length (`seq.view.direction.mixed`); one that cannot be resolved is `genomic`
 * (`seq.view.direction.unresolved`). An `along` path is left in `along` for the coding pieces.
 */
function resolveDirection(
  model: StaticModel,
  at: string,
  node: Node,
  system: string,
  follow: { leader: string; path: Path } | undefined,
  along: Map<string, Path>,
  report: Report,
  leaderDir: (leader: string) => 1 | -1,
): 1 | -1 {
  const written = param(node, 'direction');
  const value = written ?? (param(node, 'follow') !== undefined ? 'follow' : 'genomic');
  const of = (p: Path) => {
    const d = pathDirection(p);
    if (d.mixed) report('seq.view.direction.mixed', 'The path the direction is read along runs both ways', `${at}`);
    return d.dir;
  };
  if (value === 'genomic') return 1;
  if (value === 'reversed') return -1;
  if (value === 'follow') {
    if (follow) return (leaderDir(follow.leader) * of(follow.path)) as 1 | -1;
    // a follow that is unplaced says so itself; a `follow` written with no leader is unresolved
    if (written !== undefined || param(node, 'follow') === undefined)
      report('seq.view.direction.unresolved', 'The section follows no placed leader', at);
    return 1;
  }
  const ref =
    value && typeof value === 'object' && !Array.isArray(value)
      ? exactAddress((value as Record<string, Value>).along)
      : '';
  const path = ref ? pathAlong(model, ref, system) : undefined;
  if (!path) {
    report('seq.view.direction.unresolved', 'No path reaches the section from what it reads along', at);
    return 1;
  }
  along.set(at, path);
  return path.edges.length || path.step !== 0 ? of(path) : 1;
}

/** The path from what a section reads along (a system, or a map: its own linear part) to the section's system. */
/** A map's own path, from its `from` to its `to`. */
function mapPath(map: MapEntry): Path {
  return { from: map.from, to: map.to, step: map.step, linear: map.linear, edges: [] };
}

function pathAlong(model: StaticModel, ref: string, system: string): Path | undefined {
  const map = model.maps.get(ref);
  if (map) {
    const forward = mapPath(map);
    if (map.to === system) return forward;
    if (map.from === system) return inversePath(forward);
    return model.path(map.from, system).path;
  }
  if (!model.systems.has(ref) || !model.systems.has(system)) return undefined;
  return ref === system ? identityPath(system) : model.path(ref, system).path;
}

function axisOf(
  model: StaticModel,
  node: Node,
  at: string,
  system: string,
  dir: 1 | -1,
  followPath: Path | undefined,
  alongPath: Path | undefined,
  report: Report,
): Axis {
  const s = model.systems.get(system);
  const cyclic = s?.cyclic ?? false;
  const length = s?.length;
  const written = param(node, 'pieces');
  const range = param(node, 'range');
  const r = Array.isArray(range) && range.length === 2 ? (range as number[]) : undefined;
  const pieces = sectionPieces(model, written, system, followPath, alongPath);
  const form = written === undefined ? 'whole' : Array.isArray(written) ? written.join(',') : String(written);
  if (pieces === 'none') report('seq.view.section.pieces', 'The section’s pieces reach no piece of its system', at);
  if (pieces === 'ambiguous')
    report(
      'seq.view.section.pieces',
      'The section’s coding pieces are ambiguous: more than one codon map reaches its system (name one)',
      at,
    );
  const base = { system, ...(length !== undefined ? { length } : {}), cyclic, dir };
  const wholeAxis = (): Axis => {
    const whole = { ...base, pieces: 'whole' };
    // a range that wraps the origin waits for the length, as the whole system does
    if (r && (r[0] <= r[1] || !cyclic)) return { ...whole, ranges: wholePieces(0, r, dir, cyclic) };
    return length === undefined ? whole : { ...whole, ranges: wholePieces(length, r, dir, cyclic) };
  };
  if (!Array.isArray(pieces)) return wholeAxis();
  const pad = param(node, 'pad');
  const padded = typeof pad === 'number' && pad > 0 ? padPieces(pieces, pad, length) : pieces;
  // within the range, in its display order (forward), then turned for a reversed axis
  const within: number[] = [];
  const frame = r ? (r[0] <= r[1] || !cyclic ? r : [r[0], length ?? Infinity, 0, r[1]]) : [0, length ?? Infinity];
  for (let i = 0; i + 1 < frame.length; i += 2)
    for (let j = 0; j + 1 < padded.length; j += 2) {
      const a = Math.max(frame[i], padded[j]);
      const b = Math.min(frame[i + 1], padded[j + 1]);
      if (b > a) within.push(a, b);
    }
  // pieces that miss the range show the range: an axis of no pieces would wait for ever
  if (!within.length) {
    report('seq.view.section.pieces', 'The section’s pieces are all outside its range', at);
    return wholeAxis();
  }
  return { ...base, pieces: form, ranges: dir === -1 ? reversePieces(within) : within };
}

/**
 * A section's `pieces` as genomic pieces (sorted, merged): `whole` (none here: the system or its range), `coding` (the
 * follow path's coding runs, else those of the path it reads along, else the one codon map onto its system:
 * `ambiguous` when there are more), a ranges node's items on its system or a codon map's coding runs on it, or a list of
 * them (their union). `none` when what is named reaches no piece.
 */
function sectionPieces(
  model: StaticModel,
  written: Value | undefined,
  system: string,
  followPath: Path | undefined,
  alongPath: Path | undefined,
): number[] | 'whole' | 'none' | 'ambiguous' {
  if (written === undefined || written === 'whole') return 'whole';
  const out: number[] = [];
  if (written === 'coding') {
    let p = followPath?.step ? followPath : alongPath?.step ? alongPath : undefined;
    if (!p) {
      // the writer asked for the coding pieces: the one codon map onto the system names them
      let only: MapEntry | undefined;
      for (const m of model.maps.values()) {
        if (m.kind !== 'codon' || m.to !== system) continue;
        if (only) return 'ambiguous';
        only = m;
      }
      if (only) p = mapPath(only);
    }
    const c = p && codingRanges(p);
    return c?.length ? c : 'none';
  }
  for (const ref of Array.isArray(written) ? written : [written]) {
    const a = exactAddress(ref);
    if (a === undefined) return 'none';
    const map = model.maps.get(a);
    const ann = model.annotations.get(a);
    if (map?.kind === 'codon' && map.to === system) {
      const c = codingRanges(mapPath(map));
      if (!c) return 'none';
      out.push(...c);
    } else if (ann?.shape === 'ranges' && ann.target === system && ann.items) {
      for (const item of ann.items) if (item.end > item.start) out.push(item.start, item.end);
    } else return 'none';
  }
  const merged = Array.from(HalfOpenRanges.normalize(out));
  return merged.length ? merged : 'none';
}

/** Pieces widened by `pad` each side within the system, merged where they overlap. */
function padPieces(pieces: readonly number[], pad: number, length: number | undefined): number[] {
  const out: number[] = [];
  for (let i = 0; i + 1 < pieces.length; i += 2)
    out.push(Math.max(0, pieces[i] - pad), Math.min(length ?? Infinity, pieces[i + 1] + pad));
  return Array.from(HalfOpenRanges.normalize(out));
}

/** Pieces in the other display order: the last first, each still start below end. */
function reversePieces(pieces: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = pieces.length - 2; i >= 0; i -= 2) out.push(pieces[i], pieces[i + 1]);
  return out;
}

function templateOf(v: Value | undefined): TrackTemplate | undefined {
  if (!v || typeof v !== 'object' || Array.isArray(v)) return undefined;
  const r = v as Record<string, Value>;
  const source = exactAddress(r.source);
  if (source === undefined || typeof r.node_kind !== 'string') return undefined;
  const filter: Record<string, string> = {};
  if (r.filter && typeof r.filter === 'object' && !Array.isArray(r.filter))
    for (const [k, x] of Object.entries(r.filter as Record<string, Value>)) if (typeof x === 'string') filter[k] = x;
  return {
    source,
    kind: r.node_kind,
    depth: typeof r.depth === 'number' ? r.depth : 0,
    order: typeof r.order === 'string' ? r.order : 'key',
    filter,
  };
}

function layerOf(
  model: StaticModel,
  l: Node,
  index: number,
  track: string,
  system: string | undefined,
  report: Report,
  template = false,
): ResolvedLayer {
  const type = String(param(l, 'type'));
  const address = l.id === undefined ? undefined : `${track}/${l.id}`;
  const at = address ?? `${track}`;
  const fallbacks: string[] = [];
  for (let t: string | undefined = type; t !== undefined && !fallbacks.includes(t); t = LayerTypes[t]?.fallback)
    fallbacks.push(t);
  const params: Record<string, Value> = {};
  for (const [k, v] of Object.entries(l.params ?? {})) if (k !== 'type' && k !== 'data') params[k] = v;
  const theme: Record<string, Value> = {};
  for (const c of ThemeChannelNames) if (param(l, c) !== undefined) theme[c] = param(l, c)!;
  const base = {
    ...(address ? { address } : {}),
    index,
    type,
    fallbacks,
    params,
    theme,
  };
  const layerType = LayerTypes[type];
  if (!layerType) {
    report('seq.view.layer.type', 'The layer’s type is not one this version draws', at);
    return { ...base, placed: false };
  }
  // a section with no system was reported; its layers are not placed
  if (system === undefined) return { ...base, placed: false };
  // a template's layer without data is placed by each instance on its node
  if (template && param(l, 'data') === undefined && type !== 'letters' && type !== 'ruler')
    return { ...base, placed: false };

  const data = dataOf(model, layerData(model.spec, param(l, 'data')), type, system);
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

/**
 * A layer's data with a user-defined ref followed: a ref (a chain of them) naming one address stands for it
 * (_spec/README.md, Language, Refs); a set, or a ref that cannot be followed, is left as written (and names nothing).
 */
function layerData(spec: Node, value: Value | undefined): Value | undefined {
  if (typeof value !== 'string') return value;
  const f = followRef(spec, value);
  return 'addresses' in f && f.addresses.length === 1 ? f.addresses[0] : value;
}

/** A layer's data: the section's letters when it names none; an address's node; a query's first node. */
function dataOf(model: StaticModel, value: Value | undefined, type: string, section: string): LayerData | undefined {
  if (value === undefined)
    return type === 'letters' || type === 'ruler'
      ? { kind: 'address', address: section, system: section, shape: 'letters' }
      : undefined;
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
 * A leader's continuous positions `[lo, hi]` (the edges of its window, lowest first) carried to the follower's
 * positions, `[lo, hi]` lowest first: the low edge's continuous image (`point`) and the high one's (`projectEnd`, the
 * far edge of the place before it), so a follower moves smoothly and is exact at element edges. Where an edge maps
 * nowhere (an intron, a UTR), the hull of the whole positions' image there (`followRange`). Ends past each other: over a
 * cyclic follower's origin the wrap (`lo > hi`), or its hull with `arc: 'hull'`; on a linear follower the hull (the
 * leader's range wraps the coding sequence's origin). Undefined when nothing maps.
 */
export function followSpan(
  section: ResolvedSection,
  lo: number,
  hi: number,
  arc: 'shortest' | 'hull' = 'shortest',
): [number, number] | undefined {
  const follow = section.follow;
  if (!follow) return undefined;
  let a = follow.point(lo);
  let b = projectEnd(follow.path, hi);
  if (a === undefined || b === undefined) {
    const hull = followRange(section, [Math.floor(lo), Math.ceil(hi)], { arc });
    if (!hull || hull.length < 2) return undefined;
    const first = hull[0];
    const last = hull[hull.length - 1];
    if (a === undefined) a = follow.dir === 1 ? first : last;
    if (b === undefined) b = follow.dir === 1 ? last : first;
  }
  const [x, y] = follow.dir === 1 ? [a, b] : [b, a];
  if (x > y || (x === y && hi > lo)) {
    if (section.axis.cyclic && arc !== 'hull') return [x, y];
    const hull = followRange(section, [Math.floor(lo), Math.ceil(hi)], { arc: 'hull' });
    return hull && hull.length >= 2 ? [hull[0], hull[hull.length - 1]] : undefined;
  }
  return [x, y];
}

/** A follower's window in its own offsets: `start` below `end`, or over a cyclic origin `start > end` (`wraps`). */
export interface FollowWindow {
  readonly start: number;
  readonly end: number;
  readonly wraps: boolean;
}

/**
 * A follower's window from its leader's (`start`, `end`: offsets along the leader's axis, display order): the images of
 * the leader's window edges (`followSpan`) as offsets along the follower's axis, in its display order (a reversed
 * follower's lowest position is its window's end). An offset on a leader's break is read as a window's edge
 * (`positionAt`'s sides). `leader` and `axis` are the axes with their pieces (a pending one placed by the caller).
 */
export function followWindow(
  section: ResolvedSection,
  leader: AxisShape,
  axis: AxisShape,
  start: number,
  end: number,
  arc: 'shortest' | 'hull' = 'shortest',
): FollowWindow | undefined {
  const ps = positionAt(leader, start, 'start');
  const pe = positionAt(leader, end, 'end');
  if (ps === undefined || pe === undefined) return undefined;
  const span = followSpan(section, Math.min(ps, pe), Math.max(ps, pe), arc);
  if (!span) return undefined;
  const a = nearestOffsetOf(axis, span[0]);
  const b = nearestOffsetOf(axis, span[1]);
  if (a === undefined || b === undefined) return undefined;
  const [s, e] = axis.dir === -1 ? [b, a] : [a, b];
  if (s > e && axis.cyclic) return { start: s, end: e, wraps: true };
  return { start: Math.min(s, e), end: Math.max(s, e), wraps: false };
}

/**
 * A follower's range from its leader's, both as `[start, end, …]` pieces in display order (an axis's shape). A leader
 * piece that wraps its cyclic origin (`start > end`) is split on it first. The image's hull, or on a cyclic follower
 * the shortest arc that holds it (EC-12), as two pieces when it wraps (in the follower's display order: on a reversed
 * axis the piece below the origin first); `arc: 'hull'` takes the hull there too. Undefined when nothing maps.
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
  return arc.length === 4 && section.axis.dir === -1 ? [arc[2], arc[3], arc[0], arc[1]] : arc;
}
