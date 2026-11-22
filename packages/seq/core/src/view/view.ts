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
  projectPoint,
  type StaticModel,
  viaList,
} from '@molstar/seqx-core/resolve';
import { viaUnusedIssue } from '@molstar/seqx-core/resolve/maps-via';
import { wholePieces } from './axis.js';
import {
  assemblyBlocks,
  type BlockAxes,
  type BlockGroup,
  blockAxes,
  blockNodes,
  openedComplex,
  unknownOverrides,
  sectionAssemblies,
  sectionComplex,
  systemRef,
} from './blocks.js';
import type { ResolvedBridge, StackEntry } from './bridges.js';
import { ViewIssues } from './issues.js';
import { dataOf, layerData } from './layer-data.js';
import { type MemberAxis, memberAxis } from './members.js';

/** The six shapes a layer reads (an extension's type reads one of them: _spec/README.md "Context", Extensions). */
export type LayerShape = 'letters' | 'ranges' | 'values' | 'cells' | 'pairs' | 'grid' | 'matrix';
type Shape = LayerShape;

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
  // regions outlined across tracks: a reader without frames draws none (an outline drawn as bars would read as data)
  frame: { shape: 'ranges' },
  // an alignment's rows and its logo: no fallback (a reader without them draws none)
  rows: { shape: 'cells' },
  logo: { shape: 'cells' },
  // a chromosome's stained bands (an ideogram): a reader without them draws the ranges
  bands: { shape: 'ranges', fallback: 'ranges' },
  // pairs drawn as their ends where a reader lacks arcs; a grid has none (a reader without it draws none)
  pairs: { shape: 'pairs', fallback: 'ranges' },
  grid: { shape: 'grid' },
  // samples × positions: a reader without it draws the states' mean as values
  matrix: { shape: 'matrix', fallback: 'values' },
};

/**
 * The layer kinds a reader knows (_spec/README.md "Resolved view", Layer kinds): the built-ins' (`LayerTypes`) and its
 * composition's extensions' types, each its shape and fallback, and the extensions composed by id with how each reads a version
 * (what the spec's `extensions` is read against). A composition makes one (`SeqContext.kinds`); resolving without one
 * reads the built-ins'.
 */
export interface LayerKinds {
  readonly types: Readonly<Record<string, LayerType>>;
  /** Each extension composed, by id: how it reads a spec's `extensions` version for it (`old`, `newer`). */
  readonly extensions?: Readonly<Record<string, (want: string) => 'old' | 'newer' | undefined>>;
  /**
   * The navigator pictures it draws, by name (its `navigator-picture` entries): a `content` it lacks is the axis, said
   * (none when absent). The axis is always drawn (the fixed fallback).
   */
  readonly pictures?: readonly string[];
}

/**
 * The built-in layer kinds: what a view resolved without a composition reads.
 * @internal
 */
export const BuiltInLayerKinds: LayerKinds = { types: LayerTypes, pictures: ['tracks', 'bands', 'rows'] };

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
  /**
   * A dynamic node's access plan and its system's length (a template's instance: the overlay's, not the model's), or
   * data over an assembly's plan on the section's or block's sequence (the annotation's `bySystem`).
   */
  readonly access?: AccessPlan;
  readonly length?: number;
  /**
   * A query's form after its annotation (`<annotation>/(except :of (item :consequence missense))`): the filter its rows
   * pass, read row by row by the layers (_spec/README.md "Themes", Filters): the user's, a trailing `(except …)` or any
   * lone form but an item form; counted ("n filtered").
   */
  readonly filter?: Node;
  /** A leading item form (`<annotation>/(item :category Domain)/(except …)`): what the track is, its rejected rows not counted. */
  readonly base?: Node;
  /** The forms after the annotation are none the layers read row by row (`seq.render.filter`). */
  readonly unread?: boolean;
}

export interface ResolvedLayer {
  /**
   * Its track's address and its id: every layer of a valid spec has one (_spec/README.md, Language), so undefined only
   * on a tree read without the language's checks; `index` is its place in its track.
   */
  readonly address?: string;
  readonly index: number;
  readonly type: string;
  /**
   * Its type, then the fallbacks a reader without it draws, in order: its `fallback` param after its type (an
   * extension's type without one: none), then the built-ins' chain from there. The type of an extension this reader has
   * at an older major (0.x: minor) than the spec's `extensions` says is left out (its fallbacks drawn).
   */
  readonly fallbacks: readonly string[];
  /** The shape the first of its fallbacks this reader knows reads; absent when it knows none. */
  readonly shape?: LayerShape;
  /**
   * Its params other than `type` and `data`, as written: an extension's type's plain params and the namespaced keys of
   * extensions (`rondo:bundle`) included, carried untouched.
   */
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
  /** Px inside its well above and below its data, as written (_spec/README.md "Render", Tracks layout). */
  readonly padding?: { readonly top: number; readonly bottom: number };
  /** Its layers; a template's layers without data are not placed here: each instance places them on its node. */
  readonly layers: readonly ResolvedLayer[];
  readonly template?: TrackTemplate;
}

export interface ResolvedSection {
  readonly address: string;
  /**
   * Its axis; with blocks (_spec/README.md "Resolved view", Blocks) their display units in order: on one system in a
   * row, the blocks' pieces end to end (the breaks between them); else units `[0, total)` on no system (its
   * `system` the section's address, `pieces: 'blocks'`), each block's own axis in `blocks`.
   */
  readonly axis: Axis;
  readonly follow?: Follow;
  readonly tracks: readonly ResolvedTrack[];
  /** Its layout: its own `layout`, else its view's. */
  readonly layout: string;
  /**
   * Its layout's params as written (text: `group`, `numbers`, `group_labels`, `align`, `max_width`, `height`; rings:
   * `rotation`, `labels`, `height`; blocks: `arrange`, `copies`; tracks and rings: `heights`, the track sizing), and the
   * extensions' namespaced keys (`rondo:bundle`), carried for the extension that reads them.
   */
  readonly params: Readonly<Record<string, Value>>;
  /** Its blocks in order (copies merged by `copies`), when it has `block` children. */
  readonly blocks?: readonly ResolvedBlock[];
  /** A section of several assemblies: a group each, its blocks a run (EC-210; a wider gap and a label in the layouts). */
  readonly groups?: readonly BlockGroup[];
  /** A section of assemblies: the blocks past `maxAssemblyBlocks`, drawn as one closing run (absent when none). */
  readonly more?: number;
  /**
   * The categorical domains of the built-in fields (_spec/README.md "Themes", Built-in fields; EC-218, EC-191): `block`
   * its blocks' labels in order, `system` its blocks' systems' names, over assemblies their names in each assembly's
   * `order` before any `include` (a sequence one hue on every section of its assembly).
   */
  readonly domains?: { readonly block: readonly string[]; readonly system: readonly string[] };
  /** A grid section's second axis (_spec/README.md "Resolved view", Grids): down the plot. */
  readonly y?: GridAxis;
  /** A grid section's x when it is categorical (`of` a rows system, EC-35): its members' order (`member_order`). */
  readonly members?: MemberAxis;
  /**
   * The maps its `via` names (sorted; absent when none): its follow, its `direction: { along }`, its layers' data and
   * what is projected onto it take `model.path(from, to, via)` (_spec/README.md "Resolved view", Via).
   */
  readonly via?: readonly string[];
}

/**
 * A grid's y (`of_y`, the section's system when absent: the whole system, top down) and what it shows (`follow_y`):
 * `x` x's window (its path from x's system when y is another: undefined when the same), a section's region through
 * `follow` (its path from the leader's system), or its own window (`viewport.y`) when neither.
 */
export interface GridAxis {
  readonly axis: Axis;
  readonly followsX: boolean;
  readonly follow?: Follow;
  /**
   * A grid over two sections of blocks (`follow_y` naming a second section of blocks, EC-209): that section's blocks
   * down the plot, `axis` their display axis and `follow` the identity on their units (absent: y is x's blocks).
   */
  readonly blocks?: readonly ResolvedBlock[];
  /** A categorical y (`of_y` a rows system, EC-35): its members' order (`member_order_y`). */
  readonly members?: MemberAxis;
}

/**
 * A block of a section (_spec/README.md "Resolved view", Blocks): its own axis (system, range, direction), where its
 * display units start in the section's, its label, the blocks it stands for (merged copies, `×N`), and each track's
 * layers resolved onto its system (a path a block: a layer whose data reaches it is placed there).
 */
export interface ResolvedBlock {
  readonly address: string;
  readonly label: string;
  readonly axis: Axis;
  /**
   * How it is drawn: `whole` (one chip over its range: a molecule's or a unit of one's default), `residues`, `atoms`
   * (its axis on the atoms system).
   */
  readonly form: BlockForm;
  /** Its first display unit in the section's (0 in a stack; undefined while a length before it is unknown). */
  readonly offset?: number;
  /** The blocks it stands for: itself, and the copies merged into it. */
  readonly members: readonly string[];
  /** Each member's system (a merged copy's unit), in `members`' order: what an end on a copy is drawn on. */
  readonly systems: readonly string[];
  /**
   * @internal Its copy group's key, read by the pairs layer alone (each copy of a pair within one system): blocks of one
   * key are copies (units of one copy group, or the same system, with the same range and form), the key `copies:
   * 'merge'` merges by, so the rule is written once.
   */
  readonly copy: string;
  /** Each of the section's tracks' layers on this block, by the track's index. */
  readonly layers: readonly (readonly ResolvedLayer[])[];
  /** Its display units a position (`block.scale`, absent: 1): its units in the section are its positions times it. */
  readonly scale?: number;
  /** Its group's index in the section's `groups` (a section of several assemblies). */
  readonly group?: number;
  /** Its system's name (the `system` field's value: a derived sequence's printed name). */
  readonly name?: string;
}

/** A block's form (_spec/README.md "Resolved view", Blocks). */
export type BlockForm = 'whole' | 'residues' | 'atoms';

/** A y block's layers: none (its section's tracks are drawn by the x blocks' strips). */
const NoBlockLayers: readonly (readonly ResolvedLayer[])[] = [];

/** The params a section's layout reads (_spec/README.md "Resolved view"). */
const LayoutParams = [
  'group',
  'numbers',
  'group_labels',
  'align',
  'max_width',
  'height',
  'rotation',
  'labels',
  'arrange',
  'copies',
  'aspect',
  'heights',
] as const;

/**
 * A path between two sections of a view, precomposed: hover and selection are projected along it. A section of blocks
 * links through each of its blocks' systems (`fromSystem`, `toSystem`: the systems the path joins).
 */
export interface ViewLink {
  readonly from: string;
  readonly to: string;
  readonly fromSystem: string;
  readonly toSystem: string;
  readonly path: Path;
}

/** The channel slots a layer has (`color`, `ink`, `size`, `opacity`, `line`, `arrow`, `shape`). */
export const ThemeChannelNames = ['color', 'ink', 'size', 'opacity', 'line', 'arrow', 'shape'] as const;

/**
 * A navigator of the view (_spec/README.md "Resolved view", Navigators): its params as written with their defaults,
 * the sections it binds (those the view has, in order; `missing` those it names that the view lacks), the system it
 * draws and the static path from each bound section's system to it.
 */
export interface ResolvedNavigator {
  readonly address: string;
  readonly id: string;
  readonly of: readonly string[];
  readonly missing: readonly string[];
  readonly axes: 'positions' | 'rows' | 'both';
  /** The picture it draws: the one `content` names where the reader has it, else `axis`. */
  readonly content: string;
  /** Its params as written (a namespaced picture's own, plain: what its entry reads). */
  readonly params: Readonly<Record<string, Value>>;
  /** The tracks drawn small (`content: tracks`), by address; the section's but its rulers when absent. */
  readonly tracks?: readonly string[];
  readonly bands?: string;
  /** The system it draws: `system`, else its first bound section's; undefined when it binds none. */
  readonly system?: string;
  /** Whether it draws its first bound section's own axis (no `system`: its range, pieces, direction), else `system` whole. */
  readonly sectionAxis: boolean;
  readonly place: 'top' | 'bottom' | 'left' | 'right';
  readonly size?: number | readonly [number, number];
  readonly hidden: boolean;
  /** The path from each bound section's system to `system`, by section; absent where none joins them. */
  readonly paths: Readonly<Record<string, Path>>;
  /**
   * A navigator over blocks (_spec/README.md "Resolved view", Navigators): its picture a section of blocks, the bound
   * section's own (no `system`: its `of`, `include` and block children) or of the assemblies `system` names (one or a
   * list: `assemblies`). Its windows are carried by `{ system, position }` through the places of each block (identity
   * in units where the block lists are equal), never along `paths` (each such section's path is its own identity).
   */
  readonly blocks?: true;
  /** The assemblies `system` names (one or a list), the picture `of` them; absent: the bound section's blocks. */
  readonly assemblies?: readonly string[];
  /**
   * `navigator.switch: include` (a click on a block runs `section.include` on the bound section); absent: `none`, the
   * default (optional: a host's own navigator record needs none).
   */
  readonly switch?: 'include';
}

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
  /** Its navigators in spec order, hidden ones included. */
  readonly navigators: readonly ResolvedNavigator[];
  /**
   * Its bridges in spec order (_spec/README.md "Resolved view", Bridges), hidden ones and those not placed included;
   * absent when it has none (or the bridge pass is not composed).
   */
  readonly bridges?: readonly ResolvedBridge[];
  /** Its sections and bridges in the order written (what the layout stacks), beside `bridges`. */
  readonly stack?: readonly StackEntry[];
  readonly issues: readonly Issue[];
}

/**
 * The view's bridge pass (`bridges.ts`), set where a composition draws bridges (render's pairs type registers its load,
 * `withViewKind`): a model of a spec with a bridge waits for it. A view resolved without it has no bridges, its sides
 * stacked as a reader without bridges draws them (_spec/README.md "Resolved view", Bridges).
 * @internal
 */
export const BridgePass: { resolve?: typeof import('./bridges.js').resolveBridges } = {};

/**
 * A bridge in a view of a composition that registers none (the pairs layer type listed without `registerBridges`):
 * said, never stacked silently (_spec/README.md "Resolved view", Bridges).
 */
function unregisteredBridge(view: Node, address: string, report: Report): undefined {
  const b = view.children?.find((c) => c.kind === 'bridge');
  if (b)
    report(
      'seq.view.bridge.unregistered',
      'A bridge is not drawn: the composition registers no bridges; its sides are stacked',
      `${address}/${b.id}`,
    );
  return undefined;
}

/** No systems: what an unmapped system links to besides itself. */
const NoSystems: readonly string[] = [];

const caches = new WeakMap<LayerKinds, WeakMap<StaticModel, Map<string, ResolvedView>>>();

/** The resolved views of a model's spec, in spec order; kept per layer kinds, model and view. */
export function resolveViews(model: StaticModel, kinds: LayerKinds = BuiltInLayerKinds): ResolvedView[] {
  const out: ResolvedView[] = [];
  for (const c of model.spec.children ?? [])
    if (c.kind === 'views' && c.id !== undefined)
      for (const v of c.children ?? [])
        if (v.kind === 'view' && v.id !== undefined) out.push(resolveView(model, printAddress([c.id, v.id]), kinds)!);
  return out;
}

/**
 * A view resolved over the static model with the layer kinds a reader knows (the built-ins' by default), or undefined
 * when the spec has no view at the address.
 */
export function resolveView(
  model: StaticModel,
  address: string,
  kinds: LayerKinds = BuiltInLayerKinds,
): ResolvedView | undefined {
  let cache = caches.get(kinds);
  if (!cache) caches.set(kinds, (cache = new WeakMap()));
  let views = cache.get(model);
  if (!views) cache.set(model, (views = new Map()));
  const hit = views.get(address);
  if (hit) return hit;
  const node = nodeAt(model.spec, address);
  if (node?.kind !== 'view') return undefined;
  const v = build(model, node, address, kinds);
  views.set(address, v);
  return v;
}

/**
 * A view node resolved over the model as if it were at `address`, uncached: what a navigator lays its picture out
 * through (a view of its own over the spec's data, never in the spec).
 */
export function resolveViewNode(
  model: StaticModel,
  node: Node,
  address: string,
  kinds: LayerKinds = BuiltInLayerKinds,
): ResolvedView {
  return build(model, node, address, kinds);
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

/**
 * A track's `padding` as its two sides: a number is both, a record one or both (a side absent is 0); none when absent
 * or not of the language's form (its validation says so). Negative or non-finite sides are 0.
 */
function paddingOf(v: Value | undefined): { top: number; bottom: number } | undefined {
  const side = (x: Value | undefined) => (typeof x === 'number' && Number.isFinite(x) && x > 0 ? x : 0);
  if (typeof v === 'number') return { top: side(v), bottom: side(v) };
  if (v === null || typeof v !== 'object' || Array.isArray(v)) return undefined;
  const r = v as Readonly<Record<string, Value>>;
  return { top: side(r.top), bottom: side(r.bottom) };
}

type Report = (code: keyof typeof ViewIssues.codes, message: string, at: string, cause?: Issue) => void;

/**
 * What resolving a view's layers reads of the reader (its layer kinds, the extensions it has too old) and collects:
 * the layers of each extension's types it lacks, by namespace, each address with what it is drawn as.
 */
interface Kinds {
  readonly types: Readonly<Record<string, LayerType>>;
  readonly outdated: ReadonlySet<string>;
  readonly needs: Map<string, Map<string, string | undefined>>;
}

/**
 * The spec's `extensions` against the reader's (_spec/README.md "Resolved view", Layer kinds): one the reader has at an
 * older major (0.x: minor) is outdated, its types drawn as their fallbacks; at an older minor of the same major they are
 * drawn and that is said (`seq.view.extension.version`, info). Each extension reads the version itself (`reads`).
 */
function outdatedOf(wants: Value | undefined, kinds: LayerKinds, at: string, report: Report): Set<string> {
  const out = new Set<string>();
  for (const [id, want] of Object.entries(wants && typeof wants === 'object' ? wants : {})) {
    const r = typeof want === 'string' ? kinds.extensions?.[id]?.(want) : undefined;
    if (r === 'old') out.add(id);
    else if (r) report('seq.view.extension.version', `\`${id}\` ${want} is newer than this reader’s`, at);
  }
  return out;
}

function build(model: StaticModel, view: Node, address: string, layerKinds: LayerKinds): ResolvedView {
  const issues: Issue[] = [];
  const report: Report = (code, message, at, cause) =>
    issues.push(
      issueOf(ViewIssues, code, message, {
        address: at,
        hint: cause ? `${cause.message}: ${cause.hint ?? ''}`.trim() : ViewIssues.codes[code].hint,
      }),
    );
  const wants = model.spec.params?.extensions as Readonly<Record<string, Value>> | undefined;
  const kinds: Kinds = {
    types: layerKinds.types,
    outdated: outdatedOf(wants, layerKinds, address, report),
    needs: new Map(),
  };
  const sectionNodes = (view.children ?? []).filter((c) => c.kind === 'section' && c.id !== undefined);
  const written = param(view, 'layout');
  const viewLayout = typeof written === 'string' ? written : 'tracks';
  const byAddress = new Map<string, Node>();
  const systemOf = new Map<string, string>();
  for (const s of sectionNodes) byAddress.set(`${address}/${s.id}`, s);
  // a section's `via`: the maps its paths take where rivals join the same systems
  const viaOf = (at: string) => viaList(param(byAddress.get(at)!, 'via'));
  // a grid section with no system of its own that follows a section of blocks: a grid over them (`blocksGrid`)
  const overBlocks = new Map<string, string>();
  // and its y over a second section of blocks (`follow_y` naming it): that section's address
  const overBlocksY = new Map<string, string>();
  // each section's block nodes, derived once a build (an assembly's blocks are a walk of its sequences)
  const nodesOf = new Map<string, Node[]>();
  // a section of assemblies' blocks past `maxAssemblyBlocks`: its closing run's count
  const moreOf = new Map<string, number>();
  const blocksAt = (at: string): Node[] => {
    let bs = nodesOf.get(at);
    if (!bs) nodesOf.set(at, (bs = blockNodes(model, byAddress.get(at)!, at)));
    return bs;
  };
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    // a section of blocks with no `of` (or `of` a complex) is on its first block's system
    const assemblies = sectionAssemblies(model, s);
    const derived = assemblies && assemblyBlocks(model, s, assemblies, at);
    if (derived && derived.total > derived.nodes.length)
      report(
        'seq.view.section.blocks',
        `${derived.total - derived.nodes.length} of ${derived.total} blocks are past \`maxAssemblyBlocks\``,
        at,
      );
    if (derived) nodesOf.set(at, derived.nodes);
    if (derived && derived.total > derived.nodes.length) moreOf.set(at, derived.total - derived.nodes.length);
    const first = blocksAt(at)[0];
    for (const id of unknownOverrides(model, s))
      report('seq.view.block.unit', 'The block overrides no unit the source’s complex has held', `${at}/${id}`);
    // an assembly not yet open has no blocks: pending, never a missing system
    if (assemblies && !first && assemblies.some((a) => a.pending)) {
      systemOf.set(at, '');
      continue;
    }
    const own = sectionComplex(model, s) === undefined && !assemblies ? param(s, 'of') : undefined;
    const system = systemRef(model, own ?? (first ? param(first, 'of') : undefined));
    const leader = exactAddress(param(s, 'follow'));
    const lead = leader === undefined ? undefined : byAddress.get(leader);
    const grid = (param(s, 'layout') ?? viewLayout) === 'grid';
    if (grid && first) {
      // a grid over blocks takes its leader's: it has none of its own (_spec/README.md "Resolved view", A grid over blocks)
      report('seq.view.section.system', 'A grid section has no blocks of its own: it follows a section of blocks', at);
      systemOf.set(at, '');
      continue;
    }
    if (grid && param(s, 'of') === undefined && leader !== undefined && lead && blocksAt(leader).length) {
      // the leader's blocks in a row (the plot's axes are their units): not a stack, not a ring; a second section of
      // blocks named by `follow_y` the same (_spec/README.md "Resolved view", A grid over blocks)
      const yLeader = exactAddress(param(s, 'follow_y'));
      const yLead = yLeader === undefined ? undefined : byAddress.get(yLeader);
      const yBlocks = yLead && yLeader !== leader && blocksAt(yLeader!).length > 0;
      const inRow = (n: Node) => param(n, 'arrange') !== 'stack' && (param(n, 'layout') ?? viewLayout) === 'tracks';
      if (!inRow(lead) || (yBlocks && !inRow(yLead))) {
        report('seq.view.section.system', 'A grid follows a section of blocks in a row only', at);
        systemOf.set(at, '');
        continue;
      }
      overBlocks.set(at, leader);
      if (yBlocks) overBlocksY.set(at, yLeader!);
      systemOf.set(at, '');
      continue;
    }
    // a source's complex not held: pending while its source holds nothing, said once it holds other nodes
    const source = first || own === undefined ? undefined : openedComplex(model, s);
    if (source !== undefined) {
      if (model.opened?.has(source)) report('seq.view.section.complex', 'The source holds no such complex', at);
      systemOf.set(at, '');
      continue;
    }
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
    const r = model.path(ls, system, viaOf(at));
    if (r.path) follows.set(at, { leader, path: r.path });
    else report('seq.view.follow.unplaced', 'The section cannot follow its leader', at, r.issue);
  }
  const axes = new Map<string, Axis>();
  const dirs = new Map<string, 1 | -1>();
  const along = new Map<string, Path>();
  const dirOf = (at: string, depth = 0): 1 | -1 => {
    const known = dirs.get(at);
    if (known !== undefined) return known;
    const d = resolveDirection(
      model,
      at,
      byAddress.get(at)!,
      systemOf.get(at)!,
      follows.get(at),
      along,
      report,
      (l) => (depth < byAddress.size ? dirOf(l, depth + 1) : 1),
      viaOf(at),
    );
    dirs.set(at, d);
    return d;
  };
  const blocksOf = new Map<string, BlockAxes>();
  for (const [at, s] of byAddress) {
    if (overBlocks.has(at)) continue;
    const dir = dirOf(at);
    const bs = blocksAt(at);
    if (bs.length) {
      const b = blockAxes(
        model,
        s,
        bs,
        at,
        dir,
        (param(s, 'layout') ?? viewLayout) === 'rings',
        (n, a, system, d) => axisOf(model, n, a, system, d, undefined, undefined, report),
        (message, a) => report('seq.view.section.system', message, a),
      );
      blocksOf.set(at, b);
      axes.set(at, b.axis);
      continue;
    }
    const axis = axisOf(model, s, at, systemOf.get(at)!, dir, follows.get(at)?.path, along.get(at), report);
    const d = param(s, 'direction');
    const ref = d && typeof d === 'object' && !Array.isArray(d) ? exactAddress((d as Record<string, Value>).along) : '';
    axes.set(at, ref ? { ...axis, along: ref } : axis);
  }

  for (const [at, leader] of overBlocks) {
    const b = blocksOf.get(leader);
    if (!b) {
      report('seq.view.section.system', 'The grid follows no section of blocks', at);
      continue;
    }
    blocksOf.set(at, b);
    axes.set(at, b.axis);
    follows.set(at, { leader, path: identityPath(b.axis.system) });
  }

  const sections: ResolvedSection[] = [];
  for (const s of sectionNodes) {
    const at = `${address}/${s.id}`;
    const axis = axes.get(at) ?? { system: at, cyclic: false, dir: 1, pieces: 'blocks' };
    const known = model.systems.has(axis.system);
    const via = viaOf(at);
    let follow: Follow | undefined;
    const f = follows.get(at);
    const leaderAxis = f && axes.get(f.leader);
    if (f && leaderAxis) follow = followOf(f.leader, leaderAxis, f.path);
    const tracks: ResolvedTrack[] = [];
    const ba = blocksOf.get(at);
    const blockLayers: ResolvedLayer[][][] = ba ? ba.blocks.map(() => []) : [];
    for (const t of s.children ?? []) {
      if (t.kind !== 'track' || t.id === undefined) continue;
      // a hidden track is kept in the spec and left out of the view: its row collapses
      if (param(t, 'hidden') === true) continue;
      const ta = `${at}/${t.id}`;
      const layers: ResolvedLayer[] = [];
      const template = templateOf(param(t, 'for_each'));
      const perBlock: ResolvedLayer[][] = blockLayers.map(() => []);
      (t.children ?? []).forEach((l, index) => {
        if (l.kind !== 'layer') return;
        if (!ba) {
          const system = known ? axis.system : undefined;
          layers.push(layerOf(model, kinds, l, index, ta, system, report, template !== undefined, via));
          return;
        }
        // a path a block: placed on each block its data reaches; said only when it reaches none
        const said: Parameters<Report>[] = [];
        let placed: ResolvedLayer | undefined;
        for (let k = 0; k < ba.blocks.length; k++) {
          const system = ba.blocks[k].axis.system;
          const r = layerOf(model, kinds, l, index, ta, system, (...a) => said.push(a), !!template, via);
          perBlock[k].push(r);
          placed ??= r.placed ? r : undefined;
        }
        if (!placed && said.length) report(...said[0]);
        const any = placed ?? perBlock[0]?.[perBlock[0].length - 1];
        if (any) layers.push(any);
      });
      for (let k = 0; k < perBlock.length; k++) blockLayers[k].push(perBlock[k]);
      const label = param(t, 'label');
      const height = param(t, 'height');
      const padding = paddingOf(param(t, 'padding'));
      tracks.push({
        address: ta,
        label: typeof label === 'string' ? label : t.id,
        ...(typeof height === 'number' ? { height } : {}),
        ...(padding ? { padding } : {}),
        layers,
        ...(template ? { template } : {}),
      });
    }
    const own = param(s, 'layout');
    const params: Record<string, Value> = {};
    for (const k of LayoutParams) {
      const v = param(s, k);
      if (v !== undefined) params[k] = v;
    }
    // an extension's keys on the section (`rondo:bundle`): carried for the extension that reads them
    for (const [k, v] of Object.entries(s.params ?? {})) if (k.includes(':')) params[k] = v;
    const blocks = ba?.blocks.map((b, k) => ({ ...b, layers: blockLayers[k] }));
    const more = moreOf.get(at);
    const assemblies = blocks && sectionAssemblies(model, s);
    const domains = blocks && {
      block: blocks.map((b) => b.label),
      // over assemblies their names in each one's `order` before any `include` (EC-191)
      system: assemblies
        ? assemblies.flatMap((a) => a.order.map((i) => a.sequences[i].name))
        : [...new Set(blocks.map((b) => b.name ?? b.axis.system))],
    };
    const grid = (typeof own === 'string' ? own : viewLayout) === 'grid';
    const members = grid ? memberAxis(axis.system, param(s, 'member_order')) : undefined;
    const yBlocks = overBlocksY.get(at);
    const yb = yBlocks === undefined ? undefined : blocksOf.get(yBlocks);
    const y =
      (typeof own === 'string' ? own : viewLayout) === 'grid'
        ? yb
          ? {
              axis: yb.axis,
              followsX: false,
              follow: followOf(yBlocks!, yb.axis, identityPath(yb.axis.system)),
              blocks: yb.blocks.map((b) => ({ ...b, layers: NoBlockLayers })),
            }
          : overBlocks.has(at)
            ? { axis, followsX: true }
            : gridAxis(model, s, at, axis, axes, systemOf, report, via)
        : undefined;
    sections.push({
      address: at,
      axis,
      ...(follow ? { follow } : {}),
      tracks,
      layout: typeof own === 'string' ? own : viewLayout,
      params,
      ...(blocks ? { blocks } : {}),
      ...(ba?.groups ? { groups: ba.groups } : {}),
      ...(more ? { more } : {}),
      ...(domains ? { domains } : {}),
      ...(y ? { y } : {}),
      ...(members ? { members } : {}),
      ...(via.length ? { via } : {}),
    });
  }

  // every path between the view's sections, composed now so hover never composes; a section of blocks shows its
  // blocks' systems (its own axis is their units, on no system), each one linked
  const links: ViewLink[] = [];
  const shown = (s: ResolvedSection) => [...new Set(s.blocks ? s.blocks.map((b) => b.axis.system) : [s.axis.system])];
  // linear in the blocks: a system with no map joins only itself (a genome's thousands of contigs on two sections
  // are their identities, never every pair asked for a path)
  const shownOf = new Map(sections.map((s) => [s, shown(s)]));
  const shownSet = new Map(sections.map((s) => [s, new Set(shownOf.get(s))]));
  for (const a of sections)
    for (const b of sections) {
      if (a === b) continue;
      for (const from of shownOf.get(a)!)
        for (const to of model.systems.has(from) && model.edgesFrom(from).length > 0
          ? shownOf.get(b)!
          : shownSet.get(b)!.has(from)
            ? [from]
            : NoSystems) {
          if (!model.systems.has(from) || !model.systems.has(to)) continue;
          // what is projected onto a section or from it takes its `via` (both ends' named maps)
          const via = a.via && b.via ? [...new Set([...a.via, ...b.via])] : (b.via ?? a.via);
          let path = from === to ? identityPath(from) : model.path(from, to, via).path;
          // both ends naming rival maps of one pair leave it tied: the target's map carries it, and that is said
          if (!path && a.via && b.via) {
            path = model.path(from, to, b.via).path;
            if (path)
              report(
                'seq.view.via.conflict',
                `${a.address} and ${b.address} name rival maps in their via; ${b.address}'s carries what reaches it`,
                b.address,
              );
          }
          if (path) links.push({ from: a.address, to: b.address, fromSystem: from, toSystem: to, path });
        }
    }
  for (const s of sections) if (s.via) issues.push(...unusedVia(s, along.get(s.address), links));
  // each extension the layers need and this reader lacks (or has too old), said once a namespace
  for (const [id, layers] of kinds.needs) {
    const as = [...new Set(layers.values())].map((t) => (t ? `drawn as \`${t}\`` : 'not drawn'));
    report(
      'seq.view.layer.extension',
      `${layers.size} layer(s) need \`${id}\`${wants?.[id] ? ` ${wants[id]}` : ''}: ${as}`,
      address,
    );
  }
  const scales: Record<string, Readonly<Record<string, Value>>> = {};
  for (const c of model.spec.children ?? [])
    if (c.kind === 'views' && c.id !== undefined)
      for (const n of c.children ?? [])
        if (n.kind === 'scale' && n.id !== undefined) scales[`/${c.id}/${n.id}`] = n.params ?? {};
  const ex = param(view, 'export');
  return {
    address,
    layout: viewLayout,
    scales,
    ...(ex && typeof ex === 'object' && !Array.isArray(ex) ? { export: ex as Record<string, Value> } : {}),
    sections,
    links,
    navigators: navigatorsOf(model, view, address, sections, report, layerKinds),
    // its bridges and the stack, where the view has a bridge and the pass is composed (_spec/README.md "Resolved
    // view", Bridges)
    ...(BridgePass.resolve
      ? BridgePass.resolve({ view, address, sections, model, kinds, layer: layerOf }, issues)
      : unregisteredBridge(view, address, report)),
    issues,
  };
}

/** A view's navigators (`<view>/navigators/<id>`), resolved against its sections. */
function navigatorsOf(
  model: StaticModel,
  view: Node,
  address: string,
  sections: readonly ResolvedSection[],
  report: Report,
  kinds: LayerKinds,
): ResolvedNavigator[] {
  const out: ResolvedNavigator[] = [];
  for (const c of view.children ?? []) {
    if (c.kind !== 'navigators') continue;
    for (const n of c.children ?? [])
      if (n.kind === 'navigator' && n.id !== undefined) {
        const nav = resolveNavigatorNode(model, sections, n, `${address}/${c.id}/${n.id}`, kinds);
        out.push(nav);
        // a picture this reader lacks (an extension's) is drawn as the axis, and said
        const content = param(n, 'content');
        if (content !== undefined && content !== nav.content)
          report('seq.view.navigator.picture', `No picture \`${content}\` here: drawn as the axis`, nav.address);
      }
  }
  return out;
}

/**
 * A `navigator` node resolved against a view's sections, as if at `at`: a node of the spec, or a host's definition of
 * the same shape (never saved); its picture one of the reader's (`kinds.pictures`), else the axis.
 */
export function resolveNavigatorNode(
  model: StaticModel,
  sections: readonly ResolvedSection[],
  n: Node,
  at: string,
  kinds: LayerKinds = BuiltInLayerKinds,
): ResolvedNavigator {
  const written = param(n, 'of');
  const named = (Array.isArray(written) ? written : [written])
    .map((v) => exactAddress(v))
    .filter((v): v is string => v !== undefined);
  const of = named.filter((a) => sections.some((s) => s.address === a));
  const first = sections.find((s) => s.address === of[0]);
  // a navigator over blocks: `system` an assembly or a list of them, or no `system` over a section of blocks
  const sys = param(n, 'system');
  const listed = (Array.isArray(sys) ? sys : [sys]).map((v) => exactAddress(v)).filter((v) => v !== undefined);
  const assemblies = listed.length && listed.every((a) => model.assemblies.has(a)) ? listed : undefined;
  const blocks = !!assemblies || (sys === undefined && !!first?.blocks);
  const system = assemblies ? assemblies[0] : (exactAddress(sys) ?? first?.axis.system);
  const paths: Record<string, Path> = {};
  for (const a of of) {
    const s = sections.find((x) => x.address === a)!;
    // over blocks each bound section's window is carried by its places: its own identity marks it joined
    if (blocks) {
      if (!assemblies || overAssemblies(model, s, assemblies)) paths[a] = identityPath(s.axis.system);
      continue;
    }
    if (system === undefined || !model.systems.has(system) || !model.systems.has(s.axis.system)) continue;
    const path = s.axis.system === system ? identityPath(system) : model.path(s.axis.system, system, s.via).path;
    if (path) paths[a] = path;
  }
  const tracks = param(n, 'tracks');
  const axes = (param(n, 'axes') as ResolvedNavigator['axes'] | undefined) ?? 'positions';
  const size = param(n, 'size');
  const content = param(n, 'content');
  return {
    address: at,
    id: n.id ?? at.slice(at.lastIndexOf('/') + 1),
    of,
    missing: named.filter((a) => !of.includes(a)),
    axes,
    // a picture this reader lacks (an extension's `<id>:<name>`, one its composition does not list) draws the axis
    content: kinds.pictures?.includes(content as string) ? (content as string) : 'axis',
    params: n.params ?? {},
    ...(Array.isArray(tracks)
      ? { tracks: tracks.map((t) => exactAddress(t)).filter((t): t is string => t !== undefined) }
      : {}),
    ...(exactAddress(param(n, 'bands')) !== undefined ? { bands: exactAddress(param(n, 'bands'))! } : {}),
    ...(system !== undefined ? { system } : {}),
    sectionAxis: param(n, 'system') === undefined,
    place: (param(n, 'place') as ResolvedNavigator['place'] | undefined) ?? (axes === 'rows' ? 'right' : 'top'),
    ...(typeof size === 'number' || Array.isArray(size) ? { size: size as ResolvedNavigator['size'] } : {}),
    hidden: param(n, 'hidden') === true,
    paths,
    ...(blocks ? { blocks: true as const } : {}),
    ...(assemblies ? { assemblies } : {}),
    ...(param(n, 'switch') === 'include' ? { switch: 'include' as const } : {}),
  };
}

/** Whether a section shows sequences of the assemblies (its blocks', else its axis's system): a navigator over them reaches it. */
function overAssemblies(model: StaticModel, s: ResolvedSection, assemblies: readonly string[]): boolean {
  const of = (system: string) => {
    const at = system.lastIndexOf('/(sequence ');
    return at > 0 && assemblies.includes(system.slice(0, at));
  };
  if (s.blocks) return s.blocks.length > 0 && s.blocks.every((b) => of(b.axis.system));
  // a section of an assembly still pending has no blocks yet: reached once they come
  const o = s.params.of;
  if (typeof o === 'string' && assemblies.includes(o)) return true;
  if (Array.isArray(o) && o.length && o.every((x) => typeof x === 'string' && assemblies.includes(x))) return true;
  return of(s.axis.system);
}

/**
 * A section's `via` maps that join no pair on any of its paths (its follow, its along, its layers' and blocks' data,
 * its grid's y, the links onto it and from it: what hover, selection and the focus hull are projected along): each a
 * `seq.resolve.via.unused` warning at the section naming the map (the author, 2026-11-20: "Warn"); its paths were
 * chosen as if the via were absent (maps-via.ts `viaPath`).
 */
function unusedVia(s: ResolvedSection, along: Path | undefined, links: readonly ViewLink[]): Issue[] {
  const used = new Set<string>();
  const take = (p: Path | undefined) => {
    for (const e of p?.edges ?? []) used.add(e.map);
  };
  take(s.follow?.path);
  take(along);
  take(s.y?.follow?.path);
  for (const t of s.tracks) for (const l of t.layers) take(l.path);
  for (const b of s.blocks ?? []) for (const t of b.layers) for (const l of t) take(l.path);
  for (const l of links) if (l.to === s.address || l.from === s.address) take(l.path);
  return (s.via ?? []).filter((m) => !used.has(m)).map((m) => viaUnusedIssue(s.address, m));
}

/** A grid section's y (`GridAxis`): its system whole, and what it follows. */
function gridAxis(
  model: StaticModel,
  s: Node,
  at: string,
  x: Axis,
  axes: ReadonlyMap<string, Axis>,
  systemOf: ReadonlyMap<string, string>,
  report: Report,
  via: readonly string[],
): GridAxis | undefined {
  const system = systemRef(model, param(s, 'of_y')) ?? x.system;
  const entry = model.systems.get(system);
  if (!entry) {
    report('seq.view.section.system', 'The grid’s y shows no system of the spec', at);
    return undefined;
  }
  const length = entry.length;
  const axis: Axis = {
    system,
    ...(length !== undefined ? { length, ranges: [0, length] } : {}),
    cyclic: entry.cyclic,
    dir: 1,
    pieces: 'whole',
  };
  // a categorical y (a rows system, EC-35): its members' order; it follows nothing of positions, so with no `follow_y`
  // over another system than x's it keeps a window of its own
  const members = memberAxis(system, param(s, 'member_order_y'));
  const own = members ? { axis, followsX: false, members } : { axis, followsX: false };
  if (members && param(s, 'follow_y') === undefined && system !== x.system) return own;
  const written = param(s, 'follow_y') ?? 'x';
  if (written === 'none') return own;
  const leader = written === 'x' ? at : exactAddress(written);
  const from = written === 'x' ? x.system : leader === undefined ? undefined : systemOf.get(leader);
  const leaderAxis = written === 'x' ? x : leader === undefined ? undefined : axes.get(leader);
  if (leader === undefined || from === undefined || !leaderAxis || !model.systems.has(from)) {
    report('seq.view.follow.unplaced', 'The grid’s y follows no section of the view', at);
    return own;
  }
  if (written === 'x' && from === system) return { ...own, followsX: true };
  const r = model.path(from, system, via);
  if (!r.path) {
    report('seq.view.follow.unplaced', 'The grid’s y cannot follow: no path joins the systems', at, r.issue);
    return own;
  }
  return { axis, followsX: written === 'x', follow: followOf(leader, leaderAxis, r.path) };
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
  via: readonly string[],
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
  const path = ref ? pathAlong(model, ref, system, via) : undefined;
  if (!path) {
    report('seq.view.direction.unresolved', 'No path reaches the section from what it reads along', at);
    return 1;
  }
  along.set(at, path);
  return path.edges.length || path.step !== 0 ? of(path) : 1;
}

/** A map's own path, from its `from` to its `to`. */
function mapPath(map: MapEntry): Path {
  return { from: map.from, to: map.to, step: map.step, linear: map.linear, edges: [] };
}

/** The path from what a section reads along (a system, or a map: its own linear part) to the section's system. */
function pathAlong(model: StaticModel, ref: string, system: string, via: readonly string[]): Path | undefined {
  const map = model.maps.get(ref);
  if (map) {
    const forward = mapPath(map);
    if (map.to === system) return forward;
    if (map.from === system) return inversePath(forward);
    return model.path(map.from, system, via).path;
  }
  if (!model.systems.has(ref) || !model.systems.has(system)) return undefined;
  return ref === system ? identityPath(system) : model.path(ref, system, via).path;
}

/**
 * The axis of a section or a block (`node`, at `at`) on `system`: its pieces (`pieces`, by the path it follows or reads
 * along) cut to its `range`, its direction, cyclic and length from the system; issues reported at `at`.
 */
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
  kinds: Kinds,
  l: Node,
  index: number,
  track: string,
  system: string | undefined,
  report: Report,
  template = false,
  via: readonly string[] = [],
): ResolvedLayer {
  const type = String(param(l, 'type'));
  const address = l.id === undefined ? undefined : `${track}/${l.id}`;
  const at = address ?? `${track}`;
  // its chain (_spec/README.md "Resolved view", Layer kinds): its type, its `fallback` (an extension's type without
  // one: none), then the built-ins' chain from there; the first this reader knows is what it is placed as
  const types = kinds.types;
  const ns = type.includes(':') ? type.slice(0, type.indexOf(':')) : undefined;
  const written = param(l, 'fallback');
  const chain: string[] = [];
  for (let t: string | undefined = type; t !== undefined && t !== 'none' && !chain.includes(t); ) {
    chain.push(t);
    t = t === type && (typeof written === 'string' || ns) ? (written as string | undefined) : types[t]?.fallback;
  }
  const fallbacks = ns && kinds.outdated.has(ns) ? chain.slice(1) : chain;
  const drawn = fallbacks.find((t) => types[t]);
  const layerType = drawn === undefined ? undefined : types[drawn];
  const params: Record<string, Value> = {};
  for (const [k, v] of Object.entries(l.params ?? {})) if (k !== 'type' && k !== 'data') params[k] = v;
  const theme: Record<string, Value> = {};
  for (const c of ThemeChannelNames) if (param(l, c) !== undefined) theme[c] = param(l, c)!;
  const base = {
    ...(address ? { address } : {}),
    index,
    type,
    fallbacks,
    ...(layerType ? { shape: layerType.shape } : {}),
    params,
    theme,
  };
  if (drawn !== type) {
    // an extension's type is said once a namespace (`build`), a bare one (a newer built-in's) each layer
    if (ns) {
      let layers = kinds.needs.get(ns);
      if (!layers) kinds.needs.set(ns, (layers = new Map()));
      layers.set(at, drawn);
    } else report('seq.view.layer.type', 'The layer’s type is not one this version draws', at);
  }
  if (!layerType) return { ...base, placed: false };
  // a section with no system was reported; its layers are not placed
  if (system === undefined) return { ...base, placed: false };
  // a template's layer without data is placed by each instance on its node
  if (template && param(l, 'data') === undefined && drawn !== 'letters' && drawn !== 'ruler')
    return { ...base, placed: false };

  const value = layerData(model.spec, param(l, 'data'));
  // a segments map read as pairs (_spec/README.md "Resolved view", Bridges): a run a pair, `a` its `from` range
  const data = dataOf(model, value, drawn!, system);
  if (!data) {
    report('seq.view.layer.data', 'The layer’s data names nothing the static model holds', at);
    return { ...base, placed: false };
  }
  // a ranges layer over links: their ends on the section's (or the block's) system, a mark each (Render, Layer types)
  if (data.shape === 'pairs' && layerType.shape === 'ranges') {
    // its plan built here, on first ask (resolve's LinkEnds)
    if (model.annotations.get(data.address)?.ends?.plan(system))
      return { ...base, data: { ...data, system }, placed: true };
    report('seq.view.layer.unplaced', 'The links have no end on the section’s system', at);
    return { ...base, data, placed: false };
  }
  if (data.shape !== undefined && data.shape !== layerType.shape) {
    report('seq.view.layer.shape', 'The layer’s type does not draw the shape of its data', at);
    return { ...base, data, placed: false };
  }
  if (data.system === system) return { ...base, data, placed: true };
  const r = model.path(data.system, system, via);
  if (!r.path) {
    report('seq.view.layer.unplaced', 'The layer’s data does not reach its section', at, r.issue);
    return { ...base, data, placed: false };
  }
  return { ...base, data, path: r.path, placed: true };
}

export { followRange, followSpan, followWindow, type FollowWindow } from './follow.js';
export { type DataForms, splitDataForms } from './layer-data.js';
