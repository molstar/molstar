/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The static model (_spec/README.md "Resolve"): tier 1, pure, over a checked spec. Every coordinate system the spec
// declares, its maps resolved to their encodings, paths composed on demand and kept, each source bound to its runtime
// and each annotation's access plan. It implements the query forms' model, so queries run on a real spec. Entries are
// keyed by address and cached by the hashes of the nodes they read, and paths by the map entries their search read: a
// resolver given an edited spec rebuilds only what the edit touched.

import type { ParamsFor } from '@molstar/core/repr/builder';
import {
  type Atom,
  exactAddress,
  hashNode,
  type Issue,
  issueOf,
  type Node,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import {
  AlignmentSystems,
  type Coverage,
  type FormatRuntime,
  type Shape,
  type SourceEntry,
  type SystemInfo,
} from '@molstar/seqx-core/data';
import { type SeqKinds, SeqLanguage } from '@molstar/seqx-core/lang';
import {
  type ItemsRead,
  type Location,
  location,
  type Numbering,
  overlapsRanges,
  partOf,
  type QueryAnnotation,
  type QueryItem,
  type QueryModel,
  type QuerySystem,
  rangesLocation,
  unionLocations,
} from '@molstar/seqx-core/query';
import type { GeneticCode } from './codon.js';
import { type GridData, type PairsData } from './held.js';
import {
  buildLinks,
  buildMolecular,
  type Cached,
  compoundEntry,
  type Placed,
  placeChildren,
} from './model-molecules.js';
import { ResolveIssues } from './issues.js';
import type { Encoding, Linear } from './linear.js';
import {
  atomsOf,
  type CompoundEntry,
  complexesUnder,
  componentsOf,
  EmptyMolecular,
  type Molecular,
  modifiedOf,
  moleculeOf,
  unitsOf,
} from './molecule-query.js';
import { type DescribedMap, describedAddress } from './maps-described.js';
import {
  addressIds,
  alignmentEntry,
  annotationEntry,
  bindSource,
  heldEntry,
  rowsSystemOf,
  sequenceEntry,
} from './model-entries.js';
export { resolveAnnotationNode, resolveMapNode, resolveSequenceNode, rowsSystemOf } from './model-entries.js';
import { mapEntry, numberingEntry } from './maps-entry.js';
import type { StructureMapData } from './maps-segments.js';
import { pathKey, viaList, viaPath, viaRivals, withoutRivals } from './maps-via.js';
import { collapseAtoms, expandAtoms, type ModifiedIndex, type UnitEntry, unitEdges } from './molecules.js';
import { composeEdges, type Edge, findPath, inversePath, type Path, projectRanges, type Step } from './paths.js';

/** The length a system of unknown length (letters from a source, a streamed genome) answers the query forms with. */
export const OpenLength = 0x7fffffff;

/** The choices of the resolve pass, each with its recommended default (`DefaultResolveConfig`). */
export interface ResolveConfig {
  /** Counted maps a path takes at most: 4. */
  readonly maxPathSteps: number;
  /** Shortest candidates compared at most before a path is refused as over the limit: 16. */
  readonly maxPathCandidates: number;
  /** How bounded maps are held: the cheapest encoding that answers (default), or all as runs, or all expanded. */
  readonly encoding: 'cheapest' | 'runs' | 'expanded';
  /**
   * The length the query forms read on a system whose length is unknown, so labels resolve there before its source is
   * opened: int32's max. Layout never reads it: an axis of unknown length is pending.
   */
  readonly openQueryLength: number;
  /** Levels of sources nested in dynamic nodes the overlay requests at most (EC-179): 4. */
  readonly maxNestedDepth: number;
  /** Units a complex holds at most (a large capsid's assembly): 65,536; a complex past it is refused. */
  readonly maxUnits: number;
}

export const DefaultResolveConfig: ResolveConfig = {
  maxPathSteps: 4,
  maxPathCandidates: 16,
  encoding: 'cheapest',
  openQueryLength: OpenLength,
  maxNestedDepth: 4,
  maxUnits: 65536,
};

/**
 * A coordinate system the spec declares: a sequence (an alignment's columns and rows too), a molecule's residues or
 * atoms, a unit of a complex; its length when known.
 */
export interface SystemEntry {
  readonly address: string;
  readonly kind: 'sequence' | 'molecule' | 'atoms' | 'unit';
  /** A sequence's, or a unit's of one; none for a molecule's residues and atoms. */
  readonly alphabet?: 'protein' | 'dna' | 'rna';
  /** Undefined: unknown until its source is opened (pending, never a guess). */
  readonly length?: number;
  readonly cyclic: boolean;
  readonly letters?: string;
  /** The source of its letters. */
  readonly source?: string;
  /** The name its source knows it by: its `name`, else its id (an alignment's columns: `<alignment>/columns`). */
  readonly name: string;
  /** Other names its own source knows it by (an alignment's columns: `columns`, as one alignment's file says it). */
  readonly aliases?: readonly string[];
  /** Its numberings' ids: the built-ins (`index`, `ordinal`), then its `numbering` children's, in order. */
  readonly numberings: readonly string[];
  /** How its letters are read. */
  readonly access: AccessPlan;
  /** A unit's molecule and where the unit starts on it: its position i is the molecule's `start + i`. */
  readonly parent?: { readonly system: string; readonly start: number };
}

/**
 * A declared numbering (a `numbering` child of a sequence, _spec/README.md "Resolve", Numberings): inline labels as a
 * numbering at once; a source's labels are read from its rows by the store's reader (tier 2), so tier 1 holds where.
 */
export interface NumberingEntry {
  /** The numbering node's address (`/data/h/kabat`). */
  readonly address: string;
  /** Its sequence's address. */
  readonly system: string;
  /** Its id: the name forms read (`:numbering kabat`). */
  readonly id: string;
  /** The scheme its labels follow (an open word), when it says. */
  readonly scheme?: string;
  /** Inline labels, as a numbering of the system. */
  readonly numbering?: Numbering;
  /** The source whose rows hold the labels, read on the system by its name (`access`). */
  readonly source?: string;
  /** The field of the source's rows that holds the label (`label` when the node names none). */
  readonly field?: string;
  readonly access?: AccessPlan;
}

/** A map, resolved: its kind, its encoding, and its linear part (fine places of the residue side for a codon map). */
export interface MapEntry {
  readonly address: string;
  readonly kind: 'offset' | 'codon' | 'segments' | 'structure';
  readonly from: string;
  readonly to: string;
  readonly encoding: Encoding;
  readonly step: Step;
  readonly linear: Linear;
  /** A codon map's genetic code: its `table`, 1 when absent; undefined when that code is not built. */
  readonly code?: GeneticCode;
  /** Its rank among maps joining the same systems (`priority`, 0 when absent; _spec/README.md "Resolve", Paths). */
  readonly priority: number;
  /** A structure map's provenance and residue status (`structureStatus`). */
  readonly structure?: StructureMapData;
}

/** The language's params of a source: one definition, the description's. */
export type SourceParams = ParamsFor<SeqKinds, 'source'>;

/** A source as the protocol reads it, with the rest of its params. */
export type BoundSourceEntry = SourceEntry & Omit<SourceParams, 'format' | 'uri' | 'value'>;

/** A source and the runtime bound to it by its format; unavailable without one. */
export interface SourceBinding {
  readonly address: string;
  readonly entry: BoundSourceEntry;
  readonly runtime?: FormatRuntime;
}

/**
 * How a node's data is read: inline data as a `seq` source of one pair over the whole system (a ranges list also at a
 * level whose one tile is the whole system: `inlineSeq`), or a source of the spec by region. `system` is the name the
 * source knows the system by.
 */
export type AccessPlan =
  | {
      readonly kind: 'inline';
      readonly source: SourceEntry;
      readonly runtime: FormatRuntime;
      readonly system: string;
      readonly shape: Shape;
    }
  | { readonly kind: 'source'; readonly source: string; readonly system: string; readonly shape: Shape };

/**
 * An annotation: its shape, the system it is on, its fields and how it is read (none for what is held whole: `pairs`
 * and `grid`, _spec/README.md "Resolve", Pairs and grids).
 */
export interface AnnotationEntry {
  readonly address: string;
  /**
   * `cells`: an alignment, on its columns system (`target`), its rows a system of their own (`rows`); `pairs` and
   * `grid`: two systems, `target` and `ySystem` (a pair's second end, a grid's y).
   */
  readonly shape: 'ranges' | 'values' | 'cells' | 'pairs' | 'grid';
  readonly target: string;
  readonly fields: readonly string[];
  readonly access?: AccessPlan;
  /** The second system of `pairs` (`b`) and `grid` (`y`): `target`'s when the node names none. */
  readonly ySystem?: string;
  /** A `pairs` node's ends, held whole. */
  readonly pairs?: PairsData;
  /** A `grid` node's values, held whole (absent when they were dropped: a wrong length, too many). */
  readonly grid?: GridData;
  /** Inline ranges as items (a wrapped row is two rows of one item, EC-9). */
  readonly items?: readonly QueryItem[];
  /** An alignment's rows system (`rowsSystemOf`): a position a row, in data order (EC-35). */
  readonly rows?: string;
  /** An alignment's reference row, by label (what `identity` compares with). */
  readonly reference?: string;
  /** A complex's links: the units their rows' ends are on (`PairsData.aUnit`, `bUnit` index it), by address. */
  readonly units?: readonly string[];
  /** A complex's links' ends on each system as inline ranges: what a `ranges` layer over links reads (`LinkEnds`). */
  readonly ends?: LinkEnds;
}

/**
 * Links' ends on each system as inline ranges plans (_spec/README.md "Render", Layer types: ranges), each built on first
 * ask: the view's `layerOf` asks for a ranges layer over links on its system. Why lazy: most complexes' links are drawn
 * as arcs only, and the plan encodes a table. The store and the export open the plans built (`built`).
 */
export interface LinkEnds {
  /** The systems with an end. */
  readonly systems: ReadonlySet<string>;
  /** The plan of the ends on `system`, built once; undefined where none is. */
  plan(system: string): AccessPlan | undefined;
  /** The plans asked for so far. */
  built(): Iterable<AccessPlan>;
}

/** A path between two systems, or why there is none. */
export interface PathResult {
  readonly path?: Path;
  readonly issue?: Issue;
  /** With `via`: the named maps the path does not step through (`seq.resolve.via.unused`, maps-via.ts `viaPath`). */
  readonly unused?: readonly string[];
}

export interface Projection {
  readonly location: Location;
  readonly issues: readonly Issue[];
}

export interface StaticModel extends QueryModel {
  readonly spec: Node;
  /** The spec's hash, with the lengths it was resolved with: equal specs and lengths give the same model. */
  readonly hash: string;
  readonly config: ResolveConfig;
  readonly systems: ReadonlyMap<string, SystemEntry>;
  readonly maps: ReadonlyMap<string, MapEntry>;
  readonly sources: ReadonlyMap<string, SourceBinding>;
  readonly annotations: ReadonlyMap<string, AnnotationEntry>;
  /** The declared numberings, by the numbering node's address. */
  readonly numberings: ReadonlyMap<string, NumberingEntry>;
  /** Whether a declared numbering's labels are a source's (read by tier 2): what wraps a query model checks first. */
  readonly hasSourceNumberings: boolean;
  /** The molecules, complexes, units, atoms systems' group maps and modified indices (_spec/README.md "Resolve"). */
  readonly molecular: Molecular;
  /** The compounds described (`compound` nodes), by code. */
  readonly compounds: ReadonlyMap<string, CompoundEntry>;
  /** The maps with an end a source describes (a unit of its complex): the overlay's steps once held (Overlay). */
  readonly describedMaps?: ReadonlyMap<string, DescribedMap>;
  /** A sequence's modified index, or a unit's of one; undefined when it has none. */
  modified(system: string): ModifiedIndex | undefined;
  /** A component described: its `compound` node, else the host's names (`ResolveOptions.componentNames`). */
  compound(code: string): Omit<CompoundEntry, 'id'> | undefined;
  readonly issues: readonly Issue[];
  /** A system's declared numbering by its id (`kabat`), or undefined: never a built-in (`index`, `ordinal`). */
  numberingOf(system: string, id: string): NumberingEntry | undefined;
  /**
   * The shortest unambiguous path, composed once and kept (across resolves while the maps its search read hold); one
   * entry a pair and `via` set: with `via` (map addresses) the named maps' rivals are left out (maps-via.ts).
   */
  path(from: string, to: string, via?: readonly string[]): PathResult;
  /** The steps from a system along the spec's maps, both ways: what a path search reads (the overlay adds its own). */
  edgesFrom(system: string): readonly Edge[];
  /**
   * A Location carried onto a system: each part along its path (rects are not carried yet); a part with no path is
   * left out with the path's issue. With `close`, the location is first closed (`close`). The atoms group map applies
   * at the ends only: a part on atoms is collapsed onto their residues first, and onto an atoms system the image on
   * its residues is expanded last. With `via`, each part's path is `path(from, to, via)`'s.
   */
  project(loc: Location, to: string, options?: { close?: boolean; via?: readonly string[] }): Projection;
  /**
   * Closure, one step (EC-30): nucleotides gain their residues and lift to their whole codons, through each codon map
   * between the systems named (every system when none are), along the map itself.
   */
  close(loc: Location, among?: readonly string[]): Location;
}

export interface ResolveOptions {
  /** The format runtimes the application composed; inline data always reads through the `seq` runtime. */
  readonly runtimes?: readonly FormatRuntime[];
  readonly config?: Partial<ResolveConfig>;
  /**
   * The lengths of source-backed systems that their sources' coverages gave once opened (`sourceLengths`, or
   * `openModel` in render): such a system then resolves as one with inline letters of that length (_spec/README.md
   * "Resolve", Systems). Each is kept by its system and its source node, so a source edited since is unknown again.
   */
  readonly lengths?: Lengths;
  /**
   * The host's identifier resolver (its `identifierResolver` entry): the identifiers it holds the same compound as one
   * (`ccd:ATP`), never the one itself; none by default (exact match only).
   */
  readonly equivalents?: (identifier: string) => readonly string[];
  /**
   * The host's component names (its `componentNames` entry: a CCD lookup) by code, where no `compound` node describes
   * one; none by default (the code alone).
   */
  readonly componentNames?: (code: string) => Omit<CompoundEntry, 'id'> | undefined;
}

/**
 * Source-backed systems' lengths (`sourceLengths`): opaque, each keyed by its system's address and its source node's
 * hash.
 */
export type Lengths = ReadonlyMap<string, number>;

/** A length's key: the system and the source node it was read from. */
const lengthKey = (system: string, sourceHash: string) => `${system}|${sourceHash}`;

/** Each model's source nodes' hashes, by the source's address: what `sourceLengths` keys a length by. */
const SourceHashes = new WeakMap<StaticModel, ReadonlyMap<string, string>>();

export interface Resolver {
  /** The model of a spec; `lengths` as `ResolveOptions.lengths`, for this call. */
  resolve(spec: Node, lengths?: Lengths): StaticModel;
}

/**
 * The lengths the opened sources' coverages give the model's systems of unknown length (letters, an alignment's
 * columns and rows from a source), by system address: what the next resolve takes (`ResolveOptions.lengths`). A
 * system is found in its source's coverage by its name or an alias; a sequence whose opened source names systems but
 * not its own is kept as unnamed (`seq.resolve.system.unnamed` at the next resolve).
 */
export function sourceLengths(model: StaticModel, coverage: (source: string) => Coverage | undefined): Lengths {
  const out = new Map<string, number>();
  const hashes = SourceHashes.get(model);
  for (const s of model.systems.values()) {
    const hash = s.source && hashes?.get(s.source);
    if (!hash) continue;
    const known = coverage(s.source!)?.systems;
    const info = known?.find((x) => isSystemInfo(x, s));
    if (info?.length !== undefined) out.set(lengthKey(s.address, hash), info.length);
    else if (!info && known?.length && s.name !== AlignmentSystems.rows && !s.aliases)
      out.set(lengthKey(s.address, hash), Unnamed);
  }
  return out;
}

/** Whether a coverage's system is `system`: by its name, or an alias of either (`SystemEntry.aliases`). */
export function isSystemInfo(info: SystemInfo, system: SystemEntry): boolean {
  return info.name === system.name || !!info.aliases?.includes(system.name) || !!system.aliases?.includes(info.name);
}

/** A length that says the system's source opened and names no system after it (never a length). */
const Unnamed = -1;

/** The lengths' part of a model's hash: none when empty, so a spec without sources keeps its spec's hash. */
function lengthsKey(lengths: Lengths | undefined): string {
  if (!lengths?.size) return '';
  return `|${[...lengths]
    .sort((a, b) => (a[0] < b[0] ? -1 : 1))
    .map(([k, v]) => `${k}=${v}`)
    .join(',')}`;
}

/** A path search's outcome and the entries (maps, units) of the systems it read: valid while those are the same. */
interface PathEntry {
  readonly found: { readonly path: Path } | { readonly why: 'none' | 'ambiguous' | 'limit' };
  readonly deps: ReadonlyMap<string, readonly (MapEntry | UnitEntry)[]>;
}

/** A resolve pass that keeps its entries between specs, by address and the hashes of the nodes each reads. */
export function createResolver(options: ResolveOptions = {}): Resolver {
  let entries = new Map<string, unknown>();
  let last: Model | undefined;
  return {
    resolve(spec, lengths = options.lengths) {
      if (last && last.hash === `${hashNode(spec)}${lengthsKey(lengths)}`) return last;
      const next = new Map<string, unknown>();
      const cached = <T>(key: string, make: () => T): T => {
        const v = (entries.has(key) ? entries.get(key) : make()) as T;
        next.set(key, v);
        return v;
      };
      last = build(spec, { ...options, ...(lengths ? { lengths } : {}) }, cached, last?.pathEntries);
      entries = next;
      return last;
    },
  };
}

/** The static model of a checked spec, with no cache. */
export function resolveSpec(spec: Node, options: ResolveOptions = {}): StaticModel {
  return build(spec, options, (_, make) => make());
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** The kinds a data reference to a system takes: a sequence, a molecule's residues, a unit of a complex. */
const SystemRefs: readonly string[] = SeqLanguage.kinds.map.params.from.kinds;
const SystemKinds: ReadonlySet<string> = new Set(SystemRefs);
/** What `pairs.a`/`b`, `ranges.target` and `map.to` take: a system or an alignment, by its columns (the description's). */
const ColumnsRefs: readonly string[] = SeqLanguage.kinds.map.params.to.kinds;

function build(spec: Node, options: ResolveOptions, cached: Cached, carried?: ReadonlyMap<string, PathEntry>): Model {
  const config: ResolveConfig = { ...DefaultResolveConfig, ...options.config };
  const issues: Issue[] = [];
  const data: Placed[] = [];
  const byAddress = new Map<string, Placed>();
  for (const c of spec.children ?? []) {
    if (c.kind !== 'data' || c.id === undefined) continue;
    for (const n of c.children ?? []) {
      if (n.id === undefined) continue;
      const placed = { node: n, address: printAddress([c.id, n.id]) };
      data.push(placed);
      byAddress.set(placed.address, placed);
    }
  }
  const inner = placeChildren(data, byAddress);
  const systems = new Map<string, SystemEntry>();
  const ref = (at: Placed, key: string, kinds: readonly string[]): Placed | undefined => {
    const a = exactAddress(p(at.node, key));
    const t = a === undefined ? undefined : byAddress.get(a);
    // a system's node that made no system (a unit of a refused complex, or of nothing held) is named, not taken
    if (t && kinds.includes(t.node.kind) && (SystemKinds.has(t.node.kind) ? systems.has(t.address) : true)) return t;
    if (t && kinds.includes(t.node.kind)) {
      issues.push(
        issueOf(ResolveIssues, 'seq.resolve.ref.system', 'The reference names a node that is no system', {
          address: at.address,
          path: ['params', key],
          hint: ResolveIssues.codes['seq.resolve.ref.system'].hint,
        }),
      );
      return undefined;
    }
    issues.push(
      issueOf(ResolveIssues, 'seq.resolve.ref.missing', 'The reference names no node of a kind it takes', {
        address: at.address,
        path: ['params', key],
        hint: ResolveIssues.codes['seq.resolve.ref.missing'].hint,
      }),
    );
    return undefined;
  };
  // a source-backed system's length (given for its source node as it is now) is part of what an entry over it reads
  const given: Lengths = options.lengths ?? new Map();
  const sourceHashes = new Map<string, string>();
  for (const at of data) if (at.node.kind === 'source') sourceHashes.set(at.address, hashNode(at.node));
  const lengths = new Map<string, number>();
  const unnamed = new Set<string>();
  if (given.size)
    for (const at of data) {
      if (at.node.kind !== 'sequence' && at.node.kind !== 'alignment' && at.node.kind !== 'molecule') continue;
      const src = exactAddress(p(at.node, 'source'));
      const hash = src === undefined ? undefined : sourceHashes.get(src);
      if (hash === undefined) continue;
      const systems = at.node.kind === 'alignment' ? [at.address, rowsSystemOf(at.address)] : [at.address];
      for (const system of systems) {
        const n = given.get(lengthKey(system, hash));
        if (n === Unnamed) unnamed.add(system);
        else if (n !== undefined) lengths.set(system, n);
      }
    }
  const key = (...at: Placed[]) => at.map((a) => `${hashNode(a.node)}#${lengths.get(a.address) ?? ''}`).join(':');
  const runtimes = options.runtimes ?? [];

  // sources first: what reads by region names them
  const sources = new Map<string, SourceBinding>();
  for (const at of data) {
    if (at.node.kind !== 'source') continue;
    const bound = cached(`${at.address}|${key(at)}`, () => bindSource(at, runtimes));
    issues.push(...bound.issues);
    sources.set(at.address, bound.binding);
  }

  for (const at of data) {
    if (at.node.kind !== 'sequence') continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    systems.set(
      at.address,
      cached(`${at.address}|${key(at)}`, () => sequenceEntry(at, src?.address, lengths.get(at.address))),
    );
    // its source opened and names no system after it: its length stays unknown and nothing over it lays out (a spec
    // fault: said, never a pending wait)
    if (unnamed.has(at.address)) {
      const name = systems.get(at.address)!.name;
      issues.push(
        issueOf(
          ResolveIssues,
          'seq.resolve.system.unnamed',
          `The sequence's source names no system ${name}: its length is unknown and nothing over it is laid out`,
          {
            address: at.address,
            path: ['params', p(at.node, 'name') === undefined ? 'source' : 'name'],
            hint: ResolveIssues.codes['seq.resolve.system.unnamed'].hint,
          },
        ),
      );
    }
  }

  const molecular = buildMolecular(data, inner, systems, issues, config, cached, key, ref, lengths);

  // a sequence's and a unit's numberings: inline labels read at once, a source's where its rows are
  const numberings = new Map<string, NumberingEntry>();
  for (const at of [...data, ...[...inner.values()].flat()]) {
    if (at.node.kind !== 'sequence' && (at.node.kind !== 'unit' || !systems.has(at.address))) continue;
    const system = systems.get(at.address)!;
    for (const c of at.node.children ?? []) {
      if (c.kind !== 'numbering' || c.id === undefined || c.id === 'index' || c.id === 'ordinal') continue;
      const placed = { node: c, address: printAddress([...addressIds(at.address), c.id]) };
      const src = p(c, 'source') === undefined ? undefined : ref(placed, 'source', ['source']);
      numberings.set(
        placed.address,
        cached(`${placed.address}|${key(placed, at)}|${system.name}|${system.length ?? ''}`, () =>
          numberingEntry(placed, system, src?.address, system.length ?? config.openQueryLength),
        ),
      );
    }
  }

  const annotations = new Map<string, AnnotationEntry>();
  // an alignment: its columns and rows systems, its cells an annotation on the columns; before what names its columns
  // (a map's `to`, a pairs node's ends, ranges)
  for (const at of data) {
    if (at.node.kind !== 'alignment') continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    const rowsAt = rowsSystemOf(at.address);
    const e = cached(`${at.address}|${key(at)}|${lengths.get(rowsAt) ?? ''}`, () =>
      alignmentEntry(at, src?.address, lengths.get(at.address), lengths.get(rowsAt)),
    );
    systems.set(e.columns.address, e.columns);
    systems.set(e.rows.address, e.rows);
    annotations.set(at.address, e.annotation);
  }

  const maps = new Map<string, MapEntry>();
  const describedMaps = new Map<string, DescribedMap>();
  const isSource = (a: string) => byAddress.get(a)?.node.kind === 'source';
  for (const at of data) {
    if (at.node.kind !== 'map') continue;
    // an end a source describes (a unit of its complex): the overlay's once it arrives, pending until then
    const dFrom = describedAddress(p(at.node, 'from'), isSource);
    const dTo = describedAddress(p(at.node, 'to'), isSource);
    if (dFrom !== undefined || dTo !== undefined) {
      const f = dFrom ?? ref(at, 'from', SystemRefs)?.address;
      const t = dTo ?? ref(at, 'to', ColumnsRefs)?.address;
      if (f !== undefined && t !== undefined)
        describedMaps.set(at.address, { address: at.address, node: at.node, from: f, to: t });
      continue;
    }
    const from = ref(at, 'from', SystemRefs);
    const to = ref(at, 'to', ColumnsRefs);
    if (!from || !to) continue;
    const entry = cached(`${at.address}|${key(at, from, to)}|${config.encoding}`, () =>
      mapEntry(at, systems.get(from.address)!, systems.get(to.address)!, config, molecular.units.get(to.address)),
    );
    maps.set(at.address, entry.map);
    issues.push(...entry.issues);
  }

  // pairs and grids: held whole on their two systems
  for (const at of data) {
    if (at.node.kind !== 'pairs' && at.node.kind !== 'grid') continue;
    const isPairs = at.node.kind === 'pairs';
    const kinds = isPairs ? ColumnsRefs : SystemRefs;
    const x = ref(at, isPairs ? 'a' : 'x', kinds);
    if (!x) continue;
    const yKey = isPairs ? 'b' : 'y';
    const y = p(at.node, yKey) === undefined ? x : ref(at, yKey, kinds);
    if (!y) continue;
    const entry = cached(`${at.address}|${key(at, x, y)}`, () =>
      heldEntry(at, systems.get(x.address)!, systems.get(y.address)!),
    );
    annotations.set(at.address, entry.annotation);
    issues.push(...entry.issues);
  }
  for (const at of data) {
    if (at.node.kind !== 'ranges' && at.node.kind !== 'values') continue;
    const target = ref(at, 'target', at.node.kind === 'ranges' ? ColumnsRefs : SystemRefs);
    if (!target) continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    const entry = cached(`${at.address}|${key(at, target)}`, () =>
      annotationEntry(at, systems.get(target.address)!, src?.address),
    );
    annotations.set(at.address, entry.annotation);
    issues.push(...entry.issues);
  }

  buildLinks(inner, molecular, systems, annotations, cached, key);

  const compounds = new Map<string, CompoundEntry>();
  for (const at of data)
    if (at.node.kind === 'compounds')
      for (const c of at.node.children ?? [])
        if (c.kind === 'compound' && c.id !== undefined) compounds.set(c.id, compoundEntry(c));

  const model = new Model(
    spec,
    config,
    systems,
    maps,
    sources,
    annotations,
    numberings,
    issues,
    molecular,
    compounds,
    { equivalents: options.equivalents, componentNames: options.componentNames },
    carried,
    given,
    describedMaps,
  );
  SourceHashes.set(model, sourceHashes);
  return model;
}

// the entries: model-entries.ts

// ---- the model ----

/** A codon map read both ways along itself: for closure. */
interface CodonPaths {
  readonly map: MapEntry;
  readonly down: Path;
  readonly up: Path;
}

class Model implements StaticModel {
  readonly hash: string;
  /** The path searches of this model, by `from`, `to`; the next resolve carries them. */
  readonly pathEntries = new Map<string, PathEntry>();
  private readonly results = new Map<string, PathResult>();
  private readonly edges = new Map<string, Edge[]>();
  /** The entries whose steps leave each system (maps, units): what a carried search is checked against. */
  private readonly edgeMaps = new Map<string, (MapEntry | UnitEntry)[]>();
  private readonly codons: CodonPaths[] = [];
  private readonly querySystems = new Map<string, QuerySystem>();
  private readonly queryAnnotations = new Map<string, QueryAnnotation | undefined>();
  private readonly numberingsBySystem = new Map<string, Map<string, NumberingEntry>>();
  readonly hasSourceNumberings: boolean;

  constructor(
    readonly spec: Node,
    readonly config: ResolveConfig,
    readonly systems: ReadonlyMap<string, SystemEntry>,
    readonly maps: ReadonlyMap<string, MapEntry>,
    readonly sources: ReadonlyMap<string, SourceBinding>,
    readonly annotations: ReadonlyMap<string, AnnotationEntry>,
    readonly numberings: ReadonlyMap<string, NumberingEntry>,
    readonly issues: readonly Issue[],
    readonly molecular: Molecular = EmptyMolecular,
    readonly compounds: ReadonlyMap<string, CompoundEntry> = new Map(),
    private readonly host: Pick<ResolveOptions, 'equivalents' | 'componentNames'> = {},
    private readonly carried?: ReadonlyMap<string, PathEntry>,
    lengths: Lengths = new Map(),
    readonly describedMaps: ReadonlyMap<string, DescribedMap> = new Map(),
  ) {
    this.hash = `${hashNode(spec)}${lengthsKey(lengths)}`;
    let fromSource = false;
    for (const n of numberings.values()) {
      if (n.source) fromSource = true;
      let of = this.numberingsBySystem.get(n.system);
      if (!of) this.numberingsBySystem.set(n.system, (of = new Map()));
      of.set(n.id, n);
    }
    this.hasSourceNumberings = fromSource;
    const add = (e: Edge, owner: MapEntry | UnitEntry) => {
      let list = this.edges.get(e.from);
      if (!list) {
        this.edges.set(e.from, (list = []));
        this.edgeMaps.set(e.from, []);
      }
      list.push(e);
      this.edgeMaps.get(e.from)!.push(owner);
    };
    // a unit's steps up to its molecule and down from it: free, under the copy rule (_spec/README.md "Resolve", Paths)
    for (const u of molecular.units.values()) for (const e of unitEdges(u)) add(e, u);
    for (const m of maps.values()) {
      const rank = m.priority ? { priority: m.priority } : {};
      const fwd: Edge = {
        map: m.address,
        from: m.from,
        to: m.to,
        backward: false,
        step: m.step,
        linear: m.linear,
        ...rank,
      };
      const back: Edge = {
        map: m.address,
        from: m.to,
        to: m.from,
        backward: true,
        step: -m.step as Step,
        linear: m.linear.inverse(),
        ...rank,
      };
      for (const e of [fwd, back]) add(e, m);
      if (m.step !== 0)
        this.codons.push({ map: m, down: composeEdges(m.to, [back])!, up: composeEdges(m.from, [fwd])! });
    }
  }

  edgesFrom(system: string): readonly Edge[] {
    return this.edges.get(system) ?? [];
  }

  /** A carried search is valid while every system it read has the same map entries. */
  private valid(entry: PathEntry): boolean {
    for (const [system, maps] of entry.deps) {
      const now = this.edgeMaps.get(system) ?? [];
      if (now.length !== maps.length) return false;
      for (let i = 0; i < now.length; i++) if (now[i] !== maps[i]) return false;
    }
    return true;
  }

  private entry(from: string, to: string, via: readonly string[]): PathEntry {
    const k = pathKey(from, to, via);
    let e = this.pathEntries.get(k);
    if (e) return e;
    const old = this.carried?.get(k);
    if (old && this.valid(old)) e = old;
    else {
      const deps = new Map<string, readonly (MapEntry | UnitEntry)[]>();
      // a use's `via` leaves out the named maps' rivals; what the search read is every map, rivals too
      const graph = {
        edgesFrom: (s: string) => {
          if (!deps.has(s)) deps.set(s, this.edgeMaps.get(s) ?? []);
          return this.edges.get(s) ?? [];
        },
      };
      const r = findPath(via.length ? withoutRivals(graph, viaRivals(via, this.maps.values())) : graph, from, to, {
        maxSteps: this.config.maxPathSteps,
        maxCandidates: this.config.maxPathCandidates,
      });
      e = { found: r.ok ? { path: r.path } : { why: r.why }, deps };
    }
    this.pathEntries.set(k, e);
    return e;
  }

  path(from: string, to: string, via?: readonly string[]): PathResult {
    const named = via?.length ? viaList(via) : [];
    const k = pathKey(from, to, named);
    let r = this.results.get(k);
    if (r) return r;
    r = viaPath(named, (v) => this.searched(from, to, v));
    this.results.set(k, r);
    return r;
  }

  /** One search's result, a pair and `via` set: a path read the other way is the inverse of the one found. */
  private searched(from: string, to: string, named: readonly string[]): PathResult {
    const k = pathKey(from, to, named);
    const back = this.pathEntries.get(pathKey(to, from, named));
    const e = back && !this.pathEntries.has(k) ? back : this.entry(from, to, named);
    if ('path' in e.found) return { path: e === back ? inversePath(e.found.path) : e.found.path };
    return { issue: pathIssue(e.found.why, from, to) };
  }

  project(loc: Location, to: string, options: { close?: boolean; via?: readonly string[] } = {}): Projection {
    const src = options.close ? this.close(loc) : loc;
    const issues: Issue[] = [];
    const parts: { system: string; ranges: number[] }[] = [];
    // the atoms group map at the ends (never a path step): onto an atoms system through its residues
    const toAtoms = this.molecular.atoms.get(to);
    const end = toAtoms ? toAtoms.of : to;
    for (const part of src.parts) {
      const fromAtoms = this.molecular.atoms.get(part.system);
      const from = fromAtoms ? fromAtoms.of : part.system;
      const ranges = fromAtoms ? collapseAtoms(fromAtoms, part.ranges) : part.ranges;
      if (part.system === to) {
        parts.push({ system: to, ranges: Array.from(part.ranges) });
        continue;
      }
      const r = this.path(from, end, options.via);
      if (!r.path) {
        if (r.issue) issues.push(r.issue);
        continue;
      }
      const image = projectRanges(r.path, ranges);
      parts.push({ system: to, ranges: toAtoms ? expandAtoms(toAtoms, image) : image });
    }
    return { location: location({ parts }), issues };
  }

  close(loc: Location, among?: readonly string[]): Location {
    const allowed = among ? new Set(among) : undefined;
    let out = loc;
    for (const { map, down, up } of this.codons) {
      if (allowed && (!allowed.has(map.from) || !allowed.has(map.to))) continue;
      const nt = partOf(loc, map.to);
      if (!nt) continue;
      const residues = projectRanges(down, nt.ranges);
      const codons = projectRanges(up, residues);
      out = unionLocations(
        out,
        location({
          parts: [
            { system: map.from, ranges: residues },
            { system: map.to, ranges: codons },
          ],
        }),
      );
    }
    return out;
  }

  numberingOf(system: string, id: string): NumberingEntry | undefined {
    return this.numberingsBySystem.get(system)?.get(id);
  }

  // ---- the query forms' model ----

  system(address: string): QuerySystem | undefined {
    let s = this.querySystems.get(address);
    if (s) return s;
    const e = this.systems.get(address);
    if (!e) return undefined;
    // its declared numberings: the inline ones answer here; a source's labels are tier 2's (`withSourceNumberings`)
    const own = this.numberingsBySystem.get(address);
    const ids = own ? [...own.keys()] : [];
    const parent = e.parent && this.system(e.parent.system);
    s = {
      address,
      length: e.length ?? this.config.openQueryLength,
      cyclic: e.cyclic,
      numbering: (id) => own?.get(id)?.numbering,
      numberingIds: () => ids,
      // a source's labels are tier 2's: here they never arrive, so a form naming them waits
      numberingPending: (id) => own?.get(id)?.source !== undefined,
      // a unit's molecule numbers its positions where the unit has no numbering of an id
      ...(parent ? { parent: { system: parent, start: e.parent!.start } } : {}),
    };
    this.querySystems.set(address, s);
    return s;
  }

  annotation(address: string): QueryAnnotation | undefined {
    if (this.queryAnnotations.has(address)) return this.queryAnnotations.get(address);
    const e = this.annotations.get(address);
    const a =
      e?.shape === 'ranges'
        ? queryAnnotation(e, this.systems.get(e.target), this.config)
        : e?.shape === 'pairs'
          ? pairsQueryAnnotation(e, (s) => this.systems.get(s), this.config)
          : undefined;
    this.queryAnnotations.set(address, a);
    return a;
  }

  /** Dynamic nodes are the overlay's: none in the static model. */
  node(): Node | undefined {
    return undefined;
  }

  /** An inline alignment's row by its label (its rows held in the spec); a source's rows are not held here. */
  row(rows: string, key: Atom): number | undefined {
    for (const a of this.annotations.values()) {
      if (a.shape !== 'cells' || a.rows !== rows || a.access?.kind !== 'inline') continue;
      const label = (a.access.source.value as { label?: unknown } | undefined)?.label;
      if (!Array.isArray(label)) return undefined;
      const i = label.indexOf(String(key));
      return i < 0 ? undefined : i;
    }
    return undefined;
  }

  /** Chains are a complex's units (`units`): no chain addressed by its form in the static model. */
  chain(): Node | undefined {
    return undefined;
  }

  modified(system: string): ModifiedIndex | undefined {
    return modifiedOf(this.molecular, system);
  }

  components(system: string) {
    return componentsOf(this.molecular, this.systems.get(system), (s) => this.systems.get(s));
  }

  compoundClass(code: string): string | undefined {
    return this.compound(code)?.class;
  }

  compound(code: string): Omit<CompoundEntry, 'id'> | undefined {
    return this.compounds.get(code) ?? this.host.componentNames?.(code);
  }

  equivalents(identifier: string): readonly string[] {
    return this.host.equivalents?.(identifier) ?? [];
  }

  annotationsOn(system: string): readonly string[] {
    const out: string[] = [];
    for (const a of this.annotations.values()) if (a.shape === 'ranges' && a.target === system) out.push(a.address);
    return out;
  }

  units(complex: string) {
    return unitsOf(this.molecular, complex);
  }

  molecule(address: string) {
    return moleculeOf(this.molecular, address);
  }

  atoms(system: string) {
    return atomsOf(this.molecular, system);
  }

  complexes(owner: string) {
    return complexesUnder(this.molecular, owner);
  }
}

/** A path's issue, at the system it starts from, with both systems as its places. */
export function pathIssue(why: 'none' | 'ambiguous' | 'limit', from: string, to: string): Issue {
  const code =
    why === 'none'
      ? 'seq.resolve.path.none'
      : why === 'limit'
        ? 'seq.resolve.path.limit'
        : 'seq.resolve.path.ambiguous';
  const message =
    why === 'none'
      ? 'No path of maps joins the systems'
      : why === 'limit'
        ? 'More shortest paths join the systems than are compared'
        : 'Several shortest paths join the systems and map them differently';
  return {
    ...issueOf(ResolveIssues, code, message, { address: from, hint: ResolveIssues.codes[code].hint }),
    places: [{ address: from }, { address: to }],
  };
}

/** Inline items answer at once; a source's are pending over what was asked (no data in tier 1). */
export function queryAnnotation(
  e: AnnotationEntry,
  target: SystemEntry | undefined,
  config: ResolveConfig,
): QueryAnnotation {
  const length = target?.length ?? config.openQueryLength;
  return {
    address: e.address,
    target: e.target,
    fields: e.fields,
    items(over?: Location): ItemsRead {
      const part = over && partOf(over, e.target);
      // a location on none of its system asks nothing of it (no pending on another system)
      if (over && !part) return { items: [] };
      if (!e.items)
        return {
          items: [],
          pending: part ? rangesLocation(e.target, part.ranges) : rangesLocation(e.target, [0, length]),
        };
      if (!part) return { items: e.items };
      return { items: e.items.filter((it) => overlapsRanges(part.ranges, it.start, it.end, length)) };
    },
  };
}

/**
 * @internal A pairs annotation as the query forms read it (_spec/README.md "Queries", Pairs answer `items`): a pair one item, a
 * row per end on its own system (`QueryItem.system`), so `(items :overlapping <location>)` gives the pairs with an end
 * in the location, on either system. Its fields: `type`, `label` and the rows' other category columns. Held whole: the
 * rows are made once, on first read.
 */
export function pairsQueryAnnotation(
  e: AnnotationEntry,
  systemOf: (address: string) => SystemEntry | undefined,
  config: ResolveConfig,
): QueryAnnotation {
  const a = e.target;
  const b = e.ySystem ?? e.target;
  const lengthOf = (s: string) => systemOf(s)?.length ?? config.openQueryLength;
  // a complex's links: each end on its row's unit
  const units = e.units;
  const aUnit = e.pairs?.aUnit;
  const bUnit = e.pairs?.bUnit;
  let rows: QueryItem[] | undefined;
  const all = (): QueryItem[] => {
    if (rows) return rows;
    rows = [];
    const d = e.pairs;
    if (!d) return rows;
    const extra = Object.entries(d.extra ?? {});
    for (let k = 0; k < d.n; k++) {
      const fields: Record<string, Atom> = { type: d.kind[k] < 0 ? null : d.kinds[d.kind[k]] };
      if (d.labels) fields.label = d.labels[k] ?? null;
      for (const [name, col] of extra) fields[name] = col[k] ?? null;
      const key = d.keys[k];
      const sa = units && aUnit ? units[aUnit[k]] : a;
      const sb = units && bUnit ? units[bUnit[k]] : b;
      rows.push({ key, start: d.aStart[k], end: d.aEnd[k], system: sa, fields });
      rows.push({ key, start: d.bStart[k], end: d.bEnd[k], system: sb, fields });
    }
    return rows;
  };
  const out: QueryAnnotation & { readonly systems: readonly string[] } = {
    address: e.address,
    target: a,
    fields: e.fields,
    // its ends' systems, read by the query forms (internal: `QueryAnnotation` names no such field)
    systems: units ?? (b === a ? [a] : [a, b]),
    items(over?: Location): ItemsRead {
      if (!over) return { items: all() };
      return {
        items: all().filter((it) => {
          const p = partOf(over, it.system!);
          return !!p && overlapsRanges(p.ranges, it.start, it.end, lengthOf(it.system!));
        }),
      };
    },
  };
  return out;
}
