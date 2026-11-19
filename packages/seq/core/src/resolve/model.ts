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
  CategoryColumn,
  categoryColumn,
  exactAddress,
  hashNode,
  type Issue,
  issueOf,
  isColumn,
  type Node,
  printAddress,
  type Value,
} from '@molstar/core/repr/tree';
import {
  AlignmentRuntime,
  AlignmentSystems,
  type Coverage,
  type FormatRuntime,
  inlineSeq,
  type RangesField,
  type RangesList,
  SeqRuntime,
  type Shape,
  type SourceEntry,
  seqToValue,
  type Tables,
} from '@molstar/seqx-core/data';
import type { SeqKinds } from '@molstar/seqx-core/lang';
import {
  type ItemsRead,
  type Location,
  location,
  overlapsRanges,
  partOf,
  type QueryAnnotation,
  type QueryItem,
  type QueryModel,
  type QuerySystem,
  rangesLocation,
  unionLocations,
} from '@molstar/seqx-core/query';
import {
  codingLength,
  codonLinear,
  codonRuns,
  DefaultGeneticCode,
  type GeneticCode,
  GeneticCodes,
  hasUnsplitExon,
} from './codon.js';
import { type GridData, type PairsData, readGrid, readPairs } from './held.js';
import { ResolveIssues } from './issues.js';
import { composeLinear, type Encoding, expand, ImplicitLinear, isEmptyLinear, type Linear, toRuns } from './linear.js';
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
}

export const DefaultResolveConfig: ResolveConfig = {
  maxPathSteps: 4,
  maxPathCandidates: 16,
  encoding: 'cheapest',
  openQueryLength: OpenLength,
  maxNestedDepth: 4,
};

/** A coordinate system the spec declares: a sequence, its length when its letters are inline. */
export interface SystemEntry {
  readonly address: string;
  readonly kind: 'sequence';
  readonly alphabet: 'protein' | 'dna' | 'rna';
  /** Undefined: unknown until its source is opened (pending, never a guess). */
  readonly length?: number;
  readonly cyclic: boolean;
  readonly letters?: string;
  /** The source of its letters. */
  readonly source?: string;
  /** The name its source knows it by: its `name`, else its id. */
  readonly name: string;
  /** Its numberings: the built-ins (`index`, `ordinal`) until the language declares others. */
  readonly numberings: readonly string[];
  /** How its letters are read. */
  readonly access: AccessPlan;
}

/** A map, resolved: its kind, its encoding, and its linear part (fine places of the residue side for a codon map). */
export interface MapEntry {
  readonly address: string;
  readonly kind: 'offset' | 'codon';
  readonly from: string;
  readonly to: string;
  readonly encoding: Encoding;
  readonly step: Step;
  readonly linear: Linear;
  /** A codon map's genetic code: its `table`, 1 when absent; undefined when that code is not built. */
  readonly code?: GeneticCode;
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
 * How a node's data is read: inline data as a `seq` source of one tile (the whole system), or a source of the spec by
 * region. `system` is the name the source knows the system by.
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
}

/**
 * An alignment's rows system's address: its node's, then `(rows)` (rows are a system, EC-35: positions are rows in
 * data order); its columns system is the node's own address.
 */
export function rowsSystemOf(alignment: string): string {
  return `${alignment}/(rows)`;
}

/** A path between two systems, or why there is none. */
export interface PathResult {
  readonly path?: Path;
  readonly issue?: Issue;
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
  readonly issues: readonly Issue[];
  /** The source-backed systems' lengths it was resolved with (`ResolveOptions.lengths`). */
  readonly lengths: Lengths;
  /** The shortest unambiguous path, composed once and kept (across resolves while the maps its search read hold). */
  path(from: string, to: string): PathResult;
  /** The steps from a system along the spec's maps, both ways: what a path search reads (the overlay adds its own). */
  edgesFrom(system: string): readonly Edge[];
  /**
   * A Location carried onto a system: each part along its path (rects are not carried yet); a part with no path is
   * left out with the path's issue. With `close`, the location is first closed (`close`).
   */
  project(loc: Location, to: string, options?: { close?: boolean }): Projection;
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
   * The lengths of source-backed systems, by system address, that their sources' coverages gave once opened
   * (`sourceLengths`): such a system then resolves as one with inline letters of that length (_spec/README.md
   * "Resolve", Systems).
   */
  readonly lengths?: Lengths;
}

/** Source-backed systems' lengths by system address (`sourceLengths`). */
export type Lengths = ReadonlyMap<string, number>;

export interface Resolver {
  /** The model of a spec; `lengths` as `ResolveOptions.lengths`, for this call. */
  resolve(spec: Node, lengths?: Lengths): StaticModel;
}

/**
 * The lengths the opened sources' coverages give the model's systems of unknown length (letters, an alignment's
 * columns and rows from a source), by system address: what the next resolve takes (`ResolveOptions.lengths`). A
 * system is found in its source's coverage by its name or an alias.
 */
export function sourceLengths(model: StaticModel, coverage: (source: string) => Coverage | undefined): Lengths {
  const out = new Map<string, number>();
  for (const s of model.systems.values()) {
    if (!s.source) continue;
    const info = coverage(s.source)?.systems.find((x) => x.name === s.name || x.aliases?.includes(s.name));
    if (info?.length !== undefined) out.set(s.address, info.length);
  }
  return out;
}

/** Whether two lengths say the same. */
export function lengthsEqual(a: Lengths, b: Lengths): boolean {
  if (a.size !== b.size) return false;
  for (const [k, v] of a) if (b.get(k) !== v) return false;
  return true;
}

/** The lengths' part of a model's hash: none when empty, so a spec without sources keeps its spec's hash. */
function lengthsKey(lengths: Lengths | undefined): string {
  if (!lengths?.size) return '';
  return `|${[...lengths]
    .sort((a, b) => (a[0] < b[0] ? -1 : 1))
    .map(([k, v]) => `${k}=${v}`)
    .join(',')}`;
}

/** A path search's outcome and the map entries of the systems it read: valid while those are the same. */
interface PathEntry {
  readonly found: { readonly path: Path } | { readonly why: 'none' | 'ambiguous' | 'limit' };
  readonly deps: ReadonlyMap<string, readonly MapEntry[]>;
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

type Cached = <T>(key: string, make: () => T) => T;

interface Placed {
  readonly node: Node;
  readonly address: string;
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

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
  const ref = (at: Placed, key: string, kinds: readonly string[]): Placed | undefined => {
    const a = exactAddress(p(at.node, key));
    const t = a === undefined ? undefined : byAddress.get(a);
    if (t && kinds.includes(t.node.kind)) return t;
    issues.push(
      issueOf(ResolveIssues, 'seq.resolve.ref.missing', 'The reference names no node of a kind it takes', {
        address: at.address,
        path: ['params', key],
        hint: ResolveIssues.codes['seq.resolve.ref.missing'].hint,
      }),
    );
    return undefined;
  };
  // a source-backed system's length is part of what an entry over it reads
  const lengths: Lengths = options.lengths ?? new Map();
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

  const systems = new Map<string, SystemEntry>();
  for (const at of data) {
    if (at.node.kind !== 'sequence') continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    systems.set(
      at.address,
      cached(`${at.address}|${key(at)}`, () => sequenceEntry(at, src?.address, lengths.get(at.address))),
    );
  }

  const maps = new Map<string, MapEntry>();
  for (const at of data) {
    if (at.node.kind !== 'map') continue;
    const from = ref(at, 'from', ['sequence']);
    const to = ref(at, 'to', ['sequence']);
    if (!from || !to) continue;
    const entry = cached(`${at.address}|${key(at, from, to)}|${config.encoding}`, () =>
      mapEntry(at, systems.get(from.address)!, systems.get(to.address)!, config),
    );
    maps.set(at.address, entry.map);
    issues.push(...entry.issues);
  }

  const annotations = new Map<string, AnnotationEntry>();
  // an alignment: its columns and rows systems, its cells an annotation on the columns
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
  // pairs and grids: held whole on their two systems
  for (const at of data) {
    if (at.node.kind !== 'pairs' && at.node.kind !== 'grid') continue;
    const isPairs = at.node.kind === 'pairs';
    const x = ref(at, isPairs ? 'a' : 'x', ['sequence']);
    if (!x) continue;
    const yKey = isPairs ? 'b' : 'y';
    const y = p(at.node, yKey) === undefined ? x : ref(at, yKey, ['sequence']);
    if (!y) continue;
    const entry = cached(`${at.address}|${key(at, x, y)}`, () =>
      heldEntry(at, systems.get(x.address)!, systems.get(y.address)!),
    );
    annotations.set(at.address, entry.annotation);
    issues.push(...entry.issues);
  }
  for (const at of data) {
    if (at.node.kind !== 'ranges' && at.node.kind !== 'values') continue;
    const target = ref(at, 'target', ['sequence']);
    if (!target) continue;
    const src = p(at.node, 'source') === undefined ? undefined : ref(at, 'source', ['source']);
    const entry = cached(`${at.address}|${key(at, target)}`, () =>
      annotationEntry(at, systems.get(target.address)!, src?.address),
    );
    annotations.set(at.address, entry.annotation);
    issues.push(...entry.issues);
  }

  return new Model(spec, config, systems, maps, sources, annotations, issues, carried, lengths);
}

// ---- entries ----

// the transform's params (the records runtime's: _spec/README.md "Data", Records) ride on the entry
const SourceKeys = [
  'uri',
  'value',
  'index',
  'fields',
  'assembly',
  'records',
  'letters',
  'system',
  'base',
  'ends',
  'types',
] as const;

/** A source's entry from its params, bound to the first runtime that reads its format; none is an issue. */
function bindSource(at: Placed, runtimes: readonly FormatRuntime[]): { binding: SourceBinding; issues: Issue[] } {
  const format = p(at.node, 'format');
  const entry: Record<string, unknown> = { address: at.address, format: typeof format === 'string' ? format : '' };
  for (const k of SourceKeys) if (p(at.node, k) !== undefined) entry[k] = p(at.node, k);
  const bound = entry as unknown as BoundSourceEntry;
  const place = (code: 'seq.resolve.source.format' | 'seq.resolve.source.runtime') => ({
    address: at.address,
    hint: ResolveIssues.codes[code].hint,
  });
  if (typeof format !== 'string' || !format)
    return {
      binding: { address: at.address, entry: bound },
      issues: [
        issueOf(
          ResolveIssues,
          'seq.resolve.source.format',
          'The source names no format',
          place('seq.resolve.source.format'),
        ),
      ],
    };
  const runtime = runtimes.find((r) => r.formats.includes(format));
  if (runtime) return { binding: { address: at.address, entry: bound, runtime }, issues: [] };
  return {
    binding: { address: at.address, entry: bound },
    issues: [
      issueOf(
        ResolveIssues,
        'seq.resolve.source.runtime',
        'No runtime reads the format of the source',
        place('seq.resolve.source.runtime'),
      ),
    ],
  };
}

/** Inline tables as a `seq` source of one tile over the whole system, named by the system's address. */
function inlinePlan(
  address: string,
  system: { readonly address: string; readonly length?: number; readonly cyclic: boolean },
  shape: Shape,
  tables: Tables,
): AccessPlan {
  const name = system.address;
  const value = seqToValue(inlineSeq({ name, length: system.length, cyclic: system.cyclic }, tables));
  return { kind: 'inline', source: { address, format: 'seq', value }, runtime: SeqRuntime, system: name, shape };
}

/** A source read by region, the system by the name the source knows it by. */
function sourcePlan(source: string, system: string, shape: Shape): AccessPlan {
  return { kind: 'source', source, system, shape };
}

function sequenceEntry(at: Placed, source: string | undefined, known: number | undefined): SystemEntry {
  const letters = p(at.node, 'letters');
  const name = p(at.node, 'name');
  const base = {
    address: at.address,
    kind: 'sequence' as const,
    alphabet: p(at.node, 'alphabet') as SystemEntry['alphabet'],
    cyclic: p(at.node, 'cyclic') === true,
    name: typeof name === 'string' ? name : (at.node.id ?? ''),
    numberings: ['index', 'ordinal'],
  };
  if (typeof letters === 'string') {
    const tables: Tables = { letters: { form: 'letters', start: 0, letters } };
    const length = letters.length;
    return { ...base, length, letters, access: inlinePlan(at.address, { ...base, length }, 'letters', tables) };
  }
  return {
    ...base,
    // the length its source's coverage gave (`ResolveOptions.lengths`), else unknown
    ...(source && known !== undefined ? { length: known } : {}),
    ...(source ? { source } : {}),
    access: sourcePlan(source ?? '', base.name, 'letters'),
  };
}

/** An alignment's two systems and its cells, inline (read through `alignment`) or from a source by region. */
function alignmentEntry(
  at: Placed,
  source: string | undefined,
  knownColumns?: number,
  knownRows?: number,
): { columns: SystemEntry; rows: SystemEntry; annotation: AnnotationEntry } {
  const rows = p(at.node, 'rows') as { label?: Value; letters?: Value } | undefined;
  const alphabet = p(at.node, 'alphabet') as SystemEntry['alphabet'];
  const reference = p(at.node, 'reference');
  const letters = Array.isArray(rows?.letters) ? (rows.letters as string[]) : undefined;
  const label = Array.isArray(rows?.label) ? (rows.label as string[]) : undefined;
  // the table's other columns are the rows' own fields (_spec/README.md "Render", The rows layer, Row columns)
  const fields = rowsFields(rows as Record<string, Value> | undefined);
  const access: AccessPlan =
    letters && label
      ? {
          kind: 'inline',
          source: {
            address: at.address,
            format: 'alignment',
            value: {
              label,
              letters,
              ...(fields ? { fields } : {}),
              ...(firstOf(rows as Record<string, Value>, letters.length) ?? {}),
            } as unknown as Value,
          },
          runtime: AlignmentRuntime,
          system: AlignmentSystems.rows,
          shape: 'cells',
        }
      : sourcePlan(source ?? '', AlignmentSystems.rows, 'cells');
  const system = (address: string, name: string, length: number | undefined): SystemEntry => ({
    address,
    kind: 'sequence',
    alphabet,
    ...(length !== undefined ? { length } : {}),
    cyclic: false,
    name,
    numberings: ['index', 'ordinal'],
    ...(source ? { source } : {}),
    access,
  });
  const rowsAddress = rowsSystemOf(at.address);
  return {
    columns: system(
      at.address,
      AlignmentSystems.columns,
      letters ? (letters[0]?.length ?? 0) : source ? knownColumns : undefined,
    ),
    rows: system(rowsAddress, AlignmentSystems.rows, letters ? letters.length : source ? knownRows : undefined),
    annotation: {
      address: at.address,
      shape: 'cells',
      target: at.address,
      fields: ['label', 'first', ...Object.keys(fields ?? {})],
      access,
      rows: rowsAddress,
      ...(typeof reference === 'string' ? { reference } : {}),
    },
  };
}

/** An inline rows table's `first` column, as the rows' first residue numbers, when it is one number a row. */
function firstOf(rows: Record<string, Value> | undefined, n: number): { first: number[] } | undefined {
  const v = rows?.first;
  const list =
    isColumn(v) && !(v instanceof CategoryColumn)
      ? Array.from(v as ArrayLike<number>)
      : Array.isArray(v)
        ? v
        : undefined;
  return list && list.length === n && list.every((x) => typeof x === 'number')
    ? { first: list as number[] }
    : undefined;
}

/** An inline rows table's columns past `label`, `letters` and `first`, a plain list each; undefined when none. */
function rowsFields(rows: Record<string, Value> | undefined): Record<string, unknown[]> | undefined {
  let out: Record<string, unknown[]> | undefined;
  for (const key of Object.keys(rows ?? {})) {
    // `first` is the rows' first residue numbers (EC-39), never a field
    if (key === 'label' || key === 'letters' || key === 'first') continue;
    const v = rows![key];
    const list =
      v instanceof CategoryColumn
        ? Array.from(v.codes, (k) => (k < 0 ? null : v.dict[k]))
        : isColumn(v)
          ? Array.from(v as ArrayLike<number>)
          : Array.isArray(v)
            ? v
            : undefined;
    if (list) (out ??= Object.create(null) as Record<string, unknown[]>)[key] = list;
  }
  return out;
}

type MapCode =
  | 'seq.resolve.map.empty'
  | 'seq.resolve.map.codon_length'
  | 'seq.resolve.map.table'
  | 'seq.resolve.length.unknown';

function mapEntry(
  at: Placed,
  from: SystemEntry,
  to: SystemEntry,
  config: ResolveConfig,
): { map: MapEntry; issues: Issue[] } {
  const issues: Issue[] = [];
  const report = (code: MapCode, message: string) =>
    issues.push(issueOf(ResolveIssues, code, message, { address: at.address, hint: ResolveIssues.codes[code].hint }));
  const method = p(at.node, 'method');
  let linear: Linear;
  let step: Step = 0;
  let code: GeneticCode | undefined;
  if (method === 'codon') {
    step = 1;
    const exons = p(at.node, 'exons') as unknown as { start: ArrayLike<number>; end: ArrayLike<number> };
    const strand = p(at.node, 'strand') === '-' ? '-' : '+';
    const size = to.cyclic ? to.length : undefined;
    const unsplit = hasUnsplitExon(exons.start, exons.end, size);
    if (unsplit) report('seq.resolve.length.unknown', 'An exon crosses the origin of a gene of unknown length');
    const coding = codingLength(codonRuns(exons.start, exons.end, strand, size));
    const n = from.length;
    // a coding sequence missing a piece is not judged until it is whole
    if (!unsplit && (coding % 3 !== 0 || (n !== undefined && coding !== 3 * n && coding !== 3 * n + 3)))
      report('seq.resolve.map.codon_length', 'The exons are not three nucleotides a residue');
    const table = p(at.node, 'table');
    code = GeneticCodes[typeof table === 'number' ? table : DefaultGeneticCode];
    if (!code) report('seq.resolve.map.table', 'The genetic code of the map is not built');
    linear = codonLinear(exons.start, exons.end, strand, size);
    // what lies past either sequence maps nowhere (EC-21)
    if (to.length !== undefined) linear = composeLinear(linear, new ImplicitLinear(0, 1, 0, to.length));
    if (n !== undefined) linear = composeLinear(new ImplicitLinear(0, 1, 0, 3 * n), linear);
  } else {
    const by = Number(p(at.node, 'by') ?? 0);
    const lf = from.length ?? Infinity;
    if (to.cyclic && to.length !== undefined && Number.isFinite(lf))
      linear = new ImplicitLinear(by, 1, 0, lf, to.length);
    else linear = new ImplicitLinear(by, 1, Math.max(0, -by), Math.min(lf, (to.length ?? Infinity) - by));
  }
  if (isEmptyLinear(linear)) report('seq.resolve.map.empty', 'The map pairs no position');
  if (config.encoding !== 'cheapest' && linear.runs())
    linear = config.encoding === 'runs' ? toRuns(linear) : expand(linear, linear.lo, linear.hi);
  const map: MapEntry = {
    address: at.address,
    kind: method === 'codon' ? 'codon' : 'offset',
    from: from.address,
    to: to.address,
    encoding: linear.encoding,
    step,
    linear,
    ...(code ? { code } : {}),
  };
  return { map, issues };
}

/** A `map` node resolved (the overlay's item maps are built as one): its entry and its issues. */
export function resolveMapNode(
  address: string,
  node: Node,
  from: SystemEntry,
  to: SystemEntry,
  config: ResolveConfig,
): { map: MapEntry; issues: Issue[] } {
  return mapEntry({ node, address }, from, to, config);
}

/** A `ranges` or `values` node resolved on its target (the overlay's dynamic annotations are resolved so). */
export function resolveAnnotationNode(
  address: string,
  node: Node,
  target: SystemEntry,
  source: string | undefined,
): { annotation: AnnotationEntry; issues: Issue[] } {
  return annotationEntry({ node, address }, target, source);
}

/** A `pairs` or `grid` node on its two systems, held whole. */
function heldEntry(at: Placed, x: SystemEntry, y: SystemEntry): { annotation: AnnotationEntry; issues: Issue[] } {
  const base = { address: at.address, target: x.address, ySystem: y.address };
  if (at.node.kind === 'pairs') {
    const pairs = readPairs(at.node);
    const fields = ['type', ...(pairs.labels ? ['label'] : []), ...Object.keys(pairs.extra ?? {})];
    return { annotation: { ...base, shape: 'pairs', fields, pairs }, issues: [] };
  }
  const { grid, issues } = readGrid(at.address, at.node, x.length, y.length);
  return { annotation: { ...base, shape: 'grid', fields: ['value'], ...(grid ? { grid } : {}) }, issues };
}

const Strands: Readonly<Record<string, number>> = { '+': 1, '-': -1 };

function annotationEntry(
  at: Placed,
  target: SystemEntry,
  source: string | undefined,
): { annotation: AnnotationEntry; issues: Issue[] } {
  const shape = at.node.kind as 'ranges' | 'values';
  const base = { address: at.address, shape, target: target.address };
  const bySource = (s: Shape) => ({
    annotation: { ...base, fields: [], access: sourcePlan(source ?? '', target.name, s) },
    issues: [],
  });
  if (shape === 'values') {
    const values = p(at.node, 'values');
    if (!isColumn(values)) return bySource('values');
    const tables: Tables = { values: { form: 'positions', start: 0, value: values as Float32Array } };
    return {
      annotation: { ...base, fields: [], access: inlinePlan(at.address, target, 'values', tables) },
      issues: [],
    };
  }
  const rows = p(at.node, 'rows') as Record<string, Value> | undefined;
  if (!rows) return bySource('ranges');
  const { items, table, fields, wrappedUnknown } = inlineRanges(rows, target);
  const issues = wrappedUnknown
    ? [
        issueOf(ResolveIssues, 'seq.resolve.length.unknown', 'A row wraps the origin of a sequence of unknown length', {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.length.unknown'].hint,
        }),
      ]
    : [];
  return {
    annotation: { ...base, fields, items, access: inlinePlan(at.address, target, 'ranges', { ranges: table }) },
    issues,
  };
}

/**
 * Inline rows as items and a ranges table: a row that wraps a cyclic origin is two rows of one item (EC-9, EC-10); on
 * a sequence of unknown length it is left out, and said.
 */
/** The columns of inline ranges the format reads itself; the others are a table's custom columns. */
const KnownColumns = new Set(['start', 'end', 'id', 'label', 'category', 'strand']);

function inlineRanges(rows: Record<string, Value>, target: SystemEntry) {
  const start = rows.start as Int32Array;
  const end = rows.end as Int32Array;
  const id = rows.id as { get(i: number): string | null } | undefined;
  const strand = rows.strand as { get(i: number): string | null } | undefined;
  const label = Array.isArray(rows.label) ? (rows.label as (string | null)[]) : undefined;
  const category = rows.category as { get(i: number): string | null } | undefined;
  const fieldNames = [label && 'label', category && 'category', strand && 'strand'].filter((f): f is string => !!f);
  const items: QueryItem[] = [];
  const s: number[] = [];
  const e: number[] = [];
  const rowOf: number[] = [];
  const length = target.length;
  let wrappedUnknown = false;
  for (let i = 0; i < start.length; i++) {
    const fields: Record<string, Atom> = {};
    if (label) fields.label = label[i] ?? null;
    if (category) fields.category = category.get(i);
    if (strand) fields.strand = strand.get(i);
    const key = id?.get(i) ?? i;
    let pieces: [number, number][] = [[start[i], end[i]]];
    if (start[i] > end[i]) {
      if (length === undefined) {
        wrappedUnknown = true;
        continue;
      }
      pieces = [
        [start[i], length],
        [0, end[i]],
      ];
    }
    for (const [a, b] of pieces) {
      items.push({ key, start: a, end: b, ...(fieldNames.length ? { fields } : {}) });
      s.push(a);
      e.push(b);
      rowOf.push(i);
    }
  }
  const pick = <T>(f: (i: number) => T) => rowOf.map(f);
  const columns: RangesField[] = [];
  if (label) columns.push({ name: 'label', column: categoryColumn(pick((i) => label[i] ?? null)) });
  if (category) columns.push({ name: 'category', column: categoryColumn(pick((i) => category.get(i))) });
  // a table's custom columns ride along (a theme reads them by name: a score, a consequence), in their types
  for (const [name, c] of Object.entries(rows)) {
    if (KnownColumns.has(name)) continue;
    if (c instanceof CategoryColumn) columns.push({ name, column: categoryColumn(pick((i) => c.get(i))) });
    else if (Array.isArray(c))
      columns.push({
        name,
        column: categoryColumn(pick((i) => (c[i] === null || c[i] === undefined ? null : String(c[i])))),
      });
    else if (ArrayBuffer.isView(c)) {
      const typed = c as unknown as Float64Array;
      columns.push({ name, column: new (typed.constructor as Float64ArrayConstructor)(pick((i) => typed[i])) });
    }
  }
  const table: RangesList = {
    form: 'list',
    start: Int32Array.from(s),
    end: Int32Array.from(e),
    ...(strand ? { strand: Int8Array.from(pick((i) => Strands[strand.get(i) ?? ''] ?? 0)) } : {}),
    ...(id ? { id: categoryColumn(pick((i) => id.get(i))) } : { row: Int32Array.from(rowOf) }),
    ...(columns.length ? { fields: columns } : {}),
  };
  return { items, table, fields: fieldNames, wrappedUnknown };
}

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
  private readonly edgeMaps = new Map<string, MapEntry[]>();
  private readonly codons: CodonPaths[] = [];
  private readonly querySystems = new Map<string, QuerySystem>();
  private readonly queryAnnotations = new Map<string, QueryAnnotation | undefined>();

  constructor(
    readonly spec: Node,
    readonly config: ResolveConfig,
    readonly systems: ReadonlyMap<string, SystemEntry>,
    readonly maps: ReadonlyMap<string, MapEntry>,
    readonly sources: ReadonlyMap<string, SourceBinding>,
    readonly annotations: ReadonlyMap<string, AnnotationEntry>,
    readonly issues: readonly Issue[],
    private readonly carried?: ReadonlyMap<string, PathEntry>,
    readonly lengths: Lengths = new Map(),
  ) {
    this.hash = `${hashNode(spec)}${lengthsKey(lengths)}`;
    for (const m of maps.values()) {
      const fwd: Edge = { map: m.address, from: m.from, to: m.to, backward: false, step: m.step, linear: m.linear };
      const back: Edge = {
        map: m.address,
        from: m.to,
        to: m.from,
        backward: true,
        step: -m.step as Step,
        linear: m.linear.inverse(),
      };
      for (const e of [fwd, back]) {
        let list = this.edges.get(e.from);
        if (!list) {
          this.edges.set(e.from, (list = []));
          this.edgeMaps.set(e.from, []);
        }
        list.push(e);
        this.edgeMaps.get(e.from)!.push(m);
      }
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

  private entry(from: string, to: string): PathEntry {
    const k = `${from}\u0000${to}`;
    let e = this.pathEntries.get(k);
    if (e) return e;
    const old = this.carried?.get(k);
    if (old && this.valid(old)) e = old;
    else {
      const deps = new Map<string, readonly MapEntry[]>();
      const r = findPath(
        {
          edgesFrom: (s) => {
            if (!deps.has(s)) deps.set(s, this.edgeMaps.get(s) ?? []);
            return this.edges.get(s) ?? [];
          },
        },
        from,
        to,
        { maxSteps: this.config.maxPathSteps, maxCandidates: this.config.maxPathCandidates },
      );
      e = { found: r.ok ? { path: r.path } : { why: r.why }, deps };
    }
    this.pathEntries.set(k, e);
    return e;
  }

  path(from: string, to: string): PathResult {
    const k = `${from}\u0000${to}`;
    let r = this.results.get(k);
    if (r) return r;
    // a path read the other way is the inverse of the one found
    const back = this.pathEntries.get(`${to}\u0000${from}`);
    const e = back && !this.pathEntries.has(k) ? back : this.entry(from, to);
    if ('path' in e.found) r = { path: e === back ? inversePath(e.found.path) : e.found.path };
    else r = { issue: pathIssue(e.found.why, from, to) };
    this.results.set(k, r);
    return r;
  }

  project(loc: Location, to: string, options: { close?: boolean } = {}): Projection {
    const src = options.close ? this.close(loc) : loc;
    const issues: Issue[] = [];
    const parts: { system: string; ranges: number[] }[] = [];
    for (const part of src.parts) {
      const r = this.path(part.system, to);
      if (!r.path) {
        if (r.issue) issues.push(r.issue);
        continue;
      }
      parts.push({ system: to, ranges: projectRanges(r.path, part.ranges) });
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

  // ---- the query forms' model ----

  system(address: string): QuerySystem | undefined {
    let s = this.querySystems.get(address);
    if (s) return s;
    const e = this.systems.get(address);
    if (!e) return undefined;
    s = {
      address,
      length: e.length ?? this.config.openQueryLength,
      cyclic: e.cyclic,
      numbering: () => undefined,
      numberingIds: () => [],
    };
    this.querySystems.set(address, s);
    return s;
  }

  annotation(address: string): QueryAnnotation | undefined {
    if (this.queryAnnotations.has(address)) return this.queryAnnotations.get(address);
    const e = this.annotations.get(address);
    const a = e?.shape === 'ranges' ? queryAnnotation(e, this.systems.get(e.target), this.config) : undefined;
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

  /** Structures are not built yet: no chain, and the form reports the node missing. */
  chain(): Node | undefined {
    return undefined;
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
      if (!e.items) return { items: [], pending: over ?? rangesLocation(e.target, [0, length]) };
      const part = over && partOf(over, e.target);
      if (!over) return { items: e.items };
      if (!part) return { items: [] };
      return { items: e.items.filter((it) => overlapsRanges(part.ranges, it.start, it.end, length)) };
    },
  };
}
