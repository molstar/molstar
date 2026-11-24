/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn, Column, ColumnType, Issue, Node, Value } from '@molstar/core/repr/tree';
import type { AssetRangeReader } from '@molstar/core/util/assets';

// The source protocol (_spec/README.md "Data", Protocol): what a format runtime opens, what an opened source covers,
// and the request and answer of one region. Types only, but for the abort error.

/**
 * The table shapes: `cells` an alignment's rows on 1D tiles of its rows system, columns whole (_spec/README.md "Data",
 * Cells); `pairs` on 1D tiles by either end, `grid` on 2D tiles, `matrix` positions tiled and states whole, `rows`
 * keyed records (_spec/README.md "Data", Protocol: the region shapes).
 */
export type Shape = 'letters' | 'ranges' | 'values' | 'cells' | 'pairs' | 'grid' | 'matrix' | 'rows';

export const Shapes: readonly Shape[] = ['letters', 'ranges', 'values', 'cells', 'pairs', 'grid', 'matrix', 'rows'];

/** The shapes whose `list` (or `items`) form a request may ask at any level, as items tiles (`RegionRequest.items`). */
export const ItemShapes: readonly Shape[] = ['ranges', 'pairs', 'matrix'];

/** A field a source's ranges carry besides `start`, `end`, `strand` and `id`. */
export interface FieldInfo {
  readonly name: string;
  readonly type: ColumnType;
  readonly doc?: string;
}

export interface SystemInfo {
  /** The source's own name (a chromosome, a transcript). */
  readonly name: string;
  readonly length?: number;
  readonly alphabet?: string;
  readonly cyclic?: boolean;
  /** Names a spec may use for it (`chr1`, `1`). */
  readonly aliases?: readonly string[];
  /**
   * A sided pairs source's side (between two genomes): `a` the `from` side of the map the pairs serve, `b` the other; a
   * name may be listed once a side (_spec/README.md "Data", Protocol: sided pairs).
   */
  readonly side?: 'a' | 'b';
}

/** A grid's or pairs' whole-assembly system: the listed sequences end to end, in order, served at coarse levels only. */
export interface WholeSystem {
  readonly system: string;
  readonly sequences: readonly string[];
}

export interface ShapeInfo {
  readonly system: string;
  readonly shape: Shape;
  /** A ranges table's fields; a grid's or matrix's channels; a pairs field that can be summed. */
  readonly fields?: readonly FieldInfo[];
  /** A grid's second systems: a list, or `'*'` (every system of the coverage). */
  readonly y?: readonly string[] | '*';
  /** A grid served once (the tile at or above the diagonal), the reader transposing the other (Hi-C). */
  readonly symmetric?: boolean;
  /** How a grid's coarse levels were made: `sum` (counts) or `mean` (scores; when absent). */
  readonly aggregate?: 'sum' | 'mean';
  /** A category the pairs bins keep (`type`). */
  readonly binBy?: string;
  /** A matrix's state names, in order. */
  readonly states?: readonly string[];
  /** A grid's or pairs' whole-assembly system. */
  readonly whole?: WholeSystem;
}

/**
 * A data node a source brings (a transcript, a nested source): identity is (source, key); `node` is a spec node of
 * its kind, checked by the language before it is held. A type only for now.
 */
export interface DynamicNode {
  readonly key: string;
  readonly node: Node;
}

/** What an opened source serves, known after open. */
export interface Coverage {
  readonly systems: readonly SystemInfo[];
  readonly shapes: readonly ShapeInfo[];
  /** Bin widths in positions, ascending; `1` is per position. A source without levels serves `[1]`. */
  readonly levels: readonly number[];
  /** Nodes known up front. */
  readonly nodes?: readonly DynamicNode[];
  /**
   * The fields `findNames` searches (GFF3's `Name`, `ID`, `Alias`), when the source serves it (_spec/README.md "Data",
   * Names).
   */
  readonly names?: readonly string[];
}

/** An item a source found by name (EC-228): where it is, its id when it has one, what to show, the field matched. */
export interface NameHit {
  readonly system: string;
  /** Positions of the system, 0-based, end exclusive: the item's whole extent. */
  readonly start: number;
  readonly end: number;
  /** The item's id, which `(item :id …)` addresses on its source. */
  readonly id?: string;
  readonly label: string;
  /** The field that matched (one of `Coverage.names`). */
  readonly field: string;
}

export interface FindNamesOptions {
  /** At most this many hits (default 20). */
  readonly max?: number;
  /** Names starting with the text too, after the exact ones (default false). */
  readonly prefix?: boolean;
}

export interface RegionRequest {
  readonly system: string;
  /** Positions of the system, 0-based, end exclusive; at a coarser level, aligned to its bin width. */
  readonly start: number;
  readonly end: number;
  /** A bin width the coverage lists. */
  readonly level: number;
  readonly shapes: readonly Shape[];
  /**
   * The values channel asked (a `values` node's `field`: `plddt`, `b_factor`); absent, the source's default value. A
   * source serving one value ignores it (_spec/README.md "Data", Protocol).
   */
  readonly field?: string;
  /** A grid's second axis, at the same level (square levels): positions of `system`, 0-based, end exclusive. */
  readonly y?: { readonly system: string; readonly start: number; readonly end: number };
  /** The side a request on a sided pairs source reads (`SystemInfo.side`). */
  readonly side?: 'a' | 'b';
  /**
   * Items at a tile size, not bins: the `list` form of `ranges` and `pairs` (a matrix's `items`) at any listed level,
   * every item with an end in the region, uncut (_spec/README.md "Data", Protocol: items tiles).
   */
  readonly items?: boolean;
}

/** `letters`: one letter per position from `start`, a string, or codes into an alphabet. */
export interface LettersTable {
  readonly form: 'letters';
  readonly start: number;
  readonly letters: string | Uint8Array;
  /** The alphabet the codes index (required with codes). */
  readonly alphabet?: string;
}

/**
 * A VCF source's fixed fields, in order (_spec/README.md "Data", VCF, Fields): the label and `name` the VCF `ID`, the
 * alleles and the record's; every other field of its rows is an INFO key. Here, not in the reader, so a describer reads
 * them without the reader.
 */
export const VcfFixedFieldNames = ['label', 'name', 'ref', 'alt', 'qual', 'filter', 'category'] as const;

/** A ranges field: a named column, one value a row (any name: `kind` and `$x` too). */
export interface RangesField {
  readonly name: string;
  readonly column: Column;
}

/**
 * `ranges` as a list: `strand` −1, 0, 1 (absent: 0); `fields` the columns the coverage lists. A row's identity is its
 * `id`, else `row` (stable across the tiles it is cut into): readers dedupe a range over a tile's edge by it and its
 * start; rows of one id with other starts are pieces of one item.
 */
export interface RangesList {
  readonly form: 'list';
  readonly start: Int32Array;
  readonly end: Int32Array;
  readonly strand?: Int8Array;
  readonly id?: CategoryColumn;
  readonly row?: Int32Array;
  readonly fields?: readonly RangesField[];
  /**
   * Each row's item's whole extent, when rows are pieces of one item (BED12's blocks, a transcript's exons): a reader of
   * region data knows an item's ends from any tile (_spec/README.md "Data", Protocol).
   */
  readonly itemStart?: Int32Array;
  readonly itemEnd?: Int32Array;
}

/**
 * Pairs as a list (_spec/README.md "Data", Protocol): both ends of each row, every pair with an end in the region asked,
 * the other end anywhere. A pair's identity is its `row`, else its `id` (never `.` or empty), with both ends' systems
 * and starts; rows sharing `row` with other starts are the pieces of one alignment.
 */
export interface PairsList {
  readonly form: 'list';
  /** Each end's system (absent: both ends on the requested system). */
  readonly aSystem?: CategoryColumn;
  readonly bSystem?: CategoryColumn;
  readonly aStart: Int32Array;
  readonly aEnd: Int32Array;
  readonly bStart: Int32Array;
  readonly bEnd: Int32Array;
  /** −1, 0, 1 (EC-206). */
  readonly aStrand?: Int8Array;
  readonly bStrand?: Int8Array;
  readonly type?: CategoryColumn;
  readonly id?: CategoryColumn;
  readonly row?: Int32Array;
  readonly fields?: readonly RangesField[];
}

/**
 * `pairs` at a coarser level (the pairs level): each bin pair with at least one pair, its count, by relative strand (1
 * straight, −1 reversed, 0 none) and by `category` when the coverage's `binBy` lists it, a field summed when asked
 * (`sum`). Bins are of the level's width on each end's system, aligned on position 0.
 */
export interface PairsBins {
  readonly form: 'bins';
  readonly aSystem?: CategoryColumn;
  readonly bSystem?: CategoryColumn;
  readonly aBin: Int32Array;
  readonly bBin: Int32Array;
  readonly strand?: Int8Array;
  readonly category?: CategoryColumn;
  readonly count: Float32Array;
  readonly sum?: Float32Array;
}

/** `grid` dense: `nx × ny` cells from bins (`x`, `y`) at the level, x-major, NaN missing. */
export interface GridDense {
  readonly form: 'dense';
  readonly x: number;
  readonly y: number;
  readonly nx: number;
  readonly ny: number;
  readonly value: Float32Array;
}

/** `grid` sparse: the listed cells by bin (absolute, at the level); others are `fill` (0 counts, NaN missing). */
export interface GridSparse {
  readonly form: 'sparse';
  readonly x: number;
  readonly y: number;
  readonly xBin: Int32Array;
  readonly yBin: Int32Array;
  readonly value: Float32Array;
  readonly fill: number;
}

/**
 * `matrix` at its finest level: a column an item (a variant, a CpG) at its own extent × every state, item-major,
 * `value` float32 (NaN missing) or `codes` uint8 (`0` no copy, `1` some, `2` every copy, `3` partly called, `255`
 * uncalled).
 */
export interface MatrixItems {
  readonly form: 'items';
  readonly start: Int32Array;
  readonly end: Int32Array;
  readonly states: number;
  readonly value?: Float32Array;
  readonly codes?: Uint8Array;
  /**
   * With `codes`, the call a cell was written as (a VCF's `0|1`, `1/2`, `1`, `./.`): an index a cell into `callTexts`
   * (0 none). The class says how many copies of this allele; the call says what the sample has.
   */
  readonly calls?: Uint8Array;
  readonly callTexts?: readonly string[];
  readonly id?: CategoryColumn;
  readonly row?: Int32Array;
}

/**
 * `matrix` at a coarser level: a bin per column from `start` (aligned), every state, bin-major. Its own form, `summary`
 * (not `bins`): a pairs table's `bins` and a matrix's are told apart by name, not by their columns.
 */
export interface MatrixSummary {
  readonly form: 'summary';
  readonly start: number;
  readonly states: number;
  readonly mean?: Float32Array;
  readonly min?: Float32Array;
  readonly max?: Float32Array;
  readonly count?: Float32Array;
  /** Codes: the majority class (255 none) and its share (0 to 255), as cells' coarse form. */
  readonly classes?: Uint8Array;
  readonly shares?: Uint8Array;
}

/** `rows`: keyed records on the source's records system, a position a record from `start`. */
export interface RowsKeyed {
  readonly form: 'keyed';
  readonly start: number;
  readonly key: CategoryColumn;
  readonly fields: readonly RangesField[];
}

/** `ranges` at a coarse level as density: a count per bin (bin index on the system, at the level's width). */
export interface RangesDensity {
  readonly form: 'density';
  readonly bin: Int32Array;
  readonly count: Float32Array;
  readonly category?: CategoryColumn;
}

/** `values` at level 1: one value per position from `start` (NaN: missing). */
export interface ValuesPositions {
  readonly form: 'positions';
  readonly start: number;
  readonly value: Float32Array;
}

/** `values` at a coarser level: one summary per bin of the level's width from `start` (aligned to it). */
export interface ValuesBins {
  readonly form: 'bins';
  readonly start: number;
  readonly min: Float32Array;
  readonly max: Float32Array;
  readonly mean: Float32Array;
  readonly count: Float32Array;
}

/**
 * `cells` at level 1: a tile's rows from `start` × every column (`columns`), row-major letter bytes (ASCII: gaps `-`
 * and `.`), with the rows table: each row's `label` (its key), its first residue's number (`first`, EC-39) and its own
 * fields (`fields`: Stockholm's `#=GS` tags, an inline table's other columns), a column each over the tile's rows.
 */
export interface CellsCodes {
  readonly form: 'codes';
  readonly start: number;
  readonly columns: number;
  readonly codes: Uint8Array;
  readonly label?: CategoryColumn;
  readonly first?: Int32Array;
  readonly fields?: readonly RangesField[];
}

/**
 * `cells` at a coarser level: a row bin (from `start`, aligned to the level) per column, its majority letter byte
 * (`classes`, 255 none) and that letter's share (`shares`, 0 to 255), row-major.
 */
export interface CellsMajority {
  readonly form: 'majority';
  readonly start: number;
  readonly columns: number;
  readonly classes: Uint8Array;
  readonly shares: Uint8Array;
}

export type RangesTable = RangesList | RangesDensity;
export type CellsTable = CellsCodes | CellsMajority;
export type ValuesTable = ValuesPositions | ValuesBins;
export type PairsTable = PairsList | PairsBins;
export type GridTable = GridDense | GridSparse;
export type MatrixTable = MatrixItems | MatrixSummary;
export type RowsTable = RowsKeyed;

export interface Tables {
  readonly letters?: LettersTable;
  readonly ranges?: RangesTable;
  readonly values?: ValuesTable;
  readonly cells?: CellsTable;
  readonly pairs?: PairsTable;
  readonly grid?: GridTable;
  readonly matrix?: MatrixTable;
  readonly rows?: RowsTable;
}

export type Table =
  | LettersTable
  | RangesTable
  | ValuesTable
  | CellsTable
  | PairsTable
  | GridTable
  | MatrixTable
  | RowsTable;

/**
 * Why a region was not answered: `failed` is retried with backoff then a cooldown, `malformed` fails at once with a
 * cooldown, `not-served` and `too-large` (ask a coarser level) are never retried.
 */
export type FailureKind = 'failed' | 'not-served' | 'malformed' | 'too-large';

export const FailureKinds: readonly FailureKind[] = ['failed', 'not-served', 'malformed', 'too-large'];

export type RegionAnswer =
  | { readonly ok: true; readonly tables: Tables; readonly nodes?: readonly DynamicNode[] }
  | { readonly ok: false; readonly kind: FailureKind; readonly issue: Issue };

/** A source as a spec names it: its address, its format, and its data inline (`value`) or by location (`uri`). */
export interface SourceEntry {
  readonly address: string;
  readonly format: string;
  readonly uri?: string;
  /** A host's name for the data (an edited model with no file): read by a host's adapter only, never by a default one. */
  readonly ref?: string;
  /** Where the source's index is, for a format read with one (FASTA's `.fai`). */
  readonly index?: string;
  readonly value?: Value;
  // the transform: what a file of records' columns mean (_spec/README.md "Data", Records)
  /** A role or a field of the ranges by column name or dotted path. */
  readonly fields?: Readonly<Record<string, string>>;
  /** The dotted path to the records. */
  readonly records?: string;
  /** What the source serves: its `records` (absent), or a `matrix` (an array of arrays a grid reads whole). */
  readonly shape?: 'records' | 'matrix';
  /** The dotted path to the letters of the source's own system. */
  readonly letters?: string;
  /** The source's own system's name. */
  readonly system?: string;
  /** The number of the first position. */
  readonly base?: 0 | 1;
  readonly ends?: 'inclusive' | 'exclusive';
  /** A field's declared type. */
  readonly types?: Readonly<Record<string, 'int32' | 'float32' | 'float64' | 'category' | 'string'>>;
}

/** The asset port as a runtime sees it: text or bytes by location, never the network in tests. */
export interface AssetReader {
  readText(uri: string, signal: AbortSignal): Promise<string>;
  /** The bytes as stored (a MessagePack file; a text read whole whose size is bounded before it is parsed). */
  readBytes(uri: string, signal: AbortSignal): Promise<Uint8Array>;
  /**
   * Byte ranges of the asset (core's `AssetRangeReader`): an indexed format (bigWig) reads a header, an index and a few
   * blocks. Optional: a format reads the file whole through `readBytes` when the host has none.
   */
  openRange?(uri: string, signal: AbortSignal): Promise<AssetRangeReader>;
}

/** A parse's bounds, checked before it allocates (_spec/README.md "Data", Records). */
export interface OpenLimits {
  /** Rows a source of records may hold. */
  readonly maxRows: number;
  /** Bytes a file read whole and parsed may have. */
  readonly maxBytes: number;
}

/** The store's default settings' (`maxRows`, `maxBytes`): a million rows, 64 MB. */
export const DefaultOpenLimits: OpenLimits = { maxRows: 1_000_000, maxBytes: 64 * 1024 * 1024 };

export interface OpenContext {
  readonly assets: AssetReader;
  readonly signal: AbortSignal;
  /** The store's settings, passed by whoever opens the source. */
  readonly limits: OpenLimits;
}

export interface OpenedSource {
  readonly coverage: Coverage;
  /**
   * What opening said without refusing the source (the rows a file of records left out or clipped), each with its
   * address: whoever opens it reports them with the source's own issues.
   */
  readonly issues?: readonly Issue[];
  /** Rejects with an `AbortError` when `signal` aborts; any other rejection is a `failed` answer. */
  read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer>;
  /**
   * Optional: tells `fn` when the data changed (a host edited its structure), until the returned unsubscribe is called;
   * a change since the source was opened and before `fn` was given is told too, once, soon after. Whoever attached the source opens it again and reattaches it in place (the viewer's source sync, coalesced a frame:
   * _spec/README.md "Data", Protocol). A source read once (a file) has none.
   */
  changes?(fn: () => void): () => void;
  /**
   * Optional: items whose names (the fields `coverage.names` lists) match `text`, case-insensitive, the exact ones first,
   * then by prefix when asked, at most `max`; rejects with an `AbortError` when `signal` aborts (_spec/README.md "Data",
   * Names).
   */
  findNames?(text: string, options: FindNamesOptions, signal: AbortSignal): Promise<readonly NameHit[]>;
  /**
   * Optional: what reading `request` costs before it is read (compressed bytes; a tabix source's `regionBytes`), known
   * exactly from an index, never estimated; `undefined` when not known. The store sums it over a layer's tiles
   * (`maxBytesPerLayer`).
   */
  readCost?(request: RegionRequest): number | undefined;
  dispose(): void;
}

/** A format runtime, a registry entry: it opens the sources whose `format` it reads. */
export interface FormatRuntime {
  /** `'seq'`, `'fasta'`, `'bigwig'`. */
  readonly name: string;
  readonly formats: readonly string[];
  open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource>;
}

/**
 * The name field under which a `seq` source finds its dynamic nodes by key, when its coverage lists it in `names`
 * (_spec/README.md "Data", The seq format).
 */
export const SeqNodeKeyField = 'key';

/** The rejection of an aborted read: never an error state. */
export function abortError(): Error {
  const e = new Error('The read was aborted');
  e.name = 'AbortError';
  return e;
}

export function isAbortError(e: unknown): boolean {
  return typeof e === 'object' && e !== null && (e as { name?: unknown }).name === 'AbortError';
}
