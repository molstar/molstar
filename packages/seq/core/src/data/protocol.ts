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
 * Cells); `pairs`, `matrix` and `grid` by region come later.
 */
export type Shape = 'letters' | 'ranges' | 'values' | 'cells';

export const Shapes: readonly Shape[] = ['letters', 'ranges', 'values', 'cells'];

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
}

export interface ShapeInfo {
  readonly system: string;
  readonly shape: Shape;
  readonly fields?: readonly FieldInfo[];
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
}

/** `letters`: one letter per position from `start`, a string, or codes into an alphabet. */
export interface LettersTable {
  readonly form: 'letters';
  readonly start: number;
  readonly letters: string | Uint8Array;
  /** The alphabet the codes index (required with codes). */
  readonly alphabet?: string;
}

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

export interface Tables {
  readonly letters?: LettersTable;
  readonly ranges?: RangesTable;
  readonly values?: ValuesTable;
  readonly cells?: CellsTable;
}

export type Table = LettersTable | RangesTable | ValuesTable | CellsTable;

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
  dispose(): void;
}

/** A format runtime, a registry entry: it opens the sources whose `format` it reads. */
export interface FormatRuntime {
  /** `'seq'`, `'fasta'`, `'bigwig'`. */
  readonly name: string;
  readonly formats: readonly string[];
  open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource>;
}

/** The rejection of an aborted read: never an error state. */
export function abortError(): Error {
  const e = new Error('The read was aborted');
  e.name = 'AbortError';
  return e;
}

export function isAbortError(e: unknown): boolean {
  return typeof e === 'object' && e !== null && (e as { name?: unknown }).name === 'AbortError';
}
