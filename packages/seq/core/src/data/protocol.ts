/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CategoryColumn, Column, ColumnType, Issue, Node, Value } from '@molstar/core/repr/tree';

// The source protocol (_spec/README.md "Data", Protocol): what a format runtime opens, what an opened source covers,
// and the request and answer of one region. Types only, but for the abort error.

/** The table shapes of the first slice; `pairs`, `matrix`, `grid` and `rows` come later. */
export type Shape = 'letters' | 'ranges' | 'values';

export const Shapes: readonly Shape[] = ['letters', 'ranges', 'values'];

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
 * `id`, else `row` (stable across the tiles it is cut into): readers dedupe a range over a tile's edge by it.
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

export type RangesTable = RangesList | RangesDensity;
export type ValuesTable = ValuesPositions | ValuesBins;

export interface Tables {
  readonly letters?: LettersTable;
  readonly ranges?: RangesTable;
  readonly values?: ValuesTable;
}

export type Table = LettersTable | RangesTable | ValuesTable;

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
  readonly value?: Value;
}

/** The asset port as a runtime sees it: bytes by location, never the network in tests. */
export interface AssetReader {
  readText(uri: string, signal: AbortSignal): Promise<string>;
}

export interface OpenContext {
  readonly assets: AssetReader;
  readonly signal: AbortSignal;
}

export interface OpenedSource {
  readonly coverage: Coverage;
  /** Rejects with an `AbortError` when `signal` aborts; any other rejection is a `failed` answer. */
  read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer>;
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
