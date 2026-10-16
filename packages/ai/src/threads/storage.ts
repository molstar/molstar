/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Stops a listening started with `watch` or `subscribe`. */
export type Unsubscribe = () => void;

/** A row's fields, as a backend stores them. */
export type Field = string | number | undefined;
export type Row = Readonly<Record<string, Field>>;

/** One row written or deleted, as a backend stores it. */
export interface RowWrite<T extends string = string> {
  readonly table: T;
  readonly key: string;
  /** The row, or undefined for a delete. */
  readonly row: Row | undefined;
}

/**
 * Where a scope's rows are kept between runs: the threads' storage port (threads.md § Storage). A scope is
 * one thread. The host keeps a loaded scope in memory and writes every command through as one write set, applied
 * only if the scope's generation is still the one it read.
 */
export interface RowStore<T extends string = string> {
  /** A scope's rows and generation, or undefined when it has none. */
  load(scope: string): Promise<{ generation: number; rows: readonly RowWrite<T>[] } | undefined>;
  /**
   * Applies the writes atomically and moves the generation from `expected` to `expected + 1`, or throws
   * `RowConflict` when the stored generation is not `expected` (someone else wrote). A missing scope has
   * generation 0.
   */
  apply(scope: string, expected: number, writes: readonly RowWrite<T>[]): Promise<void>;
  /** Forgets a scope: its rows and its generation. */
  drop(scope: string): Promise<void>;
  /**
   * The scope's generation as stored now (0 when it has none): cheap, so a host can check its cached rows are
   * still current before using them, and see another runtime's writes.
   */
  generation(scope: string): Promise<number>;
  /**
   * Hears which scope was written (or dropped) through another handle to the same storage, so another tab's or
   * process's write is seen. A store that cannot tell has none. Writes through this very handle are not told:
   * the host that made them knows.
   */
  watch?(listener: (scope: string) => void): Unsubscribe;
}

export class RowConflict extends Error {
  constructor(scope: string) {
    super(`Rows of ${scope} were written elsewhere.`);
    this.name = 'RowConflict';
  }
}

/** The tables of a thread's scope; every row is its value as JSON (`json`). */
export type ThreadTable = 'thread' | 'entry' | 'stream' | 'deltas' | 'receipt';

export type ThreadRows = RowStore<ThreadTable>;
