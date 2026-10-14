/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Numbers, strings, issues, provenance and custom data (seqspec:model#custom-data). `-1` means none for every brand.

import type { SeqSet } from "./set.js";

export type Brand<T, B extends string> = T & { readonly "@type": B };
/** A position in a sequence, 0..length-1. */
export type SeqIndex = Brand<number, "seq-index">;
/** Between positions, 0..length. */
export type BoundaryIndex = Brand<number, "boundary">;
export type ColumnIndex = Brand<number, "column">;
/** An alignment row, 0..rowCount-1. */
export type RowIndex = Brand<number, "row">;
/** A sequence or a ligand. */
export type MoleculeId = Brand<number, "molecule">;
export type UnitId = Brand<number, "unit">;
export type ComplexId = Brand<number, "complex">;
export type GroupId = Brand<number, "copyGroup">;
export type AlignmentId = Brand<number, "alignment">;
export type AnnotationId = Brand<number, "annotation">;
export type MappingId = Brand<number, "mapping">;
export type ItemIndex = Brand<number, "item">;
/** An interned string. */
export type NameId = Brand<number, "name">;

/**
 * Each distinct string once; append-only, shared by every version of one lineage (seqspec:model#identity-and-lineage).
 */
export class StringTable {
  private readonly ids = new Map<string, number>();
  private readonly values: string[] = [];
  get size() {
    return this.values.length;
  }
  intern(s: string): NameId {
    let id = this.ids.get(s);
    if (id === undefined) {
      id = this.values.length;
      this.values.push(s);
      this.ids.set(s, id);
    }
    return id as NameId;
  }
  /** `""` for -1. */
  get(id: number): string {
    return id < 0 ? "" : (this.values[id] ?? "");
  }
  /** -1 when absent. */
  find(s: string): NameId {
    return (this.ids.get(s) ?? -1) as NameId;
  }
}

/**
 * An issue's severity (seqspec:format#validation): `error` invalid input, `warning` a loss the person should know
 * about, `info` a quiet note on how valid input was read.
 */
export type Severity = "error" | "warning" | "info";

/** The one shape every package reports. `path` is a JSON pointer into the input. */
export interface Issue {
  code: string;
  severity: Severity;
  path: string;
  message: string;
}

export const Issue = {
  error: (code: string, path: string, message: string): Issue => ({ code, severity: "error", path, message }),
  warning: (code: string, path: string, message: string): Issue => ({ code, severity: "warning", path, message }),
  info: (code: string, path: string, message: string): Issue => ({ code, severity: "info", path, message }),
};

/** seqspec:format#positions; `method` an open vocabulary, `confidence` NaN when not stated. */
export interface Provenance {
  readonly source: NameId;
  readonly method: NameId;
  readonly note: NameId;
  readonly confidence: number;
  /** NameId per DOI or URL. */
  readonly references: Int32Array;
  /** A program's JSON (an edge value). */
  readonly parameters: unknown;
}

export interface ProvenanceInput {
  source: string;
  method?: string | null;
  note?: string | null;
  confidence?: number | null;
  references?: string[];
  parameters?: unknown;
}

export const Provenance = {
  create(p: ProvenanceInput, strings: StringTable): Provenance {
    const opt = (s: string | null | undefined) => (s == null ? (-1 as NameId) : strings.intern(s));
    return {
      source: strings.intern(p.source),
      method: opt(p.method),
      note: opt(p.note),
      confidence: p.confidence ?? Number.NaN,
      references: Int32Array.from(p.references ?? [], (r) => strings.intern(r)),
      parameters: p.parameters ?? null,
    };
  },
};

/** Every model object's node `custom` and runtime properties. */
export interface CustomData {
  /** The node's `custom`, by namespace, uninterpreted. */
  readonly doc: Readonly<Record<string, unknown>>;
  /** Runtime properties: the only mutable part. */
  readonly props: Map<string, unknown>;
}

/** A runtime property, namespaced; `compute` runs lazily, once per object (a new version starts without it). */
export interface CustomProperty<O, T> {
  readonly name: string;
  compute?(object: O, set: SeqSet): T;
}

export const CustomData = {
  create: (doc?: object, props?: Record<string, unknown>): CustomData => ({
    doc: (doc ?? {}) as Record<string, unknown>,
    props: new Map(Object.entries(props ?? {})),
  }),
  /** The property's value: attached, else computed once and kept; undefined when neither. */
  get<O, T>(data: CustomData, p: CustomProperty<O, T>, object: O, set: SeqSet): T | undefined {
    if (data.props.has(p.name)) return data.props.get(p.name) as T;
    if (!p.compute) return undefined;
    const value = p.compute(object, set);
    data.props.set(p.name, value);
    return value;
  },
  attach<O, T>(data: CustomData, p: CustomProperty<O, T>, value: T): void {
    data.props.set(p.name, value);
  },
};

/**
 * The double nearest the shortest decimal that reads back as the float32 `x` (0.15, not 0.15000000596046448); -0 as 0.
 * 7 digits, else 8, else 9. Fewer than 7 can win only where float32's step exceeds 7 digits' (a 7-digit neighbor reads
 * back too); another 8 digits only at a power of two (2^-96, 2^87, 2^90). Checked on every float32.
 */
export function float32Number(x: number): number {
  if (x === 0 || !Number.isFinite(x)) return x + 0;
  const y = Number(x.toPrecision(7));
  if (Math.fround(y) === x) {
    const d = 0.999999 * 10 ** (Math.floor(Math.log10(Math.abs(y))) - 6);
    if (Math.fround(y + d) !== x && Math.fround(y - d) !== x) return y;
    for (let p = 1; p < 7; p++) if (Math.fround(Number(x.toPrecision(p))) === x) return Number(x.toPrecision(p));
    return y;
  }
  const z = Number(x.toPrecision(8));
  if (Math.fround(z) === x) return z;
  const w = Number(x.toPrecision(9));
  const u = Math.abs(x) === 2 ** Math.floor(Math.log2(Math.abs(x))) ? Number((2 * w - z).toPrecision(8)) : w;
  return Math.fround(u) === x ? u : w;
}

/**
 * A typed column's numbers as the JSON numbers they mean (seqspec:format#values): a `Float32Array`'s by
 * `float32Number`, so they key, compare and print as the JSON text does; NaN as `nan`. Made once per column where it is
 * read.
 */
export function columnNumbers<T>(a: ArrayLike<number>, nan: T): (number | T)[] {
  const f32 = a instanceof Float32Array;
  return Array.from(a, (x) => (Number.isNaN(x) ? nan : f32 ? float32Number(x) : x));
}
