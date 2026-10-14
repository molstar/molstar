/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The load stages' contract types (seqspec:sources#the-loader): the host's loaders, format requests and results, the
// stage context. Additions to the spec's shapes are marked: `FormatRequest.want` and `.geometry`, and `AlignerInput`.

import type { AlignmentInput, Issue, Molecule, MoleculeInput, ParsedTree, SeqSet } from "@molstar/seq-core/model";
import type { Geometry } from "./doc.js";
import type { SvsRegistries } from "./registry.js";

/** One file the host loads: a URI or a host file; `maxBytes` its bound. */
export interface LoadRequest {
  uri: string | null;
  hostFile: string | null;
  maxBytes: number | null;
}

/**
 * An alignment region's samples (seqspec:sources#levels-and-the-zoom-gate): `TILE` × `TILE`, row-major, by class and
 * share.
 */
export interface Tile {
  readonly cls: Uint8Array;
  readonly share: Uint8Array;
}

/** The host's port for all I/O. Tests pass fixture loaders: a URI → file map, a scripted stream provider. */
export interface Loaders {
  /** Set by `archiveLoaders`: the document is read from an archive, its main file in this directory (`ValidateOptions.base`). */
  readonly base?: string;
  load(request: LoadRequest, signal: AbortSignal): Promise<Uint8Array>;
  /** Optional: a parse the host shares (one Mol* model for the 3D view and a complex provider); null: read the bytes. */
  parsed?(request: LoadRequest, format: string, signal: AbortSignal): Promise<FormatResult | null>;
  /**
   * A streamed source's provider (seqspec:sources#streamed-sources), by its `stream` id or its `protocol`; null: none
   * for it (`svs.provider.unavailable`).
   */
  stream?(source: StreamSource): StreamProvider | null;
}

// The stream contract (seqspec:sources#dynamic-nodes): wire shapes, in the format's names.
export type StreamDataKind = "letters" | "items" | "values" | "alignment";
export interface StreamSource {
  /** The source node's ref, or a nested source's derived ref. */
  readonly ref: string;
  readonly stream: string | null;
  readonly protocol: string | null;
  readonly uri: string | null;
  readonly data_kind: StreamDataKind;
  /** The source that added it; null: the document's. */
  readonly parent: string | null;
  /** Its readers' shape params (`row_count`, `column_count`, `length`, `target`; a statistic's `measure`), for a provider with no manifest of its own. */
  readonly declared: Readonly<Record<string, unknown>>;
}
/** The provider does not serve this resource: reported once, never retried. */
export class StreamNotServed extends Error {
  constructor(message: string) {
    super(message);
    this.name = "StreamNotServed";
  }
}
export interface StreamRequest {
  source: string;
  /** A manifest key; null: the manifest. */
  system: string | null;
  start: number;
  stop: number;
  /** An alignment's rows, else 0. */
  row_start: number;
  row_stop: number;
  /** Into the manifest's `levels`; -1: the data itself. */
  level: number;
  labels: boolean;
  priority: "view" | "margin" | "prefetch";
  /**
   * A data request's bound on what a provider reads for it (ST4b, the zoom gate by size): past it, as past its own
   * item bound, it answers with a level's bins, estimated where it can without reading the data. Absent or null: none.
   */
  max_bytes?: number | null;
  /** A data request the reader asked for whole (a "Load anyway"): the provider reads it, no zoom gate (ST4b). */
  force?: boolean;
}
export interface Manifest {
  systems: { key: string; extent: [number, number]; numbering: Record<string, unknown> | null; label: string | null }[];
  /** Finest first: positions per bin, or an alignment's [rows, columns] per bin. */
  levels: (number | [number, number])[];
  domain?: [number, number];
  vocabulary?: Record<string, string[]>;
  rows?: number;
  columns?: number;
  alphabet?: "protein" | "dna" | "rna" | "generic" | null;
  inserts?: [number, number][];
  labels?: boolean;
  nodes?: StreamNode[];
}
export interface StreamNode {
  kind: "sequence" | "ligand" | "unit" | "mapping" | "annotation" | "source";
  key: string;
  params: Record<string, unknown>;
}
/** A table in the columns form (seqspec:format#positions): typed columns allowed (`share` as a `Uint8Array`). */
export type StreamTable = Readonly<Record<string, ArrayLike<unknown>>>;
export type DataPatch = {
  source: string;
  system: string | null;
  start: number;
  stop: number;
  row_start: number;
  row_stop: number;
  level: number;
  nodes?: StreamNode[];
} & (
  | { kind: "manifest"; manifest: Manifest }
  | { kind: "items"; items: StreamTable }
  | { kind: "bins"; bins: StreamTable; estimate?: boolean } // `estimate`: from an index, not counted (ST4b)
  | { kind: "letters"; letters: string }
  | { kind: "rows"; labels: string[]; rows: string[] | null }
);
export interface StreamProvider {
  /** Every request, the manifest first. Rejects with `StreamNotServed`, any other error (a failure), or the signal's reason. */
  fetch(request: StreamRequest, signal: AbortSignal): Promise<DataPatch>;
}
/** A streamed source compiled: its source, manifest and provider (`SourceResult`, and the readers' `svs:stream`). */
export interface Streamed {
  readonly source: StreamSource;
  readonly manifest: Manifest;
  readonly provider: StreamProvider;
}

export interface FormatRequest {
  record: string | number | null;
  tag: string | null;
  categoryName: string | null;
  blockHeader: string | null;
  blockIndex: number;
  flavor: "aligned" | "a3m";
  /** Added: what the reader wants (FASTA is molecules for a sequence, an alignment for an alignment); null: the format's first. */
  want: FormatResult["kind"] | null;
  /** Added: the reading annotation's geometry (a Stockholm `SS_cons` read as pairs gives base pairs). */
  geometry: Geometry | null;
}

export type FormatResult =
  | { kind: "molecules"; molecules: MoleculeInput[] }
  | { kind: "alignments"; alignments: AlignmentInput[] }
  | { kind: "tree"; tree: ParsedTree }
  /**
   * `text`: the cells are strings (TSV, CSV), typed after `fields` by their column's role
   * (seqspec:format#annotation-tables).
   */
  | { kind: "table"; columns: Record<string, ArrayLike<unknown>>; length: number; text?: boolean }
  /** An extension format's parsed source (a structure), for its complex and annotation providers. */
  | { kind: "object"; value: unknown };

export interface StageContext {
  readonly signal: AbortSignal;
  readonly registries: SvsRegistries;
  /** The set so far; providers get the set of every data node. */
  readonly set: SeqSet | null;
  /** A source loaded for this document, by the source's ref or the ref of a data node read from it: its default parse. */
  source(ref: string): FormatResult | null;
  /** Added: the source nodes' results by ref, filled by the caller from `loadNode` on each (`sourceParse` answers `source`). */
  readonly sources: ReadonlyMap<string, SourceResult>;
  /** Added: refs of unavailable nodes (a stage's result said so), filled by the caller: their readers are left out quietly. */
  readonly unavailable: ReadonlySet<string>;
  report(issue: Issue): void;
}

/** Added: a `source` node's load result (`LoadResult.source`), what its readers take through `StageContext.sources`. */
export type SourceResult =
  /** The default request's parse (null: unreadable, reported) and a parse per request, each once, kept here. */
  | {
      readonly kind: "loaded";
      /** The registered format that reads it. */
      readonly format: string;
      readonly parsed: FormatResult | null;
      read(request: FormatRequest, ctx: StageContext): Promise<FormatResult | null>;
    }
  /** A streamed source: its manifest fetched, its data by region through its provider. */
  | ({ readonly kind: "streamed" } & Streamed)
  /** The bytes could not be loaded (reported): its readers get nothing. */
  | { readonly kind: "failed" }
  /** No package reads it: a format, service or protocol not registered or served (`svs.provider.unavailable`). */
  | { readonly kind: "unavailable" };

/** What an aligner reads of a molecule; an alignment row's own residues are aligned too (as a polymer). */
export type AlignerInput = Pick<Molecule, "key" | "kind" | "alphabet" | "length" | "code">;

export const emptyRequest = (over: Partial<FormatRequest> = {}): FormatRequest => ({
  record: null,
  tag: null,
  categoryName: null,
  blockHeader: null,
  blockIndex: 0,
  flavor: "aligned",
  want: null,
  geometry: null,
  ...over,
});
