/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The genome chunk's door (_spec/README.md "Resolve", Assemblies, Loaded when used): assembly resolve, the names rule,
// derived sequences and the map between assemblies are one chunk, imported the first time a spec holds an assembly
// (and the genome query forms with it). Until it lands an assembly is pending, never wrong. Why: a reader that never
// opens a genome pays nothing for it.

import { loadGenomeForms } from '@molstar/seqx-core/query';
import type { AssemblyNames, NamedSequence } from './assembly-names.js';
import type { Genome } from './genome-chunk.js';
import type { Placed } from './model-molecules.js';

/** An assembly's sequence as derived: its printed name, the name its source knows it by, and its own settings. */
export interface AssemblySequence extends NamedSequence {
  /** The name its source lists it by (the printed name when inline only). */
  readonly sourceName?: string;
  readonly cyclic?: boolean;
}

export interface AssemblyEntry {
  readonly address: string;
  /** Its accession or name, and its other names: what `source.assembly` is checked against. */
  readonly names: readonly string[];
  readonly source?: string;
  readonly aliasSource?: string;
  /** Its source has not opened: no sequence is derived yet (a section of it is pending, EC-112). */
  readonly pending: boolean;
  /** Its sequences, in its `order` before `include`: what `system`'s categorical domain reads. */
  readonly sequences: readonly AssemblySequence[];
  readonly order: readonly number[];
  /** The systems derived (after `include`), in order; a sequence left out has none. */
  readonly systems: readonly string[];
  /** Each sequence's system by its index, or undefined when `include` leaves it out. */
  readonly systemAt: readonly (string | undefined)[];
  readonly lookup: AssemblyNames;
  /** The index of the sequence a name is (its names, aliases, the `chr` prefix), or undefined. */
  find(name: string): number | undefined;
  /** A section's `include` over it (names through `find`, patterns, a length floor); undefined: every block. */
  include(value: unknown): ((i: number, s: AssemblySequence) => boolean) | undefined;
}

let chunk: typeof Genome | undefined;
let loading: Promise<void> | undefined;
let failed: unknown;

/** The genome chunk when it has landed, else undefined (`loadGenome` asks for it). */
export function genomeLoaded(): typeof Genome | undefined {
  return chunk;
}

/** Why the last load of the genome chunk failed, or undefined: resolve says it (`seq.resolve.genome.load`). */
export function genomeLoadFailure(): unknown {
  return failed;
}

/** What lands with the genome chunk (a renderer's genome drawing, `withGenome`): each loaded with it. */
const companions: (() => Promise<unknown>)[] = [];

/**
 * A chunk that lands with the genome chunk (render's genome drawing registers one): loaded with it, so a section of
 * assemblies, pending until the chunk lands, is never drawn before it; loaded at once when the chunk already has.
 */
export function withGenome(load: () => Promise<unknown>): void {
  companions.push(load);
  if (chunk) void load();
}

/**
 * Loads the genome chunk, the genome query forms and its companions (`withGenome`), once while it is in flight. The promise never rejects: a failed
 * import is kept (`genomeLoadFailure`) and the next call tries again, so a model never waits on a rejected load, and no
 * resolve retries by itself (a reader asks again: `loadGenome`).
 */
export function loadGenome(): Promise<void> {
  if (chunk) return Promise.resolve();
  failed = undefined;
  loading ??= Promise.all([import('./genome-chunk.js'), loadGenomeForms(), ...companions.map((f) => f())]).then(
    ([m]) => {
      chunk = m.Genome;
      loading = undefined;
    },
    (e: unknown) => {
      failed = e ?? new Error('the genome chunk did not load');
      loading = undefined;
    },
  );
  return loading;
}

/** An assembly while the genome chunk loads: pending, its sources named (so they open meanwhile), nothing derived. */
export function pendingAssembly(
  at: Placed,
  source: string | undefined,
  aliasSource: string | undefined,
): { entry: AssemblyEntry; systems: never[]; issues: never[] } {
  const entry: AssemblyEntry = {
    address: at.address,
    names: [],
    ...(source ? { source } : {}),
    ...(aliasSource ? { aliasSource } : {}),
    pending: true,
    sequences: [],
    order: [],
    systems: [],
    systemAt: [],
    lookup: { index: new Map(), match: 'chr', dropped: [] },
    find: () => undefined,
    include: () => undefined,
  };
  return { entry, systems: [], issues: [] };
}
