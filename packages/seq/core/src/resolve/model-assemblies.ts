/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model's assemblies (_spec/README.md "Resolve", Assemblies): each assembly's entry and derived systems, built by
// the genome chunk once it lands (pending until then, a failed load said once), and its sources' coverage carried to
// the next resolve inside the `Lengths` (opaque: the systems ride along, the key's value their content's hash).

import { hashString } from '@molstar/core/data/util/hash-functions';
import { type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import type { Coverage, SystemInfo } from '@molstar/seqx-core/data';
import { type AssemblyEntry, genomeLoaded, genomeLoadFailure, loadGenome, pendingAssembly } from './genome.js';
import { ResolveIssues } from './issues.js';
import type { Cached, Placed } from './model-molecules.js';
import type { Lengths, SystemEntry } from './model.js';

/** A length's key: the system and the source node it was read from. */
export const lengthKey = (system: string, sourceHash: string) => `${system}|${sourceHash}`;

/** The systems an assembly's sources listed, by their key in a `Lengths` (opaque: they ride along with it). */
const Coverages = new WeakMap<Lengths, ReadonlyMap<string, readonly SystemInfo[]>>();
/** A coverage's systems hashed once, when the source opens (its names, lengths and aliases). */
const ContentHashes = new WeakMap<readonly SystemInfo[], number>();

function contentHash(systems: readonly SystemInfo[]): number {
  let h = ContentHashes.get(systems);
  if (h === undefined) {
    h = systems.length;
    for (const s of systems) h = (h * 31 + hashString(`${s.name}\t${s.length}\t${s.aliases?.join(',') ?? ''}`)) | 0;
    ContentHashes.set(systems, h);
  }
  return h;
}

/** The assemblies' sources' coverage into `out` (what `sourceLengths` returns): each keyed by its content's hash. */
export function assemblyLengths(
  assemblies: Iterable<AssemblyEntry>,
  hashes: ReadonlyMap<string, string> | undefined,
  coverage: (source: string) => Coverage | undefined,
  out: Map<string, number>,
): void {
  const listed = new Map<string, readonly SystemInfo[]>();
  for (const a of assemblies)
    for (const src of [a.source, a.aliasSource]) {
      const hash = src && hashes?.get(src);
      const systems = hash ? coverage(src!)?.systems : undefined;
      if (!systems) continue;
      const k = lengthKey(`${a.address}#${src}`, hash!);
      out.set(k, contentHash(systems));
      listed.set(k, systems);
    }
  if (listed.size) Coverages.set(out, listed);
}

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** What resolving the assemblies reads of the build. */
export interface AssemblyBuild {
  readonly data: readonly Placed[];
  readonly ref: (at: Placed, key: string, kinds: readonly string[]) => Placed | undefined;
  readonly sourceHashes: ReadonlyMap<string, string>;
  readonly given: Lengths;
  readonly cached: Cached;
  readonly key: (...at: Placed[]) => string;
  readonly systems: Map<string, SystemEntry>;
  readonly issues: Issue[];
}

/**
 * The assemblies by address, their derived systems added to `systems`, and the genome chunk's load when a spec holds an
 * assembly and the chunk has not landed (each assembly pending meanwhile; a failed load said once, never retried by a
 * resolve).
 */
export function resolveAssemblies(b: AssemblyBuild): {
  assemblies: Map<string, AssemblyEntry>;
  loading?: Promise<unknown>;
} {
  const assemblies = new Map<string, AssemblyEntry>();
  let loading: Promise<unknown> | undefined;
  let said = false;
  const listed = Coverages.get(b.given);
  for (const at of b.data) {
    if (at.node.kind !== 'assembly') continue;
    const src = p(at.node, 'source') === undefined ? undefined : b.ref(at, 'source', ['source'])?.address;
    const alias = p(at.node, 'alias_source') === undefined ? undefined : b.ref(at, 'alias_source', ['source'])?.address;
    const coverageOf = (s: string | undefined) => {
      const h = s && b.sourceHashes.get(s);
      const k = h ? lengthKey(`${at.address}#${s}`, h) : undefined;
      return k ? { systems: listed?.get(k), hash: b.given.get(k) } : undefined;
    };
    const c = coverageOf(src);
    const t = coverageOf(alias);
    const genome = genomeLoaded();
    const failure = genome ? undefined : genomeLoadFailure();
    if (!genome && failure === undefined) loading ??= loadGenome();
    if (failure !== undefined && !said) {
      said = true;
      b.issues.push(
        issueOf(ResolveIssues, 'seq.resolve.genome.load', `The genome chunk did not load: ${String(failure)}`, {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.genome.load'].hint,
        }),
      );
    }
    const e = genome
      ? b.cached(`${at.address}|${b.key(at)}|${c?.hash ?? ''}|${t?.hash ?? ''}`, () =>
          genome.assemblyEntry(at, src, alias, c?.systems, t?.systems),
        )
      : pendingAssembly(at, src, alias);
    assemblies.set(at.address, e.entry);
    for (const s of e.systems) b.systems.set(s.address, s);
    b.issues.push(...e.issues);
  }
  return loading ? { assemblies, loading } : { assemblies };
}
