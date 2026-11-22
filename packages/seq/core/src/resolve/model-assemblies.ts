/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model's assemblies, the eager part (_spec/README.md "Resolve", Assemblies): with the genome chunk landed, its
// `resolveAssemblies`; before it, each assembly pending (its sources named, nothing derived) and the chunk's load the
// model waits on; a failed load said once, never retried by a resolve.

import { exactAddress, type Issue, issueOf, type Node, type Value } from '@molstar/core/repr/tree';
import type { Coverage } from '@molstar/seqx-core/data';
import { type AssemblyEntry, genomeLoaded, genomeLoadFailure, loadGenome, pendingAssembly } from './genome.js';
import { ResolveIssues } from './issues.js';
import type { Cached, Placed } from './model-molecules.js';
import type { AnnotationEntry, Lengths, SystemEntry } from './model.js';
import { statesSystemOf } from './model-region.js';

/** A length's key: the system and the source node it was read from. */
export const lengthKey = (system: string, sourceHash: string) => `${system}|${sourceHash}`;

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

const p = (n: Node, key: string): Value | undefined => n.params?.[key];

/** The assemblies by address (their systems added to the build's), and the genome chunk's load while it is awaited. */
export function resolveAssemblies(b: AssemblyBuild): {
  assemblies: Map<string, AssemblyEntry>;
  loading?: Promise<unknown>;
} {
  const genome = genomeLoaded();
  if (genome) return { assemblies: genome.resolveAssemblies(b) };
  const assemblies = new Map<string, AssemblyEntry>();
  const failure = genomeLoadFailure();
  for (const at of b.data) {
    if (at.node.kind !== 'assembly') continue;
    const src = p(at.node, 'source') === undefined ? undefined : b.ref(at, 'source', ['source'])?.address;
    const alias = p(at.node, 'alias_source') === undefined ? undefined : b.ref(at, 'alias_source', ['source'])?.address;
    if (failure !== undefined && !assemblies.size)
      b.issues.push(
        issueOf(ResolveIssues, 'seq.resolve.genome.load', `The genome chunk did not load: ${String(failure)}`, {
          address: at.address,
          hint: ResolveIssues.codes['seq.resolve.genome.load'].hint,
        }),
      );
    assemblies.set(at.address, pendingAssembly(at, src, alias).entry);
  }
  return assemblies.size && failure === undefined ? { assemblies, loading: loadGenome() } : { assemblies };
}

/**
 * The assemblies' sources' coverage into `out` (`sourceLengths`), once the genome chunk has landed; with `data`, the
 * sources of data over an assembly too.
 */
export function assemblyLengths(
  assemblies: Iterable<AssemblyEntry>,
  hashes: ReadonlyMap<string, string> | undefined,
  coverage: (source: string) => Coverage | undefined,
  out: Map<string, number>,
  data?: Iterable<{ readonly address: string; readonly source: string }>,
): void {
  genomeLoaded()?.assemblyLengths(assemblies, hashes, coverage, out, data);
}

/** The kinds of data read over an assembly by sequence (_spec/README.md "Resolve", Assemblies). */
const OverAssemblyKinds = new Set(['ranges', 'values', 'pairs', 'grid', 'matrix']);

/** An assembly entry's identity as a cache key's part (a new entry, a new key). */
const AssemblyKeys = new WeakMap<AssemblyEntry, string>();
let assemblyKeys = 0;
function assemblyKey(a: AssemblyEntry | undefined): string {
  if (!a) return '';
  let k = AssemblyKeys.get(a);
  if (k === undefined) AssemblyKeys.set(a, (k = String(++assemblyKeys)));
  return k;
}

/**
 * Data over an assembly into `annotations` (and a matrix's states into the build's systems), once the genome chunk has
 * landed: a plan a derived sequence (`genome-data.ts`), the names its source lists carried in the lengths as an
 * assembly's own source's are. Kept while the node, its assemblies' entries and what its source listed are the same.
 */
export function assemblyAnnotations(
  b: AssemblyBuild,
  assemblies: ReadonlyMap<string, AssemblyEntry>,
  lengths: ReadonlyMap<string, number>,
  annotations: Map<string, AnnotationEntry>,
): void {
  const genome = genomeLoaded();
  if (!genome) return;
  const at = (n: Node, k: string) => assemblies.get(exactAddress(p(n, k)) ?? '');
  for (const d of b.data) {
    const n = d.node;
    if (!OverAssemblyKinds.has(n.kind)) continue;
    const over = at(n, 'target') ?? at(n, 'a') ?? at(n, 'x');
    if (!over) continue;
    const second = at(n, 'b') ?? at(n, 'y');
    const src = p(n, 'source') === undefined ? undefined : b.ref(d, 'source', ['source'])?.address;
    const listed = genome.listedSystems(b.given, b.sourceHashes, d.address, src);
    const states = lengths.get(statesSystemOf(d.address));
    const e = b.cached(
      `${d.address}|${b.key(d)}|${assemblyKey(over)}|${assemblyKey(second)}|${listed?.hash ?? ''}|${states ?? ''}`,
      () =>
        genome.assemblyAnnotation(d, {
          assemblies,
          systems: b.systems,
          ...(src !== undefined ? { source: src } : {}),
          ...(listed?.systems ? { listed: listed.systems } : {}),
          ...(states !== undefined && states >= 0 ? { states } : {}),
        }),
    );
    if (!e) continue;
    annotations.set(d.address, e.annotation);
    if (e.states) b.systems.set(e.states.address, e.states);
    b.issues.push(...e.issues);
  }
}
