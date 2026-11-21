/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The model's assemblies in the genome chunk (_spec/README.md "Resolve", Assemblies): each assembly's entry and
// derived systems, and its sources' coverage carried to the next resolve inside the `Lengths` (opaque: the systems ride
// along, the key's value their content's hash, computed once when the source opens).

import { hashString } from '@molstar/core/data/util/hash-functions';
import type { Node, Value } from '@molstar/core/repr/tree';
import type { Coverage, SystemInfo } from '@molstar/seqx-core/data';
import { assemblyEntry } from './assemblies.js';
import type { AssemblyEntry } from './genome.js';
import { type AssemblyBuild, lengthKey } from './model-assemblies.js';
import type { Lengths } from './model.js';

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

/** The assemblies by address, their derived systems added to the build's. */
export function resolveAssemblies(b: AssemblyBuild): Map<string, AssemblyEntry> {
  const assemblies = new Map<string, AssemblyEntry>();
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
    const e = b.cached(`${at.address}|${b.key(at)}|${c?.hash ?? ''}|${t?.hash ?? ''}`, () =>
      assemblyEntry(at, src, alias, c?.systems, t?.systems),
    );
    assemblies.set(at.address, e.entry);
    for (const s of e.systems) b.systems.set(s.address, s);
    b.issues.push(...e.issues);
  }
  return assemblies;
}
