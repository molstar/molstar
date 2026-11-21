/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Fields } from '@molstar/core/repr/schema';
import type { Issue, Node } from '@molstar/core/repr/tree';
import { SeqLanguageVersion } from '@molstar/seqx-core/lang';
import { versionParts } from '@molstar/seqx-core/view';
import type { SeqConfigItem, SeqRegistryEntry } from './context.js';

// Extensions (_spec/README.md "Context", Extensions): a descriptor, data plus providers, given at composition beside
// the registry. Its id is its namespace: every entry it ships is named `<id>:<name>`; its `requires` covers the
// language's version; registration is the registry's, the descriptors whole. Kept small: every composition bundles it.

/** An extension: what an outside package ships (layer types, layouts, sizings, schemes, formats, query entries). */
export interface SeqExtension {
  /** Its namespace: a lowercase word (`rondo`), never `seq` (reserved). */
  readonly id: string;
  /** The package's version (semver): what a spec's `extensions` is read against. */
  readonly version: string;
  /** The language versions its contract targets: `0.5` (0.x: the minor is the major), `1`, `1.2`. */
  readonly requires: string;
  readonly doc?: string;
  /** Its entries, each named `<id>:<name>`. */
  readonly registry: readonly SeqRegistryEntry[];
  /**
   * The params it reads, by the type that reads them: its own types' plain params (`'rondo:ribbons': { bundle }`), and
   * under a built-in kind the namespaced keys it reads on others' nodes (`section: { 'rondo:bundle' }`); checked by
   * `check` (`defineExtension` of `@molstar/seqx-core/extension` writes it from these).
   */
  readonly params?: Readonly<Record<string, Fields>>;
  /** Its params' issues on a spec: what the context's `check` adds where the extension is composed. */
  readonly check?: (spec: Node) => Issue[];
  /** Its settings, keys under its namespace (`rondo.maxRibbons`). */
  readonly config?: readonly SeqConfigItem<any>[];
}

/** Whether `requires` covers a language version: its major (0.x: its minor), at least its minor (0.x: its patch). */
export function requiresCovers(requires: string, version = SeqLanguageVersion): boolean {
  const [a, b] = versionParts(requires);
  const [la, lb] = versionParts(version);
  return la === a && lb >= b;
}

/**
 * The extensions checked (_spec/README.md "Context", Extensions): the same descriptor twice is one, another under one
 * id an error; every entry and setting named in its namespace (`seq.context.extension.namespace`, every offender listed
 * in one error); its `requires` covering the language (`seq.context.extension.version`). Throws before anything is
 * registered.
 */
export function checkExtensions(extensions: readonly SeqExtension[]): SeqExtension[] {
  const byId = new Map<string, SeqExtension>();
  const bad: string[] = [];
  const old: string[] = [];
  for (const e of extensions) {
    if ((byId.get(e.id) ?? e) !== e) bad.push(`${e.id} (twice)`);
    byId.set(e.id, e);
  }
  for (const e of byId.values()) {
    const ns = (name: string) => name.startsWith(`${e.id}:`) || bad.push(`${e.id}: ${name}`);
    if (!/^[a-z][a-z0-9_-]*$/.test(e.id) || e.id === 'seq') bad.push(`${e.id} (its id)`);
    // each entry's provider by its name (a format runtime's, a layer type's, a scheme's; a query entry's own)
    for (const r of e.registry as any[])
      ns((r.runtime ?? r.layout ?? r.sizing ?? r.scheme ?? (r.kind === 'layer-type' ? r.type : r)).name);
    for (const c of e.config ?? []) ns(c.key.replace('.', ':'));
    if (!requiresCovers(e.requires)) old.push(`${e.id} requires ${e.requires}`);
  }
  const fail = (code: string, list: string[]) => {
    if (list.length) throw Object.assign(new Error(`${code}: ${list.join(', ')}`), { code });
  };
  fail('seq.context.extension.namespace', bad);
  fail('seq.context.extension.version', old);
  return [...byId.values()];
}
