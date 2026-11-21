/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Fields } from '@molstar/core/repr/schema';
import type { Issue, Node } from '@molstar/core/repr/tree';
import type { SchemeRegistry } from '@molstar/seqx-core/theme';
import type { LayerKinds } from '@molstar/seqx-core/view';
import type { SeqConfigItem, SeqRegistryEntry } from './context.js';

// Extensions (_spec/README.md "Context", Extensions): a descriptor, data plus providers, given at composition beside
// the registry. Its id is its namespace: every entry it ships is named `<id>:<name>`; its `requires` covers the
// language's version; registration is the registry's, the descriptors whole. Its checks are `defineExtension`'s (an extension's entry); every composition bundles this module.

/** The tables an extension adds to past its registry entries: the layer kinds, the schemes. */
export interface ExtensionTables {
  readonly kinds: LayerKinds;
  readonly schemes?: SchemeRegistry;
}

/** An extension: what an outside package ships (layer types, layouts, sizings, navigator pictures, schemes, formats, query entries). */
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
   * The params it reads, by the type that reads them: its own types' and pictures' plain params
   * (`'rondo:ribbons': { bundle }`), and under a built-in kind the namespaced keys it reads on others' nodes
   * (`section: { 'rondo:bundle' }`); checked by `check` (`defineExtension` of `@molstar/seqx-core/extension` writes it from these).
   */
  readonly params?: Readonly<Record<string, Fields>>;
  /** Its params' issues on a spec: what the context's `check` adds where the extension is composed. */
  readonly check?: (spec: Node) => Issue[];
  /** Its settings, keys under its namespace (`rondo.maxRibbons`). */
  readonly config?: readonly SeqConfigItem<any>[];
  /**
   * Written by `defineExtension` (`@molstar/seqx-core/extension`), which checked the descriptor: what it adds to the
   * context's tables past its registry entries, so a composition with no extension bundles none of it.
   */
  readonly defined?: {
    /** The tables with its parts added: its layer types' kinds (and its version reading), its schemes. */
    compose(tables: ExtensionTables): ExtensionTables;
  };
}

/**
 * The extensions checked at composition (_spec/README.md "Context", Extensions): each made by `defineExtension` (which
 * checked its namespace and its `requires`), the same descriptor twice one, another under one id an error. Small: every
 * composition bundles it; the checks of a descriptor are its maker's.
 */
export function checkExtensions(extensions: readonly SeqExtension[]): SeqExtension[] {
  const out = [...new Set(extensions)];
  const bad = out.filter((e, i) => !e.defined || out.findIndex((x) => x.id === e.id) < i).map((e) => e.id);
  if (bad.length) throw new Error(`seq.context.extension.namespace: not made by defineExtension, or twice: ${bad}`);
  return out;
}
