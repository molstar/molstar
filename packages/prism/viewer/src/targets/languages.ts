/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Expression } from '@molstar/query-language/language/expression';
import type { QueryLanguage } from './types.js';

// The query languages target text may be written in (targets.md §1): each an entry a viewer's spec lists
// (`PrismViewerSpec.languages`), so a viewer bundles only the transpilers it names (`targets/languages/*`, one module
// each). Registered when a viewer is made; a language no viewer registered is not read (keywords and address lists
// always are).

export interface QueryLanguageEntry {
  readonly id: QueryLanguage;
  /** Its name for people ("PyMOL"). */
  readonly name: string;
  /** What its selections look like, for a failure's hint. */
  readonly hint: string;
  /** The text as MolQL; throws with the parser's message. */
  transpile(source: string): Expression;
}

const registered = new Map<QueryLanguage, QueryLanguageEntry>();
// The registered entries as a list, made again only when one is added: params read it per call (`LanguageParam`).
let list: readonly QueryLanguageEntry[] = [];

/** Makes languages readable (a viewer's spec's, when it is made); again is harmless. */
export function registerQueryLanguages(languages: readonly QueryLanguageEntry[]): void {
  for (const l of languages) registered.set(l.id, l);
  if (registered.size !== list.length) list = [...registered.values()];
}

export function queryLanguage(id: QueryLanguage): QueryLanguageEntry | undefined {
  return registered.get(id);
}

/** The languages read, in the order the spec lists them. */
export function queryLanguages(): readonly QueryLanguageEntry[] {
  return list;
}

/** The language of a text that names none: the first the viewer reads (PyMOL in the default spec). */
export function defaultQueryLanguage(): QueryLanguage {
  return list[0]?.id ?? 'pymol';
}
