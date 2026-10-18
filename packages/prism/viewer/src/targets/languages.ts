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
// each). Registered when a viewer is made, on that viewer's own set; a language a viewer does not list is not read
// (keywords and address lists always are).

export interface QueryLanguageEntry {
  readonly id: QueryLanguage;
  /** Its name for people ("PyMOL"). */
  readonly name: string;
  /** What its selections look like, for a failure's hint. */
  readonly hint: string;
  /** The text as MolQL; throws with the parser's message. */
  transpile(source: string): Expression;
}

/**
 * The languages one viewer reads (README decision 29): its own set, so two viewers on a page can read different ones.
 * Seeded from the page-wide defaults (`registerQueryLanguages`, below), then the spec's.
 */
export class QueryLanguages {
  readonly #registered = new Map<QueryLanguage, QueryLanguageEntry>();
  // The registered entries as a list, made again only when one is added: params read it per call (`LanguageParam`).
  #list: readonly QueryLanguageEntry[] = [];

  constructor(languages: readonly QueryLanguageEntry[] = []) {
    this.register(languages);
  }

  /** Makes languages readable; again is harmless. */
  register(languages: readonly QueryLanguageEntry[]): void {
    for (const l of languages) this.#registered.set(l.id, l);
    if (this.#registered.size !== this.#list.length) this.#list = [...this.#registered.values()];
  }

  get(id: QueryLanguage): QueryLanguageEntry | undefined {
    return this.#registered.get(id);
  }

  /** The languages read, in the order they were registered (the spec lists them). */
  get list(): readonly QueryLanguageEntry[] {
    return this.#list;
  }

  /** The language of a text that names none: the first the viewer reads (PyMOL in the default spec). */
  get default(): QueryLanguage {
    return this.#list[0]?.id ?? 'pymol';
  }
}

// The page-wide defaults: what a viewer starts with, and what target functions called outside a viewer read. A viewer
// registers its spec's languages on its own set, never here.
const pageWide = new QueryLanguages();

/** Makes languages readable page-wide: the seed of every viewer made after, and of target functions called with no viewer. */
export function registerQueryLanguages(languages: readonly QueryLanguageEntry[]): void {
  pageWide.register(languages);
}

/** A viewer's set: the page-wide defaults, then its spec's languages. */
export function createQueryLanguages(languages: readonly QueryLanguageEntry[] = []): QueryLanguages {
  return new QueryLanguages([...pageWide.list, ...languages]);
}

/** The page-wide set: what a target function reads when no viewer's set is given. */
export function pageQueryLanguages(): QueryLanguages {
  return pageWide;
}

export function queryLanguage(id: QueryLanguage): QueryLanguageEntry | undefined {
  return pageWide.get(id);
}

/** The page-wide languages, in registration order. */
export function queryLanguages(): readonly QueryLanguageEntry[] {
  return pageWide.list;
}

/** The page-wide default language. */
export function defaultQueryLanguage(): QueryLanguage {
  return pageWide.default;
}
