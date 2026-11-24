/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `search-matcher` registry entries (_spec/README.md "Context", Registry entries; "Search", Matchers): a matcher, or
// one whose code loads on first use (the built-ins: a lazy chunk, so a composition that never searches loads none of
// it). The built-ins are unnamespaced; an application's are `<id>:<name>`.

import type { SearchMatcher } from './types.js';

/** A matcher whose code is a lazy chunk: its name and title at once, the matcher when loaded. */
export interface LazySearchMatcher {
  readonly name: string;
  readonly title: string;
  readonly alphabets?: SearchMatcher['alphabets'];
  load(): Promise<SearchMatcher<any>>;
}

/** What a `search-matcher` entry holds. */
export type SearchMatcherEntry = SearchMatcher<any> | LazySearchMatcher;

/** The built-ins' names: never shadowed, never namespaced. */
export const BuiltInMatcherNames: readonly string[] = ['degenerate', 'exact', 'regex'];

const lazy = (name: 'exact' | 'degenerate' | 'regex', title: string): LazySearchMatcher => ({
  name,
  title,
  load: () => import('./matchers.js').then((m) => m.BuiltInSearchMatchers[name]),
});

/** The built-in matchers as registry entries, in the order a choice lists them (the default first). */
export const BuiltInSearchRegistry: readonly {
  readonly kind: 'search-matcher';
  readonly matcher: LazySearchMatcher;
}[] = [
  { kind: 'search-matcher', matcher: lazy('degenerate', 'IUPAC') },
  { kind: 'search-matcher', matcher: lazy('exact', 'Exact') },
  { kind: 'search-matcher', matcher: lazy('regex', 'Regular expression') },
];

/** An entry's matcher, loaded when lazy. */
export function loadMatcher(entry: SearchMatcherEntry): Promise<SearchMatcher<any>> {
  return 'load' in entry ? entry.load() : Promise.resolve(entry);
}

const Namespaced = /^[a-z][a-z0-9-]*:[a-zA-Z][\w-]*$/;

/**
 * The registry's matchers checked: an application's namespaced (`<id>:<name>`), so no built-in is shadowed. Throws
 * naming every offender.
 */
export function checkMatcherNames(matchers: readonly SearchMatcherEntry[]): void {
  const builtIns = new Set<SearchMatcherEntry>(BuiltInSearchRegistry.map((e) => e.matcher));
  const bad = matchers.filter((m) => !builtIns.has(m) && !Namespaced.test(m.name)).map((m) => m.name);
  if (bad.length)
    throw new Error(`A search matcher of an application is namespaced (\`<id>:<name>\`): ${bad.join(', ')}`);
}
