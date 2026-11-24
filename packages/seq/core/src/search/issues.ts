/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { defineIssues } from '@molstar/core/repr/tree';

/** The search's codes (_spec/README.md "Search", Bounds): each bound refused or said with its own. */
export const SearchIssues = defineIssues('seq.search', {
  'seq.search.pattern.size': {
    severity: 'error',
    description:
      'A pattern is past its bound: an exact or IUPAC pattern’s letters, a regular expression’s characters, positions or nesting',
    hint: 'Search for a shorter pattern',
  },
  'seq.search.pattern.unsupported': {
    severity: 'error',
    description:
      'A regular expression uses what the search does not run: backreferences, lookaround, lazy or possessive repeats, \\b, \\d, \\w, flags',
    hint: 'A search runs in linear time: no backreferences or lookaround',
  },
  'seq.search.pattern.syntax': {
    severity: 'error',
    description: 'A regular expression does not read: an unclosed group or class, a repeat of nothing, a bad count',
  },
  'seq.search.pattern.empty': {
    severity: 'error',
    description: 'A pattern matches nothing but the empty string, or is empty',
    hint: 'Type the letters to search for',
  },
  'seq.search.pattern.anchor': {
    severity: 'warning',
    description: 'A regular expression anchored with ^ or $ is not searched on a cyclic system, which has no ends',
  },
  'seq.search.alphabet': {
    severity: 'warning',
    description: 'A system’s alphabet cannot hold the pattern’s letters: it is skipped, the others searched',
  },
  'seq.search.letters.limit': {
    severity: 'warning',
    description: 'A system would take the letters a search reads past their bound: it is skipped',
    hint: 'Search within a region (`within`)',
  },
  'seq.search.letters.failed': {
    severity: 'warning',
    description: 'A system’s letters could not be read: it is skipped',
  },
  'seq.search.steps.limit': {
    severity: 'warning',
    description: 'A search took the work it may take (its word steps): it ended there, its result incomplete',
    hint: 'Search for a narrower pattern, or within a region',
  },
  'seq.search.hits.limit': {
    severity: 'warning',
    description: 'A search found more hits than it holds: the count goes on, the first ones are held',
    hint: 'Search for a narrower pattern, or within a region',
  },
  'seq.search.matcher.missing': {
    severity: 'error',
    description: 'A search names a matcher the composition does not register',
    hint: 'Register it as a `search-matcher` entry',
  },
  'seq.search.save.limit': {
    severity: 'error',
    description: 'A result has more hits than a saved set holds: not saved',
    hint: 'Search for a narrower pattern, or within a region',
  },
  'seq.search.save.incomplete': {
    severity: 'error',
    description:
      'A result is incomplete (a bound reached, a system skipped, cancelled): a saved set must mean what it says',
  },
  'seq.search.save.rects': {
    severity: 'error',
    description: 'A result on an alignment’s rows is rectangles, which a saved set cannot name yet',
  },
});

/** The codes that leave a result incomplete. */
export const IncompleteSearchCodes: ReadonlySet<string> = new Set([
  'seq.search.letters.limit',
  'seq.search.letters.failed',
  'seq.search.steps.limit',
  'seq.search.hits.limit',
]);
