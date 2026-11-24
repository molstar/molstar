/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/search (_spec/README.md "Search"): the matcher port, the kernels and `kernelMatcher`, the built-in
// matchers, the bounds and the codes. Named re-exports only.

export {
  BuiltInMatcherNames,
  BuiltInSearchRegistry,
  checkMatcherNames,
  type LazySearchMatcher,
  loadMatcher,
  type SearchMatcherEntry,
} from './entries.js';
export { HitList, NoHits } from './hits.js';
export { IncompleteSearchCodes, SearchIssues } from './issues.js';
export { type KernelMatcherSpec, kernelMatcher } from './kernel-matcher.js';
export { classKernel, classMasks, exactKernel } from './kernels.js';
export { BuiltInSearchMatchers, DegenerateMatcher, ExactMatcher, RegexMatcher } from './matchers.js';
export { glushkov, parseRegex, type RegexParsed, regexKernel, regexTables } from './regex.js';
export {
  compileFailed,
  DefaultSearchBounds,
  EdgeEnd,
  EdgeStart,
  type HitVisitor,
  type LetterReader,
  type SearchBatch,
  type SearchBounds,
  type SearchCompiled,
  type SearchHits,
  type SearchKernel,
  type SearchMatcher,
  type SearchOptions,
  type SearchRunContext,
  type SearchStrands,
  type SearchTarget,
  type StepBudget,
  searchBounds,
} from './types.js';
