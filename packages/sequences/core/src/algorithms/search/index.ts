/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-core/algorithms/search`: matchers and a text index for applications to compose their own search from
// (seqspec:libraries#search). Platform-free; no dependencies. The window helpers of `matchers.ts` stay internal.

export { type ApproximateEnds, type ApproximateOptions, approximateMatcher } from './approximate.js';
export {
  classMatcher,
  compileMotif,
  exactMatcher,
  type HitVisitor,
  MAX_CLASS_PATTERN,
  type Matcher,
  toCodes,
} from './matchers.js';
export { regexMatcher } from './regex.js';
export { type TextIndex, type TextIndexOptions, textIndex } from './text-index.js';
