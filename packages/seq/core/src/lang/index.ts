/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/lang: SeqViewSpec as a description on `@molstar/core/repr` (the kinds, their params), the
// language builder, the `seq.lang.*` issue catalog and the language-level checks (_spec/README.md).

export {
  type Address,
  type Addressed,
  createSeqSpec,
  type DataBuilder,
  type Params,
  type RefValue,
  type SectionBuilder,
  type SeqSpecBuilder,
  type SeqSpecOptions,
  type TrackBuilder,
  type ViewBuilder,
  type ViewsBuilder,
} from './builder.js';
export { checkDataNode, checkSeqLang, checkSeqSpec, type DataNodeTargets } from './check.js';
export { LangIssues } from './issues.js';
export { type SeqKinds, SeqLanguage, type SeqLanguageType, SeqLanguageVersion } from './language.js';
