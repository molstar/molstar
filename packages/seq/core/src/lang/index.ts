/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/lang: SeqViewSpec as a description on `@molstar/core/repr` (the kinds, their params), the
// language builder, the `seq.lang.*` issue catalog, the language-level checks and user-defined refs (_spec/README.md).

export {
  type Address,
  type Addressed,
  createSeqSpec,
  type DataBuilder,
  ideogramNavigator,
  miniatureNavigator,
  type NavigatorOf,
  type NavigatorParams,
  overviewNavigator,
  type Params,
  type RefsBuilder,
  type RefTo,
  type RefValue,
  type SectionBuilder,
  type SequenceBuilder,
  type SeqSpecBuilder,
  type SeqSpecOptions,
  type SnapshotsBuilder,
  stripNavigator,
  type TrackBuilder,
  type ViewBuilder,
  type ViewsBuilder,
} from './builder.js';
export {
  checkSpecBounds,
  checkSubtreeDepth,
  DefaultSpecBounds,
  measureNode,
  queryDepth,
  type SpecBounds,
  specBounds,
} from './bounds.js';
export { checkDataNode, checkSeqLang, checkSeqSpec, type DataNodeTargets } from './check.js';
export { BoundIssues, LangIssues } from './issues.js';
export {
  type SeqKinds,
  SeqLanguage,
  type SeqLanguageType,
  SeqLanguageVersion,
  SeqPatchLanguage,
  SeqValidateOptions,
  seqLanguages,
} from './language.js';
export { type FieldOptions, measure, theme } from './theme-builder.js';
export { ScaleParamsDescription, ThemeLanguage, ThemeLanguageId } from './theme-language.js';
export {
  addRefPatch,
  DefaultRefDepth,
  type FollowedRef,
  followRef,
  type RefInput,
  type RefRefusal,
  RefsId,
  refAddresses,
  refPrefix,
} from './refs.js';
