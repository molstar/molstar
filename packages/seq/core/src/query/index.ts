/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/query: locations, element keys, the three result types and the first query forms as `repr` form
// resolvers over a model of systems, numberings and annotations (_spec/README.md "Queries").

export {
  type ElementKey,
  type Elements,
  EmptyElements,
  compareElementKeys,
  elements,
  elementsEqual,
  hashElements,
  intersectElements,
  subtractElements,
  unionElements,
} from './elements.js';
export { type Evaluation, evaluateQuery, type QueryOptions, QueryFormKinds, queryResolvers } from './forms.js';
export { type QueryIssueCode, QueryIssues } from './issues.js';
export {
  compareLabels,
  hasLabelSuffix,
  isNoInsertionCode,
  type LabelRead,
  type ParsedLabel,
  parseLabel,
  readLabel,
} from './labels.js';
export {
  compareRects,
  EmptyLocation,
  hashLocation,
  intersectLocations,
  isEmptyLocation,
  isPosition,
  type Location,
  type LocationPart,
  type LocationRect,
  location,
  locationsEqual,
  normalizeRanges,
  overlapsRanges,
  type PartInput,
  partOf,
  positionCount,
  type RectInput,
  rangesLocation,
  subtractLocations,
  unionLocations,
} from './location.js';
export type { ItemsRead, Numbering, QueryAnnotation, QueryItem, QueryModel, QuerySystem } from './model.js';
export {
  anyLabelIn,
  BuiltInNumberings,
  DefaultNumbering,
  findNumbering,
  HintScanLimit,
  labelIndex,
  type NumberingOf,
  nearestLabels,
  numberingIds,
  sequentialNumbering,
  snapLabel,
  type Tick,
} from './numbering.js';
export { fromNamed, namedType, type QueryResult, type QueryResultType, toNamed } from './result.js';
