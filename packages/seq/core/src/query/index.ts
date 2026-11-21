/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/seqx-core/query: locations, element keys, the three result types and the first query forms as `repr` form
// resolvers over a model of systems, numberings and annotations (_spec/README.md "Queries"); addresses resolved and printed back
// ("Addressing").

export {
  type AddressBounds,
  type AddressForm,
  type AddressForOptions,
  type AddressResolution,
  type AddressStatus,
  addressFor,
  DefaultAddressBounds,
  DefaultAddressPreference,
  type Located,
  locationOfResult,
  type ResolveAddressOptions,
  refResults,
  resolveAddressIn,
} from './address.js';
export { AddressIssues } from './address-issues.js';
export {
  DefaultLinkOptions,
  type Link,
  type LinkAction,
  type LinkActionName,
  LinkActions,
  type LinkOptions,
  parseLink,
  printLink,
} from './links.js';
export { EmptyKeys, type Key, type Keys, keySet, keysEqual, unionKeys } from './keys.js';
export { loadGenomeForms } from './genome-door.js';
export { type Evaluation, evaluateQuery, type QueryOptions, QueryFormKinds, queryResolvers } from './forms.js';
export type { ExprValue } from './expressions.js';
export {
  type QueryExtensionEntry,
  type QueryExtensions,
  type QueryFormContext,
  type QueryFormEntry,
  type QueryFormRole,
  type QueryOperatorEntry,
  type QueryPropertyEntry,
  queryExtensions,
} from './extensions.js';
// a selection's items, renamed apart from render's `readItems` and the overlay's `itemsOf` (the builder's JSON reads
// them; addressing.md, Schema selections)
export {
  type Items as SelectionItems,
  type ItemsInput as SelectionItemsInput,
  itemsOf as selectionItemsOf,
  readItems as readSelectionItems,
} from './items.js';
export type { PropertyDef, PropertyLevel, PropertyTable, PropertyType } from './properties.js';
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
export type {
  ItemsRead,
  Numbering,
  QueryAnnotation,
  PropertyColumn,
  PropertyRead,
  QueryAtoms,
  QueryComponents,
  QueryItem,
  QueryModel,
  QueryMolecule,
  QuerySystem,
  QueryUnit,
} from './model.js';
export {
  anyLabelIn,
  BuiltInNumberings,
  DefaultNumbering,
  findNumbering,
  HintScanLimit,
  labelIndex,
  labelsNumbering,
  type NumberingOf,
  nearestLabels,
  numberingIds,
  numberingPending,
  sequentialNumbering,
  snapLabel,
  type Tick,
} from './numbering.js';
export { fromNamed, namedType, type QueryResult, type QueryResultType, toNamed } from './result.js';
