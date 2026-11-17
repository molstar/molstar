/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @molstar/core/repr/tree: the node and its values (typed columns among them), identity (the Merkle hash, equality),
// the spelling of atoms, addresses with query segments and the form-resolver table, and the coded issue
// (_spec/README.md). Codec-independent: the codecs are `repr/codecs/*`.

export {
  type Address,
  type AddressPattern,
  addressOf,
  addressWordEnd,
  CancelEvery,
  exactAddress,
  type FormResolver,
  type FormResolverInput,
  formEnd,
  isCanonicalAddress,
  isExactAddress,
  matchesAddress,
  matchesForm,
  type Named,
  parseAddress,
  parseAddressIn,
  printAddress,
  printAddressAtom,
  resolveAddress,
  type Segment,
  type Selected,
  type Selection,
  type SelectOptions,
  selectAddress,
  unknownFormKey,
} from './address.js';
export {
  Delimiter,
  isBareString,
  KeywordPattern,
  printAtom,
  printKey,
  readQuoted,
  readWord,
  WordPattern,
} from './atoms.js';
export {
  CategoryColumn,
  type Column,
  type ColumnType,
  ColumnTypes,
  categoryColumn,
  canonicalColumn,
  categoryFromCodes,
  columnConstructor,
  columnMisfit,
  columnsEqual,
  columnType,
  columnValues,
  hashColumn,
  integerRange,
  isColumn,
  isFloatType,
  type NumericColumn,
  type NumericColumnType,
  toColumn,
} from './column.js';
export { hashBytes, hashNode, nodesEqual, valuesEqual } from './hash.js';
export {
  defineIssues,
  hasErrors,
  type Issue,
  type IssueCatalog,
  type IssueCodeInfo,
  type IssuePlace,
  issueOf,
  pathToPointer,
  type Severity,
  TreeIssues,
  type TreePath,
} from './issues.js';
export {
  type Atom,
  assertKeys,
  assertPrintable,
  IdHint,
  IdPattern,
  isNode,
  isRecord,
  isReservedKey,
  KindHint,
  KindPattern,
  kindNamespace,
  MaxDepth,
  type Node,
  own,
  type Params,
  type Value,
  type ValueRecord,
  visitTree,
} from './node.js';
