/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seqx-core/data`: the source protocol, the `seq` format and its runtime, and the data store
// (_spec/README.md "Data"). Named re-exports only.

export { SourceIssues, StoreIssues } from './issues.js';
export { DefaultTileBins, dataKey, fitLevel, nodeKey, tileRegion, tileSpan, tileWidth } from './keys.js';
export {
  type AssetReader,
  abortError,
  type Coverage,
  type DynamicNode,
  type FailureKind,
  type FieldInfo,
  type FormatRuntime,
  isAbortError,
  type LettersTable,
  type OpenContext,
  type OpenedSource,
  type RangesDensity,
  type RangesRows,
  type RangesTable,
  type RegionAnswer,
  type RegionRequest,
  type Shape,
  type ShapeInfo,
  Shapes,
  type SourceEntry,
  type SystemInfo,
  type Table,
  type Tables,
  type ValuesBins,
  type ValuesPositions,
  type ValuesTable,
} from './protocol.js';
export {
  decodeSeq,
  encodeSeq,
  inlineSeq,
  readSeq,
  type SeqData,
  SeqKind,
  type SeqPair,
  seqToNode,
  seqToValue,
} from './seq-format.js';
export { openSeq, SeqRuntime, SourceOpenError } from './seq-runtime.js';
export {
  type DataEntry,
  DataStore,
  DefaultStoreSettings,
  type EntryStatus,
  type NodeEntry,
  type Priority,
  RealClock,
  type StoreClock,
  type StoreSettings,
  type TileRegion,
  type Want,
} from './store.js';
export { answerBytes, checkTables, isEmptyAnswer, NodeBytes, sliceTable, tableRows } from './tables.js';
