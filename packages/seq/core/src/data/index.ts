/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seqx-core/data`: the source protocol, the `seq` format and its runtime, and the data store
// (_spec/README.md "Data"). Named re-exports only.

export { asyncSeqSource, type DelayOptions, delayedRuntime, delayedSource } from './async-source.js';
export { checkFai, FastaRuntime, type FaiRecord, faidx, openFasta, parseFai } from './fasta.js';
export { SourceIssues, StoreIssues } from './issues.js';
export { DefaultTileBins, dataKey, fitLevel, nodeKey, tileRegion, tileSpan, tileWidth } from './keys.js';
export {
  type AssetReader,
  type CellsCodes,
  type CellsMajority,
  type CellsTable,
  abortError,
  type Coverage,
  type DynamicNode,
  type FailureKind,
  FailureKinds,
  type FieldInfo,
  type FormatRuntime,
  isAbortError,
  type LettersTable,
  DefaultOpenLimits,
  type OpenContext,
  type OpenLimits,
  type OpenedSource,
  type RangesDensity,
  type RangesField,
  type RangesList,
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
  checkCoverage,
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
  type Wanter,
} from './store.js';
export {
  binValues,
  coarsenBins,
  DefaultSummarySettings,
  type SummarySettings,
  withSummaries,
} from './summaries.js';
export {
  answerBytes,
  type Cut,
  type AnswerBounds,
  checkAnswer,
  checkTableBounds,
  checkTables,
  DefaultAnswerBounds,
  DefaultMaxBlockValues,
  DefaultMaxColumnLength,
  DefaultMaxRows,
  isEmptyAnswer,
  joinTables,
  NodeBytes,
  sliceTable,
  tableRows,
} from './tables.js';
export {
  AlignmentRuntime,
  type AlignmentRows,
  AlignmentSystems,
  cellsSource,
  openAlignment,
} from './alignment-source.js';
