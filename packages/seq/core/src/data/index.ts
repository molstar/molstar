/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seqx-core/data`: the source protocol, the `seq` format and its runtime, and the data store
// (_spec/README.md "Data"). Named re-exports only.

export {
  type Archive,
  ArchiveEntryMissing,
  ArchiveIndexName,
  type ArchiveLimits,
  type ArchiveRead,
  ArchiveSpecName,
  type ArchiveWritten,
  archiveReader,
  DefaultArchiveLimits,
  readArchive,
  writeArchive,
} from './archive.js';
export { asyncSeqSource, type DelayOptions, delayedRuntime, delayedSource } from './async-source.js';
export { checkFai, FastaRuntime, type FaiRecord, faidx, openFasta, parseFai } from './fasta.js';
export { ArchiveIssues, SourceIssues, StoreIssues } from './issues.js';
export { DefaultMaxNameHits, findNames, type NamedSource, type SourceNameHit } from './names.js';
export {
  DefaultGridTileBins,
  DefaultTileBins,
  dataKey,
  fitLevel,
  foldTile,
  type KeyParts,
  nodeKey,
  tileRegion,
  tileSpan,
  tileWidth,
} from './keys.js';
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
  type FindNamesOptions,
  type FormatRuntime,
  type GridDense,
  type GridSparse,
  type GridTable,
  isAbortError,
  ItemShapes,
  type LettersTable,
  type MatrixSummary,
  type MatrixItems,
  type MatrixTable,
  type NameHit,
  DefaultOpenLimits,
  type OpenContext,
  type OpenLimits,
  type OpenedSource,
  type PairsBins,
  type PairsList,
  type PairsTable,
  type RangesDensity,
  type RangesField,
  type RangesList,
  type RangesTable,
  type RegionAnswer,
  type RegionRequest,
  type RowsKeyed,
  type RowsTable,
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
  type WholeSystem,
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
export { type GridPlane, gridPlane, transposeGrid } from './summaries-grid.js';
export { type BinPairsOptions, binPairs, coarsenPairBins } from './summaries-pairs.js';
export { FrameSlices, type SliceClock } from './summaries-region.js';
export type { LayerLoad, WantTile } from './store-wants.js';
export { sidedName } from './store-wants.js';
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
