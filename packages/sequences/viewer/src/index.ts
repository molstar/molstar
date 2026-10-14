/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The headless root entry: the host contract, every name listed (seqspec:viewer#the-viewer-the-host-contract). A host
// makes a viewer from a spec (`createSeqViewer`, `SeqViewerSpec`, the registry entry types and the defaults it
// registers), reads and drives its context (the managers' interfaces, with the states, events and records they hand
// out), hands it input records, adds behaviors (the providers and `defineBehavior`), sends commands and actions, and
// sets config. What a renderer uses is `/draw`'s (draw/api.ts), a test's helpers `/testing`'s; the default UIs' model
// (chrome.ts, search-box.ts, hover-card.ts, header-text.ts) is internal to `/dom` and `/react`, and the managers'
// implementations to the context. Tests reach internals by source path.
export {
  type ActionEntry,
  type ActionManager,
  type ColorRangeTarget,
  type ColorSchemeTarget,
  DefaultActions,
  type LegendTarget,
  nodePath,
  trackDataRefs,
} from './actions.js';
export type { AxisFrame, SectionAxisView } from './axis.js';
export {
  type BehaviorHandle,
  type BehaviorInfo,
  type BehaviorManager,
  type SeqBehavior,
  type SeqBehaviorProvider,
  SeqBehaviors,
} from './behaviors.js';
export { Channel, ChannelMarks, type ChannelSources } from './channels.js';
// The default UIs' compositions and what a host's own composition reads (seqspec:interaction#slots-and-chrome).
export {
  ChromeCompositions,
  CleanParts,
  descriptionOf,
  type SectionInfo,
  type SectionSwitch,
  sectionInfo,
  shownSections,
} from './chrome.js';
export {
  type CommandManager,
  CommandsBehavior,
  SeqCommand,
  SeqCommands,
  type SeqTarget,
  TargetRangeError,
} from './commands.js';
export {
  type ConfigManager,
  type CustomTrackAction,
  type LegendMode,
  type PresetDefaults,
  SeqViewerConfig,
  SeqViewerConfigItem,
  type TrackActionId,
  type TrackActionList,
} from './config.js';
export {
  type ActivateEvent,
  createSeqViewer,
  DarkSkin,
  DefaultSkin,
  type GraphResult,
  type LogEntry,
  type PatchEvent,
  type SeqIssue,
  type SeqViewerContext,
  type SeqViewerEvents,
  type SeqViewerPorts,
} from './context.js';
export type { SeqDataInput, SetDataOptions } from './data.js';
export {
  type BehaviorKit,
  CursorBehavior,
  defineBehavior,
  ContextMenuBehavior,
  FocusViewportBehavior,
  HeaderHoverBehavior,
  HoverBehavior,
  HoverLabelBehavior,
  HoverRowBehavior,
  ScrollGroupsBehavior,
  SnapshotControlsBehavior,
} from './default-behaviors.js';
export type {
  DocumentManager,
  DocumentPort,
  HistoryManager,
  LocalPatch,
  ProposeResult,
  RouteContext,
  Routed,
  ViewInfo,
  ViewManager,
} from './document.js';
// The ideogram navigator a host registers (`navigators`), apart from the renderer API (VL-16).
export { type CytoBands, ideogramNavigator } from './draw/ideogram.js';
export type { DynamicSource, HeldNode } from './dynamic.js';
export type { FocusChange, FocusEntry, FocusManager } from './focus.js';
export { FollowBehavior } from './follow.js';
export {
  DoubleClickZoomBehavior,
  DragPanBehavior,
  HeadersBehavior,
  NavigatorBehavior,
  RowZoomBehavior,
  ScrollbarBehavior,
  SelectBehavior,
  VerticalOverviewBehavior,
  WheelBehavior,
} from './gestures.js';
export type {
  Assembled,
  DataEntry,
  DataMode,
  InputCollector,
  NodeKind,
  NodeStatus,
  PreparedLayer,
  RunContext,
  SeqGraph,
  Skin,
  Stage,
  ViewModel,
} from './graph.js';
export {
  DefaultHeaderHoverParams,
  type HeaderHoverParams,
  type HeaderLight,
  type HeaderSet,
  type LociSize,
  lociSize,
} from './header-hover.js';
export type { HighlightEntry, HighlightManager, HighlightOptions, HighlightStyle } from './highlight.js';
export { Hit, HitRegion, HoverKind } from './hit.js';
export type { HoverManager, HoverState, PointerState } from './hover.js';
export type { ItemTarget } from './hover-card.js';
export {
  Binding,
  Button,
  type ClickEvt,
  DragEvt,
  DragPhase,
  type InputSink,
  Inputs,
  type InteractionEvents,
  KeyBinding,
  type KeyEvt,
  type KeyInput,
  type KeyTrigger,
  Modifier,
  PointerEvt,
  type PointerInput,
  PointerType,
  type Trigger,
  type WheelEvt,
  type WheelInput,
} from './input.js';
export { KeyboardBehavior } from './keyboard.js';
export { DefaultLabelProviders } from './label-providers.js';
export type {
  HoverLabel,
  LabelContext,
  LabelManager,
  LabelProvider,
  LabelSink,
  ReachedNumbering,
} from './labels.js';
export {
  DefaultRenderProps,
  EntryKind,
  type LayoutManager,
  type RenderProps,
  type RowModel,
  type SectionLayout,
  type VisibleRows,
} from './layout.js';
export type { HeaderLegend, LegendItem, LegendRamp, LegendSize } from './legend.js';
export type { SeqLociAny } from './lineage.js';
export {
  BookmarkActions,
  type BookmarkAdd,
  type BookmarkItem,
  type BookmarkStep,
  bookmarkAt,
  bookmarksOf,
  type HighlightRemove,
  type ResultsSave,
} from './location-sets.js';
export { MarkBit, type SectionMarks, type ShownRuns } from './marks.js';
export {
  DefaultNavigatorParts,
  DefaultNavigators,
  type NavigatorAxis,
  type NavigatorEntry,
  type NavigatorParts,
  navigatorFrom,
} from './navigators.js';
export { DefaultNodeKinds } from './node-kinds.js';
export { DefaultPresets, type PresetEntry } from './presets.js';
export { type Coordinates, isCoordinates } from './projection.js';
export type { RangeArgs, RefManager } from './refs.js';
export {
  type ColumnRenderer,
  type LayerInput,
  type LayerRenderer,
  type NamedEntry,
  type PrepareContext,
  SeqRegistrationError,
  type SeqViewerRegistry,
  type SeqViewerRegistryEntry,
  type ViewerEntries,
  type ViewerEntryOf,
  type ViewerKind,
} from './registry.js';
export { KernelColumns, KernelRenderers } from './renderer-list.js';
export type { Result, ResultsManager, ResultsOptions, ResultsState } from './results.js';
export { DefaultThemes } from './rules.js';
export { type FrameHost, type ManualScheduler, manualScheduler, rafScheduler, type Scheduler } from './scheduler.js';
export type { SearchFn } from './search-box.js';
export type { PendingSelection, Selection, SelectionManager } from './selection.js';
export { showExtents } from './show.js';
export type { SnapshotInfo, SnapshotManager } from './snapshots.js';
export {
  DefaultRegistry,
  DefaultSeqViewerSpec,
  SeqViewerBehaviors,
  SeqViewerSpec,
  ViewerIssueCodes,
  ViewerLayerTypes,
} from './spec.js';
export type {
  Live,
  Origin,
  SeqEvent,
  SeqEventEmitter,
  SeqState,
  Unsubscribable,
  WritableLive,
  WritableSeqState,
} from './state.js';
// The stream manager as a host reads it; what drawing reads is `/draw`'s, the rest internal.
export { RegionState, type StreamManager, type StreamStats } from './streams.js';
export { SystemRegistry } from './systems.js';
export { monospaceMeasurer, type TextManagerApi, type TextMeasurer } from './text.js';
export {
  KernelThemes,
  type Theme,
  type ThemeContext,
  ThemeLocation,
  type ThemeParams,
  type ThemeParts,
  type ThemeProvider,
  type ThemeValueFn,
  themeProvider,
} from './themes.js';
export type {
  DragState,
  IssueManager,
  KeyCursor,
  RowHoverOf,
  RowHoverParams,
  TaskManager,
  UiManager,
  WindowHover,
} from './ui.js';
export {
  type CommittedViewport,
  Lod,
  type MoveOpts,
  type ScrollGroupManager,
  type SectionViewport,
  type ViewportLive,
  type ViewportManager,
  type ViewportProps,
  type ViewportState,
} from './viewport.js';
