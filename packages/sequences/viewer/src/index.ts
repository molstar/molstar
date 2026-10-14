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
} from "./actions.ts";
export type { AxisFrame, SectionAxisView } from "./axis.ts";
export {
  type BehaviorHandle,
  type BehaviorInfo,
  type BehaviorManager,
  type SeqBehavior,
  type SeqBehaviorProvider,
  SeqBehaviors,
} from "./behaviors.ts";
export { Channel, ChannelMarks, type ChannelSources } from "./channels.ts";
// The default UIs' compositions and what a host's own composition reads (seqspec:interaction#slots-and-chrome).
export {
  ChromeCompositions,
  CleanParts,
  descriptionOf,
  type SectionInfo,
  type SectionSwitch,
  sectionInfo,
  shownSections,
} from "./chrome.ts";
export {
  type CommandManager,
  CommandsBehavior,
  SeqCommand,
  SeqCommands,
  type SeqTarget,
  TargetRangeError,
} from "./commands.ts";
export {
  type ConfigManager,
  type CustomTrackAction,
  type LegendMode,
  type PresetDefaults,
  SeqViewerConfig,
  SeqViewerConfigItem,
  type TrackActionId,
  type TrackActionList,
} from "./config.ts";
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
} from "./context.ts";
export type { SeqDataInput, SetDataOptions } from "./data.ts";
export {
  type BehaviorKit,
  CursorBehavior,
  defineBehavior,
  FocusViewportBehavior,
  HeaderHoverBehavior,
  HoverBehavior,
  HoverLabelBehavior,
  HoverRowBehavior,
  ScrollGroupsBehavior,
  SnapshotControlsBehavior,
} from "./default-behaviors.ts";
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
} from "./document.ts";
// The ideogram navigator a host registers (`navigators`), apart from the renderer API (VL-16).
export { type CytoBands, ideogramNavigator } from "./draw/ideogram.ts";
export type { DynamicSource, HeldNode } from "./dynamic.ts";
export type { FocusChange, FocusEntry, FocusManager } from "./focus.ts";
export { FollowBehavior } from "./follow.ts";
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
} from "./gestures.ts";
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
} from "./graph.ts";
export {
  DefaultHeaderHoverParams,
  type HeaderHoverParams,
  type HeaderLight,
  type HeaderSet,
  type LociSize,
  lociSize,
} from "./header-hover.ts";
export type { HighlightEntry, HighlightManager, HighlightOptions, HighlightStyle } from "./highlight.ts";
export { Hit, HitRegion, HoverKind } from "./hit.ts";
export type { HoverManager, HoverState, PointerState } from "./hover.ts";
export type { ItemTarget } from "./hover-card.ts";
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
} from "./input.ts";
export { KeyboardBehavior } from "./keyboard.ts";
export { DefaultLabelProviders } from "./label-providers.ts";
export type {
  HoverLabel,
  LabelContext,
  LabelManager,
  LabelProvider,
  LabelSink,
  ReachedNumbering,
} from "./labels.ts";
export {
  DefaultRenderProps,
  EntryKind,
  type LayoutManager,
  type RenderProps,
  type RowModel,
  type SectionLayout,
  type VisibleRows,
} from "./layout.ts";
export type { HeaderLegend, LegendItem, LegendRamp, LegendSize } from "./legend.ts";
export type { SeqLociAny } from "./lineage.ts";
export {
  BookmarkActions,
  type BookmarkAdd,
  type BookmarkItem,
  type BookmarkStep,
  bookmarkAt,
  bookmarksOf,
  type HighlightRemove,
  type ResultsSave,
} from "./location-sets.ts";
export { MarkBit, type SectionMarks, type ShownRuns } from "./marks.ts";
export {
  DefaultNavigatorParts,
  DefaultNavigators,
  type NavigatorAxis,
  type NavigatorEntry,
  type NavigatorParts,
  navigatorFrom,
} from "./navigators.ts";
export { DefaultNodeKinds } from "./node-kinds.ts";
export { DefaultPresets, type PresetEntry } from "./presets.ts";
export { type Coordinates, isCoordinates } from "./projection.ts";
export type { RangeArgs, RefManager } from "./refs.ts";
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
} from "./registry.ts";
export { KernelColumns, KernelRenderers } from "./renderer-list.ts";
export type { Result, ResultsManager, ResultsOptions, ResultsState } from "./results.ts";
export { DefaultThemes } from "./rules.ts";
export { type FrameHost, type ManualScheduler, manualScheduler, rafScheduler, type Scheduler } from "./scheduler.ts";
export type { SearchFn } from "./search-box.ts";
export type { PendingSelection, Selection, SelectionManager } from "./selection.ts";
export { showExtents } from "./show.ts";
export type { SnapshotInfo, SnapshotManager } from "./snapshots.ts";
export {
  DefaultRegistry,
  DefaultSeqViewerSpec,
  SeqViewerBehaviors,
  SeqViewerSpec,
  ViewerIssueCodes,
  ViewerLayerTypes,
} from "./spec.ts";
export type {
  Live,
  Origin,
  SeqEvent,
  SeqEventEmitter,
  SeqState,
  Unsubscribable,
  WritableLive,
  WritableSeqState,
} from "./state.ts";
// The stream manager as a host reads it; what drawing reads is `/draw`'s, the rest internal.
export { RegionState, type StreamManager, type StreamStats } from "./streams.ts";
export { SystemRegistry } from "./systems.ts";
export { monospaceMeasurer, type TextManagerApi, type TextMeasurer } from "./text.ts";
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
} from "./themes.ts";
export type {
  DragState,
  IssueManager,
  KeyCursor,
  RowHoverOf,
  RowHoverParams,
  TaskManager,
  UiManager,
  WindowHover,
} from "./ui.ts";
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
} from "./viewport.ts";
