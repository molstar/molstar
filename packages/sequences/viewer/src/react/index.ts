/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `@molstar/seq-viewer/react`: the default UI over the root and canvas entries. The provider and hooks, every part (each
// also a slot, and a component a host places anywhere under the provider; `Issues` is exported as `IssueBadge`), the
// default primitives, the root and the default composition, `SeqViewer`. Its styles are
// `@molstar/seq-viewer/skin.css`, shared with `/dom`. The React and DOM entries implement one slot contract and share
// its names (`ViewTabs`, `PartProps`, …): a host imports one of them (the export guard allows exactly that overlap).

export type { MenuItem } from '../chrome.js';
export { withDownloadExport } from '../dom/export.js';
export type { AnnotationPlace } from '../annotations.js';
export { AnnotationsPanel, type AnnotationsPanelProps } from './annotations.js';
export {
  type DerivedOptions,
  type HintProps,
  type MenuProps,
  type MountedSurface,
  type PartProps,
  type ResolvedUISpec,
  type RowProps,
  type ScrollAreaProps,
  type SectionProps,
  type SeqViewerUISpec,
  type SlotComponents,
  type SlotName,
  type SurfaceFactory,
  type TrackProps,
  useDerived,
  useLanding,
  useLiveVersion,
  useSeqState,
  useSeqViewer,
} from './hooks.js';
export { HoverCard } from './hover-card.js';
export {
  Busy,
  Description,
  EmptyState,
  ErrorState,
  Export,
  HiddenChips,
  Hint,
  IssueBadge,
  Legend,
  Menu,
  ScrollArea,
  Search,
  type SearchProps,
  SectionHeader,
  SectionSwitches,
  SnapshotPlayer,
  Tooltip,
  TrackActions,
  ViewTabs,
} from './parts.js';
export { ViewSurface } from './section.js';
export { SeqViewer, SeqViewerLayout, type SeqViewerProps, SeqViewerProvider, SeqViewerRoot } from './viewer.js';
