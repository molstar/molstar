/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The default composition (seqspec:viewer#the-viewer-spec): every built-in registry entry and behavior, in the shape
// of Mol* v6's DefaultPluginSpec. Nothing in the viewer imports it; a host that wants less lists its own parts.

import {
  KernelRegistry,
  list,
  MvsRegistry,
  OptionalField,
  ref,
  SimpleParamsSchema,
  type SvsRegistryEntry,
} from '@molstar/seq-core/format';
import { DefaultActions } from './actions.js';
import { CommandsBehavior, SeqCommands } from './commands.js';
import {
  RevealHighlightBehavior,
  CursorBehavior,
  FocusViewportBehavior,
  HeaderHoverBehavior,
  HoverBehavior,
  HoverLabelBehavior,
  HoverRowBehavior,
  ScrollGroupsBehavior,
  SnapshotControlsBehavior,
} from './default-behaviors.js';
import { DotPlotSelectBehavior } from './dot-plot.js';
import { FollowBehavior } from './follow.js';
import {
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
import { KeyboardBehavior } from './keyboard.js';
import { DefaultLabelProviders } from './label-providers.js';
import { BookmarkActions } from './location-sets.js';
import { DefaultNavigators } from './navigators.js';
import { DefaultNodeKinds } from './node-kinds.js';
import { DefaultPresets } from './presets.js';
import type { SeqViewerRegistryEntry } from './registry.js';
import { KernelColumns, KernelRenderers } from './renderer-list.js';
import { RingWindowBehavior } from './ring-window.js';
import { DefaultThemes } from './rules.js';

import { SeqViewerSpec, ViewerIssueCodes, viewerOwn } from './spec.js';

export const ViewerLayerTypes: SvsRegistryEntry = viewerOwn({
  layerTypes: [
    {
      name: 'seqv:map',
      description: 'other tracks in miniature over the whole axis, the view as a window',
      options: SimpleParamsSchema({ tracks: OptionalField(list(ref('track')), [], 'the tracks; empty: all others') }),
      accepts: 'none',
      layouts: ['tracks'],
      axes: ['blocks', 'columns'],
      sizing: 'free',
      fallback: 'omit',
      channels: [],
    },
  ],
});

/**
 * Every built-in, both halves; the phases that build renderers, themes and the rest add theirs. The viewer owns `seqv`:
 * its registry passes `ViewerIssueCodes` and `ViewerLayerTypes` to the format's `register` as its own (registry.ts).
 */
// Built in a pure call, so a host that does not register it (a partial spec) bundles none of it.
export const DefaultRegistry: SeqViewerRegistryEntry = /* @__PURE__ */ (() => ({
  svs: [KernelRegistry, MvsRegistry, ViewerIssueCodes, ViewerLayerTypes],
  nodeKinds: DefaultNodeKinds,
  navigators: DefaultNavigators,
  renderers: KernelRenderers,
  columns: KernelColumns,
  themes: DefaultThemes,
  presets: DefaultPresets,
  labels: DefaultLabelProviders,
  commands: Object.values(SeqCommands),
  actions: [...DefaultActions, ...BookmarkActions],
}))();

/** The default behaviors, in `DefaultSeqViewerSpec` order (seqspec:viewer#behaviors). */
export const SeqViewerBehaviors = {
  Scrollbar: ScrollbarBehavior,
  Navigator: NavigatorBehavior,
  RingWindow: RingWindowBehavior,
  VerticalOverview: VerticalOverviewBehavior,
  Headers: HeadersBehavior,
  Select: SelectBehavior,
  DotPlotSelect: DotPlotSelectBehavior,
  DragPan: DragPanBehavior,
  Wheel: WheelBehavior,
  DoubleClickZoom: DoubleClickZoomBehavior,
  RowZoom: RowZoomBehavior,
  Keyboard: KeyboardBehavior,
  Hover: HoverBehavior,
  HeaderHover: HeaderHoverBehavior,
  HoverRow: HoverRowBehavior,
  HoverLabel: HoverLabelBehavior,
  Cursor: CursorBehavior,
  FocusViewport: FocusViewportBehavior,
  ScrollGroups: ScrollGroupsBehavior,
  Follow: FollowBehavior,
  SnapshotControls: SnapshotControlsBehavior,
  Commands: CommandsBehavior,
};

/**
 * Built-in behaviors not in `DefaultSeqViewerSpec`: a host adds them (`SeqViewerSpec.Behavior(…)`). `RevealHighlight`:
 * a linked view's highlight, rested, comes into view (seqspec:interaction#highlight).
 */
export const SeqViewerOptionalBehaviors = {
  RevealHighlight: RevealHighlightBehavior,
};

/** A new object each call. */
export function DefaultSeqViewerSpec(): SeqViewerSpec {
  return {
    registry: [DefaultRegistry],
    behaviors: Object.values(SeqViewerBehaviors).map((b) => SeqViewerSpec.Behavior(b)),
  };
}
