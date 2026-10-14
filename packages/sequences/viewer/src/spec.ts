/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer spec (seqspec:viewer#the-viewer-spec), in the shape of Mol* v6's PluginSpec: registry entries, behaviors,
// config.

import {
  KernelRegistry,
  list,
  MvsRegistry,
  OptionalField,
  ref,
  SimpleParamsSchema,
  type SvsRegistryEntry,
} from '@molstar/seq-core/format';
import type { Severity } from '@molstar/seq-core/model';
import { DefaultActions } from './actions.js';
import type { SeqBehaviorProvider } from './behaviors.js';
import { CommandsBehavior, SeqCommands } from './commands.js';
import type { SeqViewerConfigItem } from './config.js';
import {
  ContextMenuBehavior,
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
import type { RenderProps } from './layout.js';
import { BookmarkActions } from './location-sets.js';
import { DefaultNavigators } from './navigators.js';
import { DefaultNodeKinds } from './node-kinds.js';
import { DefaultPresets } from './presets.js';
import type { SeqViewerRegistryEntry } from './registry.js';
import { KernelColumns, KernelRenderers } from './renderer-list.js';
import { RingWindowBehavior } from './ring-window.js';
import { DefaultThemes } from './rules.js';

export interface SeqViewerSpec {
  /** Registered by createSeqViewer, in order, before any behavior is created. */
  readonly registry?: readonly SeqViewerRegistryEntry[];
  readonly behaviors: readonly SeqViewerSpec.Behavior[];
  readonly config?: readonly (readonly [SeqViewerConfigItem, unknown])[];
  /** The analog of Mol*'s canvas3d props: the initial `ctx.renderProps` over the defaults. */
  readonly render?: Partial<RenderProps>;
}
export declare namespace SeqViewerSpec {
  interface Behavior {
    readonly behavior: SeqBehaviorProvider<any>;
    readonly defaultParams?: object;
  }
}
export const SeqViewerSpec = {
  Behavior<P extends object>(behavior: SeqBehaviorProvider<P>, defaultParams?: Partial<P>): SeqViewerSpec.Behavior {
    return { behavior, defaultParams };
  },
};

/**
 * The brand of the viewer's own format entries, the only ones `register` lets name `seqv`. A registered symbol, not the
 * objects' identity: a page holding two copies of the viewer (an HMR update, a host that bundles it twice) makes its
 * spec with one and its context with the other. Not enumerable, so a spread copy (an extension's) is not the viewer's.
 */
const OWN = Symbol.for('seqstar.viewer.seqv');
const own = (e: SvsRegistryEntry): SvsRegistryEntry => Object.defineProperty(e, OWN, { value: true });
/** Whether `e` is one of the viewer's own format entries, from any copy of the viewer. */
export const isViewerOwn = (e: SvsRegistryEntry): boolean => (e as { [OWN]?: boolean })[OWN] === true;

/**
 * The viewer's own issue codes (seqspec:viewer#registration), declared in the format's half so validation knows them.
 */
const codes = (severity: Severity, names: string) => names.split(' ').map((c) => ({ code: `seqv.${c}`, severity }));
export const ViewerIssueCodes: SvsRegistryEntry = own({
  issueCodes: [
    ...codes('warning', 'patch.dropped patch.skipped history.refused port.refused issues.truncated ui.slot'),
    // Ring labels that do not fit are hidden: that depends on the width, not the input (seqspec:format#validation).
    ...codes('info', 'render.labels'),
    ...codes('error', 'renderer.error behavior.error command.error state.error target.range'),
  ],
});

/** The viewer's own layer type (seqspec:viewer#extension-points): other tracks' overviews over the whole axis. */
export const ViewerLayerTypes: SvsRegistryEntry = own({
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
export const DefaultRegistry: SeqViewerRegistryEntry = {
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
};

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
 * Built-in behaviors not in `DefaultSeqViewerSpec`: a host adds them (`SeqViewerSpec.Behavior(…)`). `ContextMenu`: a
 * secondary press asks for a track's menu at the pointer (seqspec:interaction#the-track-menu).
 */
export const SeqViewerOptionalBehaviors = {
  ContextMenu: ContextMenuBehavior,
};

/** A new object each call. */
export function DefaultSeqViewerSpec(): SeqViewerSpec {
  return {
    registry: [DefaultRegistry],
    behaviors: Object.values(SeqViewerBehaviors).map((b) => SeqViewerSpec.Behavior(b)),
  };
}
