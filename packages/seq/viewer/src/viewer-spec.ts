/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ViewLayout } from '@molstar/seqx-core/render';
import type { SeqConfigItem, SeqPorts, SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { BindingTable, GestureBinding } from '@molstar/seqx-viewer/bindings';
import type { ViewerConfig } from '@molstar/seqx-viewer/config';
import type { SeqAnimationPort, SeqHistoryPort } from '@molstar/seqx-viewer/ports';
import type { LegendModel } from './legends.js';
import type { SeqViewer } from './viewer.js';
import type {
  AddressModel,
  HighlightModel,
  HoverModel,
  SelectionModel,
  SpecModel,
  StoreModel,
  ViewModel,
} from '@molstar/seqx-viewer/models';

// The viewer's spec (_spec/README.md "Composition"): what a viewer contains, listed in one place, as a Mol* plugin's
// spec does. Nothing registers on import; `DefaultSeqViewerSpec()` lists the built-in behaviors and ports, and an
// application spreads it and changes what it needs.

/** What a behavior is given: the models, the bindings and the configuration of the viewer it is part of. */
export interface BehaviorContext {
  readonly spec: SpecModel;
  readonly view: ViewModel;
  readonly hover: HoverModel;
  readonly store: StoreModel;
  readonly selection: SelectionModel;
  readonly highlight: HighlightModel;
  /** Addresses resolved in the spec shown, the actions on them, the refs, links, the session as addresses. */
  readonly addresses: AddressModel;
  readonly config: ViewerConfig;
  readonly bindings: BindingTable;
  /** The legends as derived, the hovered entry, the reader's color-blind substitution. */
  readonly legends: LegendModel;
  /** The last frame's layout (a rows band's state: its clades, its rows in view), undefined before the first. */
  layout(): ViewLayout | undefined;
}

/** A behavior's provider: creates the policy over the models; what it makes is its handle (`dispose` ends it). */
export interface BehaviorProvider<Params = unknown, Handle = unknown> {
  readonly name: string;
  create(context: BehaviorContext, params: Params): Handle & { dispose?(): void };
}

export interface SeqViewerBehavior<Params = any> {
  readonly provider: BehaviorProvider<Params, any>;
  readonly params?: Params;
}

export interface SeqViewerBindings {
  /** The preset (`default` when absent, or the `bindings.preset` setting). */
  readonly preset?: readonly GestureBinding[];
  /** Layers over the preset, each replacing the bindings of the same gesture, modifiers and target. */
  readonly overrides?: readonly (readonly GestureBinding[])[];
  /** The ids of commands bindings may name besides the viewer's own (an extension's). */
  readonly actions?: readonly string[];
}

/** The mounted UI's parts: the toolbar and the canvas stack (the section labels are drawn by the layout). */
/** Mounts a toolbar for a viewer into a target (before the canvas); the function it returns takes it away. */
export type ToolbarMount = (viewer: SeqViewer, target: HTMLElement) => () => void;

export interface SeqViewerUI {
  readonly canvas?: true;
  /**
   * The toolbar, a part the host places (like the legends): `mountToolbar` of `@molstar/seqx-viewer/mount-toolbar`, or
   * its own, mounted into the target above the canvas. Not listed, not in the bundle: the toolbar is `@molstar/ui`'s
   * and weighs far more than the two buttons it holds.
   */
  readonly toolbar?: ToolbarMount;
  /**
   * A mounted viewer with no `history` given gets a `@molstar/ui/history` stack of its own (`uiHistory()`): the
   * standalone viewer's Undo. A headless viewer never does.
   */
  readonly history?: true;
  /**
   * The legends, a part the host places (`mountLegend` of `@molstar/seqx-viewer/mount-legend`, or its own): mounted
   * after the canvas in the standalone viewer. Seq* never draws a legend on its canvas.
   */
  readonly legends?: ToolbarMount;
  /**
   * The navigators, parts the host places (`navigatorsMount` of `@molstar/seqx-viewer/mount-navigator` places the
   * view's declared navigators by their `place` around the canvas through `mountNavigator`; a host calls
   * `mountNavigator` itself). On in the standalone viewer.
   */
  readonly navigators?: ToolbarMount;
  /**
   * A button in the canvas's corner that hides and shows the ruler rails (`railToggleMount`): the reader's session,
   * never the spec. For a viewer embedded in a host's cell, where its room is the reader's.
   */
  readonly railToggle?: ToolbarMount;
}

export interface SeqViewerSpec {
  /** Registry entries, in order: what the viewer can read, compute and draw. */
  readonly registry?: readonly SeqRegistryEntry[];
  /** Behaviors with their params, initialized after the registry. */
  readonly behaviors: readonly SeqViewerBehavior[];
  readonly bindings?: SeqViewerBindings;
  /** Settings as `[item, value]` pairs. */
  readonly config?: readonly (readonly [SeqConfigItem<any>, unknown])[];
  /** Port adapters, each a factory with ownership; unset ports take their default adapter. */
  readonly ports?: SeqPorts;
  /**
   * The host's history (Seq*'s history port, _spec/README.md "Ports"): the spec's edits as steps of the host's stack.
   * None: edits apply at once and nothing undoes them (a headless viewer, a test).
   */
  readonly history?: SeqHistoryPort;
  /**
   * The application's clock (Seq*'s animation port, _spec/README.md "Ports", Animation): Seq*'s own motions (a jump, a
   * story step) are plans it runs. None: the page's frames (`frameAnimation`); headless, a fake clock
   * (`headlessAnimation`).
   */
  readonly animation?: SeqAnimationPort;
  /** The UI's parts; none for a headless viewer (a target given without it is left empty). */
  readonly ui?: SeqViewerUI;
}
