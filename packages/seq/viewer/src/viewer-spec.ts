/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqConfigItem, SeqPorts, SeqRegistryEntry } from '@molstar/seqx-core/context';
import type { BindingTable, GestureBinding } from '@molstar/seqx-viewer/bindings';
import type { ViewerConfig } from '@molstar/seqx-viewer/config';
import type { HoverModel, SpecModel, StoreModel, ViewModel } from '@molstar/seqx-viewer/models';

// The viewer's spec (_spec/README.md "Composition"): what a viewer contains, listed in one place, as a Mol* plugin's
// spec does. Nothing registers on import; `DefaultSeqViewerSpec()` lists the built-in behaviors and ports, and an
// application spreads it and changes what it needs.

/** What a behavior is given: the models, the bindings and the configuration of the viewer it is part of. */
export interface BehaviorContext {
  readonly spec: SpecModel;
  readonly view: ViewModel;
  readonly hover: HoverModel;
  readonly store: StoreModel;
  readonly config: ViewerConfig;
  readonly bindings: BindingTable;
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

/** The mounted UI's parts: the canvas stack in this version (the section labels are drawn by the layout). */
export interface SeqViewerUI {
  readonly canvas?: true;
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
  /** The UI's parts; none for a headless viewer (a target given without it is left empty). */
  readonly ui?: SeqViewerUI;
}
