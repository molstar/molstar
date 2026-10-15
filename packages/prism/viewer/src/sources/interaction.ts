/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { merge } from 'rxjs';
import { selectionSummary } from './selection.js';
import { fromBehavior, fromChanged, type SourceScope } from './state.js';

// The selection mode, the picking level and the selection as states (state.md §3: they are the plugin's state; what was
// missing is their source): the mode switch, the picking level, the selection bar and the agent's account read them.

export interface InteractionSources {
  /** Clicks select (true) or focus (false). */
  readonly selectionMode: ReadonlyState<boolean>;
  /** What a click picks: Mol*'s granularity (`residue`, `chain`, …). */
  readonly granularity: ReadonlyState<string>;
  /** What is selected, in Mol*'s words ("12 residues"); undefined: nothing. */
  readonly selection: ReadonlyState<string | undefined>;
  /** The focused part, in Mol*'s words; undefined: none. */
  readonly focus: ReadonlyState<string | undefined>;
}

export function createInteractionSources(scope: SourceScope, plugin: PluginContext): InteractionSources {
  const selectionMode = fromBehavior(scope, plugin.behaviors.interaction.selectionMode);
  const granularity = fromChanged(
    scope,
    () => String(plugin.managers.interactivity.props.granularity),
    plugin.managers.interactivity.events.propsUpdated,
  );
  const selected = () =>
    [...plugin.managers.structure.selection.entries.values()].map((e) => e.selection as StructureElement.Loci);
  const selection = fromChanged(
    scope,
    () => selectionSummary(selected()),
    merge(plugin.managers.structure.selection.events.changed, plugin.state.data.events.changed),
  );
  const focus = fromChanged(
    scope,
    () => {
      const loci = plugin.managers.structure.focus.current?.loci;
      return loci ? selectionSummary([loci]) : undefined;
    },
    plugin.managers.structure.focus.behaviors.current,
  );
  return { selectionMode, granularity, selection, focus };
}
