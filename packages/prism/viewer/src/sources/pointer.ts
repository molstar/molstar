/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { mapState, type ReadonlyState } from '@molstar/prism-ui/core';
import { map, type Observable, share } from 'rxjs';
import { summarize } from '../targets/index.js';
import { hoverLines } from './hover-lines.js';
import { structureNames } from './naming.js';
import { fromBehavior, type SourceScope } from './state.js';

// What is under the pointer (api.md §8: v5's `structureInteractivity`; README §3.6): the hover label's lines, and the
// hover and click events embedding hosts read, each with the target in the address language (targets.md §1).

export interface ViewerPointerEvent {
  /** The plugin's label of it, as lines; empty: nothing. */
  readonly lines: readonly string[];
  /** What it is, as a target ("A:120", "A:120@CA"); undefined: not part of a structure. */
  readonly target?: string;
  /** The structure it is in, by name. */
  readonly structure?: string;
  /** The mouse buttons and modifier keys. */
  readonly button: number;
  readonly modifiers: {
    readonly shift: boolean;
    readonly control: boolean;
    readonly alt: boolean;
    readonly meta: boolean;
  };
}

export interface PointerSources {
  /** The 3D hover label, as lines (`labels.highlight`). */
  readonly hoverLines: ReadonlyState<readonly string[]>;
  readonly hover: Observable<ViewerPointerEvent>;
  readonly click: Observable<ViewerPointerEvent>;
}

export function createPointerSources(scope: SourceScope, plugin: PluginContext): PointerSources {
  const labels = fromBehavior(scope, plugin.behaviors.labels.highlight);
  const lines = mapState(labels, (l) => hoverLines(l.labels));
  const toEvent = (e: {
    current: { loci: unknown };
    button: number;
    modifiers: { shift: boolean; control: boolean; alt: boolean; meta: boolean };
  }): ViewerPointerEvent => {
    const loci = e.current.loci;
    const base = {
      lines: hoverLines(plugin.behaviors.labels.highlight.value.labels),
      button: e.button,
      modifiers: e.modifiers,
    };
    if (!StructureElement.Loci.is(loci) || StructureElement.Loci.isEmpty(loci)) return base;
    const selected = summarize(loci, { maxRangesPerChain: 4 });
    const list = plugin.managers.structure.hierarchy.current.structures;
    const ref = list.find((s) => s.cell.obj?.data.root === loci.structure.root);
    return {
      ...base,
      target: selected.addresses.join(', '),
      ...(ref ? { structure: structureNames(list).get(ref.cell.transform.ref) } : {}),
    };
  };
  return {
    hoverLines: lines,
    hover: plugin.behaviors.interaction.hover.pipe(map(toEvent), share()),
    click: plugin.behaviors.interaction.click.pipe(map(toEvent), share()),
  };
}
