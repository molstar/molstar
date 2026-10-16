/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Representation } from '@molstar/graphics/repr/representation';
import { StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { distinctUntilChanged, map, type Observable, share } from 'rxjs';
import { summarize } from '../targets/index.js';
import { hoverLines } from './hover-lines.js';
import { structureNames } from './naming.js';
import { fromChanged, type SourceScope } from './state.js';

// What is under the pointer (api.md §8: v5's `structureInteractivity`; README decision 22): the hover label's lines,
// and the hover and click events embedding hosts read, each with the target in the address language (targets.md §1).

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
  // The lines are derived on every emission and compared themselves: Mol*'s label manager reuses one `labels` array
  // (`LociLabelManager.getLabels`), so the event's value looks unchanged to a shallow compare.
  const lines = fromChanged(
    scope,
    () => hoverLines(plugin.behaviors.labels.highlight.value.labels),
    plugin.behaviors.labels.highlight,
    sameLines,
  );
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
    // Mol* emits the hover after every pick while the pointer moves: the same element under the same keys is one event.
    hover: plugin.behaviors.interaction.hover.pipe(distinctUntilChanged(sameHover), map(toEvent), share()),
    click: plugin.behaviors.interaction.click.pipe(map(toEvent), share()),
  };
}

type Modifiers = ViewerPointerEvent['modifiers'];

function sameHover(
  a: { current: Representation.Loci; button: number; modifiers: Modifiers },
  b: { current: Representation.Loci; button: number; modifiers: Modifiers },
): boolean {
  return (
    a.button === b.button &&
    a.modifiers.shift === b.modifiers.shift &&
    a.modifiers.control === b.modifiers.control &&
    a.modifiers.alt === b.modifiers.alt &&
    a.modifiers.meta === b.modifiers.meta &&
    Representation.Loci.areEqual(a.current, b.current)
  );
}

function sameLines(a: readonly string[], b: readonly string[]): boolean {
  return a.length === b.length && a.every((l, i) => l === b[i]);
}
