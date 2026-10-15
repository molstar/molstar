/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, StructureElement } from '@molstar/model/model/structure';
import type { Representation } from '@molstar/graphics/repr/representation';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { Color } from '@molstar/core/util/color';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { themeForColor } from '@molstar/prism-ui/model/canvas-theme';
import { merge } from 'rxjs';
import { fromChanged, type SourceScope } from './state.js';

// Marks (marks.md): named residue ranges a view puts on the scene for a while, shown at once in the 3D view and in the
// sequence, each in a tone: `highlight` (yellow) or `accent` (the UI's accent). The PAE plot's box marks its aligned
// range in highlight and its scored range in accent; another view (a contact map, an alignment) can mark its own.
// They are view state: never a history step, never in the state tree, a snapshot or a session. An owner sets its marks
// whole and clears them. Beside them, the hover: the residues a view's pointer is on (a plot's cell), lit in the
// sequence as the 3D view's hover is (the 3D view by `highlight.set`).
//
// The 3D view: the marks' overpaint on every structure representation, over the representation's own (a style layer's
// overpaint stays under them), later marks over earlier; an owner may grey the rest of the structure (`dimRest`, as
// AlphaFold DB greys what its PAE box does not cover), so the marks stand out of any coloring. Mol* sets a
// representation's overpaint again whenever its overpaint cells change, so the marks are put back after it. The
// sequence: a highlight channel each (sequence/link.ts), later over earlier there too.

export type MarkTone = 'highlight' | 'accent';

export interface Mark {
  /** Who set it (`pae`): its marks are set and cleared together. */
  readonly owner: string;
  readonly id: string;
  /** What it is, for a legend or the sequence's tooltip ("Aligned A 92–190"). */
  readonly label: string;
  readonly tone: MarkTone;
  /** The structure (its state ref) and its residues there. */
  readonly structure: string;
  readonly loci: StructureElement.Loci;
  /** The rest of its structure greyed in the 3D view while it is shown, so the marks stand out of any coloring. */
  readonly dimRest?: boolean;
}

export interface MarkHover {
  readonly owner: string;
  readonly structure: string;
  readonly loci: StructureElement.Loci;
}

/** The tones' colors for a canvas theme: the sequence's highlight and Prism's accent (sequence/skin.ts). */
export const MarkTones: Readonly<Record<'light' | 'dark', Readonly<Record<MarkTone, number>>>> = {
  light: { highlight: 0xf2b400, accent: 0xae5d04 },
  dark: { highlight: 0xfacc15, accent: 0x51a2fb },
};
/** The grey of the rest of a structure under marks that dim it (`dimRest`), per canvas theme. */
export const MarkRest: Readonly<Record<'light' | 'dark', number>> = { light: 0xc4c2bf, dark: 0x5a5856 };

export type MarkColors = Readonly<Record<MarkTone | 'rest', number>>;

export interface MarksSource {
  readonly marks: ReadonlyState<readonly Mark[]>;
  readonly hover: ReadonlyState<MarkHover | undefined>;
  /** The tones' colors on the canvas as it is (light or dark): what a view draws its own marks in. */
  readonly colors: ReadonlyState<MarkColors>;
  /** An owner's marks, replacing the ones it had. */
  set(owner: string, marks: readonly Omit<Mark, 'owner'>[]): void;
  clear(owner: string): void;
  /** The residues a view's pointer is on; undefined: none (only the owner's own hover is cleared). */
  setHover(owner: string, hover: Omit<MarkHover, 'owner'> | undefined): void;
}

export function createMarksSource(scope: SourceScope, plugin: PluginContext): MarksSource {
  const marks = observableState<readonly Mark[]>([]);
  const hover = observableState<MarkHover | undefined>(undefined);
  const themeOf = () => {
    const bg = plugin.canvas3d?.props.renderer.backgroundColor;
    const theme = bg === undefined ? 'light' : themeForColor(bg as number);
    return { ...MarkTones[theme], rest: MarkRest[theme] };
  };
  const colors = fromChanged(scope, themeOf, plugin.events.canvas3d.settingsUpdated);
  scope.add(paintMarks(plugin, marks.asReadonly, colors));
  return {
    marks: marks.asReadonly,
    hover: hover.asReadonly,
    colors,
    set(owner, list) {
      const others = marks.value.filter((m) => m.owner !== owner);
      marks.set([...others, ...list.map((m) => ({ ...m, owner }))]);
    },
    clear(owner) {
      if (marks.value.some((m) => m.owner === owner)) marks.set(marks.value.filter((m) => m.owner !== owner));
    },
    setHover(owner, h) {
      if (h) hover.set({ ...h, owner });
      else if (hover.value?.owner === owner) hover.set(undefined);
    },
  };
}

/** The overpaint a representation has of its own (its overpaint cells'), without the marks. */
function ownOverpaint(plugin: PluginContext, reprRef: string): Overpaint {
  const tree = plugin.state.data.tree;
  for (const child of tree.children.get(reprRef).values()) {
    const obj = plugin.state.data.cells.get(child)?.obj;
    if (!SO.Molecule.Structure.Representation3DState.is(obj)) continue;
    const state = (obj.data as SO.Molecule.Structure.Representation3DStateData).state as { overpaint?: Overpaint };
    if (state.overpaint) return state.overpaint;
  }
  return Overpaint.Empty;
}

/**
 * Paints the marks on the structure representations, over their own overpaint, and puts that back when the marks go;
 * again after Mol* sets a representation's state (its overpaint cells changed) or makes one. Once per frame at most.
 * Returns the stop.
 */
function paintMarks(
  plugin: PluginContext,
  marks: ReadonlyState<readonly Mark[]>,
  colors: ReadonlyState<MarkColors>,
): () => void {
  /** What it painted, to put back. */
  let painted = new Map<Representation<any>, string>();
  let frame = 0;
  const paint = () => {
    frame = 0;
    const list = marks.value;
    const next = new Map<Representation<any>, string>();
    if (list.length) {
      const tone = colors.value;
      for (const s of plugin.managers.structure.hierarchy.current.structures) {
        const own = list.filter((m) => m.structure === s.cell.transform.ref);
        if (!own.length) continue;
        const dim = own.some((m) => m.dimRest);
        for (const c of s.components) {
          for (const r of c.representations) {
            const repr = r.cell.obj?.data.repr;
            const source = r.cell.obj?.data.sourceData;
            if (!repr || !source) continue;
            const layers: Overpaint.Layer<StructureElement.Loci>[] = [];
            for (const m of own) {
              const loci = StructureElement.Loci.remap(m.loci, source);
              if (!StructureElement.Loci.isEmpty(loci)) layers.push({ loci, color: Color(tone[m.tone]), clear: false });
            }
            // Dimmed: grey under them (all of it: the marks paint over), a representation with none of them too.
            if (!layers.length && !dim) continue;
            if (dim)
              layers.unshift({ loci: Structure.toStructureElementLoci(source), color: Color(tone.rest), clear: false });
            const base = ownOverpaint(plugin, r.cell.transform.ref);
            repr.setState({ overpaint: Overpaint('element-loci', [...base.layers, ...layers] as never) });
            plugin.canvas3d?.update(repr);
            next.set(repr, r.cell.transform.ref);
          }
        }
      }
    }
    // What it painted before and does not now: its own overpaint back.
    for (const [repr, ref] of painted) {
      if (next.has(repr)) continue;
      repr.setState({ overpaint: ownOverpaint(plugin, ref) });
      plugin.canvas3d?.update(repr);
    }
    painted = next;
  };
  const schedule = () => {
    if (frame) return;
    if (!marks.value.length && !painted.size) return;
    frame = typeof requestAnimationFrame === 'function' ? requestAnimationFrame(paint) : (queueMicrotask(paint), -1);
  };
  const events = plugin.state.data.events;
  const subs = [
    marks.subscribe(schedule),
    colors.subscribe(schedule),
    merge(events.object.created, events.object.updated).subscribe(schedule),
  ];
  return () => {
    for (const s of subs) s.unsubscribe();
    if (frame > 0) cancelAnimationFrame(frame);
    for (const [repr, ref] of painted) repr.setState({ overpaint: ownOverpaint(plugin, ref) });
    painted.clear();
  };
}
