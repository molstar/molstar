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
import { cancelFrameCallback, requestFrameCallback } from '@molstar/prism-ui/reactive/frame-state';

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
  /** Brought into view in the sequence when set: it focuses on its owner's marks together (both ranges of a box). */
  readonly focusSequence?: boolean;
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

/** The marks to read (`viewer.sources.marks`): a view sets its marks through `highlight.mark` and `unmark`. */
export interface MarksView {
  readonly marks: ReadonlyState<readonly Mark[]>;
  readonly hover: ReadonlyState<MarkHover | undefined>;
  /** The tones' colors on the canvas as it is (light or dark): what a view draws its own marks in. */
  readonly colors: ReadonlyState<MarkColors>;
}

/** The marks and their writers, which only the API's `highlight` uses. */
export interface MarksSource extends MarksView {
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
  // Equal tones are not emitted: `settingsUpdated` fires for every canvas setting and camera focus, not only the theme.
  const colors = fromChanged(
    scope,
    themeOf,
    plugin.events.canvas3d.settingsUpdated,
    (a, b) => a.highlight === b.highlight && a.accent === b.accent && a.rest === b.rest,
  );
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
  // A representation that is gone has no children in the tree.
  for (const child of tree.children.get(reprRef)?.values() ?? []) {
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
 * When the marks or their colors change, every representation is painted again; otherwise only those whose own cell
 * or overpaint cells changed (a volume streamed in, or another representation's update, paints nothing).
 * A removed representation is let go once the task ends, not kept until the next paint. Returns the stop; `held`
 * counts the representations it still refers to (tests).
 */
export function paintMarks(
  plugin: PluginContext,
  marks: ReadonlyState<readonly Mark[]>,
  colors: ReadonlyState<MarkColors>,
): (() => void) & { readonly held: () => number } {
  /** What it painted, to put back. */
  let painted = new Map<Representation<any>, string>();
  /** The representations of marked structures it left alone (none of the marks on them, nothing dimmed), by ref. */
  let unpainted = new Map<Representation<any>, string>();
  /** Paint every representation again (the marks or their colors changed), not only the changed ones. */
  let all = false;
  /** The cells made or updated since the last paint: representations and their overpaint cells. */
  const changed = new Set<string>();
  /** Whether a representation's cell or one of its overpaint cells changed. */
  const touched = (ref: string) => {
    if (changed.has(ref)) return true;
    const children = plugin.state.data.tree.children.get(ref);
    if (children) for (const child of children.values()) if (changed.has(child)) return true;
    return false;
  };
  let frame = 0;
  const paint = () => {
    frame = 0;
    const list = marks.value;
    const next = new Map<Representation<any>, string>();
    const nextUnpainted = new Map<Representation<any>, string>();
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
            const ref = r.cell.transform.ref;
            // As it was painted (or left) and nothing of it changed: it keeps that.
            if (!all && !touched(ref)) {
              if (painted.has(repr)) {
                next.set(repr, ref);
                continue;
              }
              if (unpainted.has(repr)) {
                nextUnpainted.set(repr, ref);
                continue;
              }
            }
            const layers: Overpaint.Layer<StructureElement.Loci>[] = [];
            for (const m of own) {
              const loci = StructureElement.Loci.remap(m.loci, source);
              if (!StructureElement.Loci.isEmpty(loci)) layers.push({ loci, color: Color(tone[m.tone]), clear: false });
            }
            // Dimmed: grey under them (all of it: the marks paint over), a representation with none of them too.
            if (!layers.length && !dim) {
              nextUnpainted.set(repr, ref);
              continue;
            }
            if (dim)
              layers.unshift({ loci: Structure.toStructureElementLoci(source), color: Color(tone.rest), clear: false });
            const base = ownOverpaint(plugin, ref);
            repr.setState({ overpaint: Overpaint('element-loci', [...base.layers, ...layers] as never) });
            plugin.canvas3d?.update(repr);
            next.set(repr, ref);
          }
        }
      }
    }
    // What it painted before and does not now: its own overpaint back (nothing for a representation that is gone).
    for (const [repr, ref] of painted) {
      if (next.has(repr) || !plugin.state.data.tree.transforms.has(ref)) continue;
      repr.setState({ overpaint: ownOverpaint(plugin, ref) });
      plugin.canvas3d?.update(repr);
    }
    painted = next;
    unpainted = nextUnpainted;
    all = false;
    changed.clear();
  };
  const schedule = () => {
    if (frame) return;
    frame = typeof requestAnimationFrame === 'function' ? requestFrameCallback(paint) : (queueMicrotask(paint), -1);
  };
  const repaintAll = () => {
    if (!marks.value.length && !painted.size) return;
    all = true;
    schedule();
  };
  const events = plugin.state.data.events;
  const subs = [
    marks.subscribe(repaintAll),
    colors.subscribe(repaintAll),
    // A structure's update updates its representations, so theirs are the events that matter.
    merge(events.object.created, events.object.updated).subscribe((e) => {
      if (!marks.value.length && !painted.size) return;
      if (!SO.Molecule.Structure.Representation3D.is(e.obj) && !SO.Molecule.Structure.Representation3DState.is(e.obj)) {
        return;
      }
      changed.add(e.ref);
      schedule();
    }),
    // A removed representation (its structure and theme with it) is let go once the removal's task ends: a bulk removal
    // prunes once. Nothing to put back for it.
    events.cell.removed.subscribe(() => {
      if (pruning || (!painted.size && !unpainted.size)) return;
      pruning = true;
      queueMicrotask(prune);
    }),
  ];
  let pruning = false;
  const prune = () => {
    pruning = false;
    const transforms = plugin.state.data.tree.transforms;
    for (const [repr, ref] of painted) if (!transforms.has(ref)) painted.delete(repr);
    for (const [repr, ref] of unpainted) if (!transforms.has(ref)) unpainted.delete(repr);
  };
  const stop = () => {
    for (const s of subs) s.unsubscribe();
    if (frame > 0) cancelFrameCallback(frame);
    for (const [repr, ref] of painted) {
      if (plugin.state.data.tree.transforms.has(ref)) repr.setState({ overpaint: ownOverpaint(plugin, ref) });
    }
    painted.clear();
    unpainted.clear();
  };
  return Object.assign(stop, { held: () => painted.size + unpainted.size });
}

/** The read side of the marks, for `viewer.sources` (read-only: api.md §6). */
export function marksView(m: MarksView): MarksView {
  return { marks: m.marks, hover: m.hover, colors: m.colors };
}
