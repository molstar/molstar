/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type State, type StateObjectCell, StateTransform } from '@molstar/core/state';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ClipVolumesTag } from '@molstar/plugin/state/transforms/shape/clip-volumes';

// The shapes in the scene (api.md §2 `shape`, the Structure panel's Shapes group): Mol*'s shape providers loaded on their
// own (a mesh file: OBJ, PLY, VTP; a primitives node), not those of a structure, a map or particles (a unit cell, a
// cut's drawing), each with the representations made of it. Shapes have no hierarchy manager: read from the state tree.

export interface ShapeEntry {
  /** The shape provider's cell. */
  readonly cell: StateObjectCell;
  /** What it is called: the file it came from (the provider's own label is "Shape"). */
  readonly label: string;
  /** Its representations (`ShapeRepresentation3D` cells under it). */
  readonly representations: readonly StateObjectCell[];
}

/** Whether a cell lies under a structure, a map or particles (their shapes are theirs: a unit cell, a cut's drawing). */
function underData(state: State, cell: StateObjectCell): boolean {
  for (let ref = cell.transform.parent; ; ) {
    const c = state.cells.get(ref);
    if (!c || ref === StateTransform.RootRef) return false;
    const obj = c.obj;
    if (
      SO.Molecule.Structure.is(obj) ||
      SO.Molecule.Model.is(obj) ||
      SO.Volume.Data.is(obj) ||
      SO.Particle.List.is(obj)
    )
      return true;
    if (c.transform.parent === ref) return false;
    ref = c.transform.parent;
  }
}

/** The shape's name: the label of the data it was parsed from (a file's name), else its own. */
function shapeLabel(state: State, cell: StateObjectCell): string {
  let label = cell.obj?.label ?? 'Shape';
  for (let ref = cell.transform.parent; ref !== StateTransform.RootRef; ) {
    const c = state.cells.get(ref);
    if (!c) break;
    if (c.obj?.type.typeClass === 'Data' && c.obj.label) label = c.obj.label;
    if (c.transform.parent === ref) break;
    ref = c.transform.parent;
  }
  return label;
}

/** The shapes loaded on their own, in tree order. */
export function shapeEntries(state: State): ShapeEntry[] {
  const out: ShapeEntry[] = [];
  const visit = (ref: string) => {
    for (const child of state.tree.children.get(ref).values()) {
      const cell = state.cells.get(child);
      if (!cell) continue;
      // a ghost (a parse step the tree hides: OBJ, PLY data) is no shape, but holds one
      if (
        !cell.state.isGhost &&
        SO.Shape.Provider.is(cell.obj) &&
        !cell.transform.tags?.includes(ClipVolumesTag) &&
        !underData(state, cell)
      ) {
        const representations: StateObjectCell[] = [];
        for (const r of state.tree.children.get(child).values()) {
          const rc = state.cells.get(r);
          if (rc && !rc.state.isGhost && SO.Shape.Representation3D.is(rc.obj)) representations.push(rc);
        }
        out.push({ cell, label: shapeLabel(state, cell), representations });
        continue;
      }
      // what holds no shape of its own is not looked into (a structure's subtree)
      if (SO.Molecule.Structure.is(cell.obj) || SO.Volume.Data.is(cell.obj) || SO.Particle.List.is(cell.obj)) continue;
      visit(child);
    }
  };
  visit(StateTransform.RootRef);
  return out;
}
