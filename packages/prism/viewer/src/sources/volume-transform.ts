/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, type Vec3 } from '@molstar/core/math/linear-algebra';
import { type State, type StateObjectCell, StateSelection, StateTransform } from '@molstar/core/state';
import { Grid, type Volume } from '@molstar/model/model/volume';
import { PluginStateObject } from '@molstar/plugin/state/objects';
import { getTransformFromParams, type TransformParam } from '@molstar/plugin/state/transforms/helpers';
import { VolumeTransform } from '@molstar/plugin/state/transforms/volume/ops';

// A volume moved by hand (Mol*'s `VolumeTransform`, a decorator on the volume, as a structure's conformation
// transform): its decorator and the matrix it applies, from where the volume was loaded.

/**
 * Whether a volume's `cell` is under a structure or a model: density streamed for it, which moves with it (its
 * transform), not by itself
 */
export function isStreamedForStructure(state: State, cell: StateObjectCell): boolean {
  for (let ref = cell.transform.parent; ref !== StateTransform.RootRef; ) {
    const c = state.cells.get(ref);
    if (!c) return false;
    if (PluginStateObject.Molecule.Structure.is(c.obj) || PluginStateObject.Molecule.Model.is(c.obj)) return true;
    if (c.transform.parent === ref) return false;
    ref = c.transform.parent;
  }
  return false;
}

/** A volume's transform decorator (`VolumeTransform` below it), if any */
export function volumeTransformCell(state: State, volumeRef: string): StateObjectCell | undefined {
  return state.select(StateSelection.Generators.byRef(volumeRef).subtree().withTransformer(VolumeTransform))[0];
}

/** A volume's center as loaded: what its decorator turns about */
export function volumeCenter(data: Volume): Vec3 {
  return Grid.getBoundingSphere(data.grid).center;
}

/** The matrix of a transform param, turning about the volume's own center where it says so (as the decorator) */
export function volumeMatrixOf(t: TransformParam, data: Volume): Mat4 {
  return getTransformFromParams(t, volumeCenter(data));
}

/** A volume's transform now (identity without one); `data` is the volume as loaded (the hierarchy's, before it) */
export function currentVolumeTransform(state: State, volumeRef: string, data: Volume): Mat4 {
  const t = volumeTransformCell(state, volumeRef)?.transform.params?.transform as TransformParam | undefined;
  return t ? volumeMatrixOf(t, data) : Mat4.identity();
}
