/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra';
import { type State, type StateObjectCell, StateSelection } from '@molstar/core/state';
import { Particle, type ParticleList } from '@molstar/model/model/particles/particle-list';
import { getTransformFromParams, type TransformParam } from '@molstar/plugin/state/transforms/helpers';
import { TransformParticleList } from '@molstar/plugin/state/transforms/particles/transform';

// A particle list moved by hand (Mol*'s `TransformParticleList`, a decorator on the list, as a structure's conformation
// transform): its decorator and the matrix it applies, from where the list was loaded.

/** A list's transform decorator (`TransformParticleList` below it), if any */
export function particleTransformCell(state: State, listRef: string): StateObjectCell | undefined {
  return state.select(StateSelection.Generators.byRef(listRef).subtree().withTransformer(TransformParticleList))[0];
}

/** The matrix of a transform param, turning about the list's own center where it says so (as the decorator) */
export function particleMatrixOf(t: TransformParam, data: ParticleList): Mat4 {
  return getTransformFromParams(t, Particle.getBoundary(data).sphere.center);
}

/** A list's transform now (identity without one); `data` is the list as loaded (the hierarchy's, before it) */
export function currentParticleTransform(state: State, listRef: string, data: ParticleList): Mat4 {
  const t = particleTransformCell(state, listRef)?.transform.params?.transform as TransformParam | undefined;
  return t ? particleMatrixOf(t, data) : Mat4.identity();
}
