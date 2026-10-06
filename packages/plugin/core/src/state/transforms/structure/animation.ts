/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureUnitTransforms } from '@molstar/model/model/structure/structure/util/unit-transforms';
import {
  unwindStructureAssembly,
  explodeStructure,
  SpinStructureParams,
  getSpinStructureAxisAndOrigin,
  spinStructure,
} from '@molstar/plugin/state/animation/helpers';
import type { Structure } from '@molstar/model/model/structure';
import { StateTransformer } from '@molstar/core/state';

export { UnwindStructureAssemblyRepresentation3D };
type UnwindStructureAssemblyRepresentation3D = typeof UnwindStructureAssemblyRepresentation3D;
const UnwindStructureAssemblyRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'unwind-structure-assembly-representation-3d',
  display: 'Unwind Assembly 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: { t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }) },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const structure = a.data.sourceData;
    const unitTransforms = new StructureUnitTransforms(structure);
    unwindStructureAssembly(structure, unitTransforms, params.t);
    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { unitTransforms },
        initialState: { unitTransforms: new StructureUnitTransforms(structure) },
        info: structure,
        repr: a.data.repr,
      },
      { label: `Unwind T = ${params.t.toFixed(2)}` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const structure = b.data.info as Structure;
    if (a.data.sourceData !== structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    if (oldParams.t === newParams.t) return StateTransformer.UpdateResult.Unchanged;
    const unitTransforms = b.data.state.unitTransforms!;
    unwindStructureAssembly(structure, unitTransforms, newParams.t);
    b.label = `Unwind T = ${newParams.t.toFixed(2)}`;
    b.data.repr = a.data.repr;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { ExplodeStructureRepresentation3D };
type ExplodeStructureRepresentation3D = typeof ExplodeStructureRepresentation3D;
const ExplodeStructureRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'explode-structure-representation-3d',
  display: 'Explode 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: { t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }) },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const structure = a.data.sourceData;
    const unitTransforms = new StructureUnitTransforms(structure);
    explodeStructure(structure, unitTransforms, params.t, structure.root.boundary.sphere);
    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { unitTransforms },
        initialState: { unitTransforms: new StructureUnitTransforms(structure) },
        info: structure,
        repr: a.data.repr,
      },
      { label: `Explode T = ${params.t.toFixed(2)}` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const structure = a.data.sourceData;
    if (b.data.info !== structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    if (oldParams.t === newParams.t) return StateTransformer.UpdateResult.Unchanged;
    const unitTransforms = b.data.state.unitTransforms!;
    explodeStructure(structure, unitTransforms, newParams.t, structure.root.boundary.sphere);
    b.label = `Explode T = ${newParams.t.toFixed(2)}`;
    b.data.repr = a.data.repr;
    return StateTransformer.UpdateResult.Updated;
  },
});

export { SpinStructureRepresentation3D };
type SpinStructureRepresentation3D = typeof SpinStructureRepresentation3D;
const SpinStructureRepresentation3D = PluginStateTransform.BuiltIn({
  name: 'spin-structure-representation-3d',
  display: 'Spin 3D Representation',
  from: SO.Molecule.Structure.Representation3D,
  to: SO.Molecule.Structure.Representation3DState,
  params: {
    t: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
    ...SpinStructureParams,
  },
})({
  canAutoUpdate() {
    return true;
  },
  apply({ a, params }) {
    const structure = a.data.sourceData;
    const unitTransforms = new StructureUnitTransforms(structure);

    const { axis, origin } = getSpinStructureAxisAndOrigin(structure.root, params);
    spinStructure(structure, unitTransforms, params.t, axis, origin);
    return new SO.Molecule.Structure.Representation3DState(
      {
        state: { unitTransforms },
        initialState: { unitTransforms: new StructureUnitTransforms(structure) },
        info: structure,
        repr: a.data.repr,
      },
      { label: `Spin T = ${params.t.toFixed(2)}` },
    );
  },
  update({ a, b, newParams, oldParams }) {
    const structure = a.data.sourceData;
    if (b.data.info !== structure) return StateTransformer.UpdateResult.Recreate;
    if (a.data.repr !== b.data.repr) return StateTransformer.UpdateResult.Recreate;

    if (oldParams.t === newParams.t && oldParams.axis === newParams.axis && oldParams.origin === newParams.origin)
      return StateTransformer.UpdateResult.Unchanged;

    const unitTransforms = b.data.state.unitTransforms!;
    const { axis, origin } = getSpinStructureAxisAndOrigin(structure.root, newParams);
    spinStructure(structure, unitTransforms, newParams.t, axis, origin);
    b.label = `Spin T = ${newParams.t.toFixed(2)}`;
    b.data.repr = a.data.repr;
    return StateTransformer.UpdateResult.Updated;
  },
});
