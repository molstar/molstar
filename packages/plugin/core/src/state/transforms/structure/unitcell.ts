/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { UnitcellParams, getUnitcellData, UnitcellRepresentation } from '@molstar/graphics/repr/shape/model/unitcell';
import { ModelSymmetry } from '@molstar/model/formats/structure/property/symmetry';
import type { PluginContext } from '@molstar/plugin/context';
import { Task } from '@molstar/core/task';
import { StateObject, StateTransformer } from '@molstar/core/state';

export { ModelUnitcell3D };
type ModelUnitcell3D = typeof ModelUnitcell3D;
const ModelUnitcell3D = PluginStateTransform.BuiltIn({
  name: 'model-unitcell-3d',
  display: 'Model Unit Cell',
  from: SO.Molecule.Model,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...UnitcellParams,
  }),
})({
  isApplicable: (a) => !!ModelSymmetry.Provider.get(a.data),
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Model Unit Cell', async (ctx) => {
      const symmetry = ModelSymmetry.Provider.get(a.data);
      if (!symmetry) return StateObject.Null;
      const data = getUnitcellData(a.data, symmetry.spacegroup.cell, params);
      const repr = UnitcellRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => UnitcellParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D(
        { repr, sourceData: data },
        { label: `Unit Cell`, description: symmetry.spacegroup.name },
      );
    });
  },
  update({ a, b, newParams }) {
    return Task.create('Model Unit Cell', async (ctx) => {
      const symmetry = ModelSymmetry.Provider.get(a.data);
      if (!symmetry) return StateTransformer.UpdateResult.Null;
      const props = { ...b.data.repr.props, ...newParams };
      const data = getUnitcellData(a.data, symmetry.spacegroup.cell, props);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      return StateTransformer.UpdateResult.Updated;
    });
  },
});
