/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateTransform, PluginStateObject as SO } from '@molstar/plugin/state/objects';
import {
  UnitcellParams,
  getUnitcellDataFromSymmetry,
  UnitcellRepresentation,
} from '@molstar/graphics/repr/shape/model/unitcell';
import type { PluginContext } from '@molstar/plugin/context';
import { Task } from '@molstar/core/task';
import { StateObject } from '@molstar/core/state/object';
import { Particle } from '@molstar/model/model/particles/particle-list';
import { Cell } from '@molstar/core/math/geometry/spacegroup/cell';
import { StateTransformer } from '@molstar/core/state';

export { ParticleListUnitcell3D };
type ParticleListUnitcell3D = typeof ParticleListUnitcell3D;
const ParticleListUnitcell3D = PluginStateTransform.BuiltIn({
  name: 'particle-list-unitcell-3d',
  display: 'Particle List Unit Cell',
  from: SO.Particle.List,
  to: SO.Shape.Representation3D,
  params: () => ({
    ...UnitcellParams,
  }),
})({
  isApplicable: (a) => !!a.data.cell,
  canAutoUpdate({ oldParams, newParams }) {
    return true;
  },
  apply({ a, params }, plugin: PluginContext) {
    return Task.create('Particle List Unit Cell', async (ctx) => {
      const { cell } = a.data;
      if (!cell) return StateObject.Null;

      const center = Particle.getBoundary(a.data).sphere.center;
      const data = getUnitcellDataFromSymmetry(cell, center, params);
      const repr = UnitcellRepresentation(
        { webgl: plugin.canvas3d?.webgl, ...plugin.representation.structure.themes },
        () => UnitcellParams,
      );
      await repr.createOrUpdate(params, data).runInContext(ctx);
      return new SO.Shape.Representation3D(
        { repr, sourceData: data },
        { label: 'Unit Cell', description: Cell.getLabel(cell) },
      );
    });
  },
  update({ a, b, newParams }) {
    return Task.create('Particle List Unit Cell', async (ctx) => {
      const { cell } = a.data;
      if (!cell) return StateTransformer.UpdateResult.Null;

      const props = { ...b.data.repr.props, ...newParams };
      const center = Particle.getBoundary(a.data).sphere.center;
      const data = getUnitcellDataFromSymmetry(cell, center, props);
      await b.data.repr.createOrUpdate(props, data).runInContext(ctx);
      b.data.sourceData = data;
      b.description = Cell.getLabel(cell);
      return StateTransformer.UpdateResult.Updated;
    });
  },
});
