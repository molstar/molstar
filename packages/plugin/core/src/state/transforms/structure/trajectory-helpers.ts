/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import type { RuntimeContext } from '@molstar/core/task';
import type { StateObject } from '@molstar/core/state';
import { type Coordinates, type Topology, Model } from '@molstar/model/model/structure';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';

export async function getTrajectory(ctx: RuntimeContext, obj: StateObject, coordinates: Coordinates) {
  if (obj.type === SO.Molecule.Topology.type) {
    const topology = obj.data as Topology;
    return await Model.trajectoryFromTopologyAndCoordinates(topology, coordinates).runInContext(ctx);
  } else if (obj.type === SO.Molecule.Model.type) {
    const model = obj.data as Model;
    return Model.trajectoryFromModelAndCoordinates(model, coordinates);
  }
  throw new Error('no model/topology found');
}
