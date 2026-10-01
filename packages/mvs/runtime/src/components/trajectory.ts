/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '@molstar/plugin/state/objects';
import { getTrajectory } from '@molstar/plugin/state/transforms/model';
import { Task } from '@molstar/core/task';
import { ParamDefinition } from '@molstar/core/util/param-definition';
import { getMVSReferenceObject } from '@molstar/mvs/helpers/utils';
import { MVSTransform } from './annotation-structure-component.js';

export const MVSTrajectoryWithCoordinates = MVSTransform({
    name: 'trajectory-with-coordinates',
    display: { name: 'Trajectory with Coordinates', description: 'Create a trajectory from existing model and the provided coordinates.' },
    from: [PluginStateObject.Molecule.Model, PluginStateObject.Molecule.Topology],
    to: PluginStateObject.Molecule.Trajectory,
    params: {
        coordinatesRef: ParamDefinition.Text('', { isHidden: true }),
    }
})({
    apply({ a, params, dependencies }) {
        return Task.create('Create trajectory from model/topology and coordinates', async ctx => {
            const coordinates = getMVSReferenceObject([PluginStateObject.Molecule.Coordinates], dependencies, params.coordinatesRef);

            if (!coordinates) {
                throw new Error('Coordinates not found.');
            }

            const trajectory = await getTrajectory(ctx, a, coordinates.data);
            const props = { label: 'Trajectory', description: `${trajectory.frameCount} model${trajectory.frameCount === 1 ? '' : 's'}` };
            return new PluginStateObject.Molecule.Trajectory(trajectory, props);
        });
    }
});