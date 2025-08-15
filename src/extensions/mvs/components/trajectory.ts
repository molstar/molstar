/**
 * Copyright (c) 2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '../../../mol-plugin-state/objects';
import { getTrajectory } from '../../../mol-plugin-state/transforms/model';
import { Task } from '../../../mol-task';
import { ParamDefinition } from '../../../mol-util/param-definition';
import { getMVSReferenceObject } from '../helpers/utils';
import { MVSTransform } from './annotation-structure-component';

export const MVSTrajectoryWithCoordinates = MVSTransform({
    name: 'trajectory-with-coordinates',
    display: { name: 'Trajectory with Coordinates', description: 'Create a trajectory from existing model and the provided coordinates.' },
    from: PluginStateObject.Molecule.Model,
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