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

export const MVSTrajectoryFromModelAndCoordinates = MVSTransform({
    name: 'trajectory-from-model-and-coordinates',
    display: { name: 'Trajectory from Topology & Coordinates', description: 'Create a trajectory from existing model/topology and coordinates.' },
    from: PluginStateObject.Root,
    to: PluginStateObject.Molecule.Trajectory,
    params: {
        modelRef: ParamDefinition.Text('', { isHidden: true }),
        coordinatesRef: ParamDefinition.Text('', { isHidden: true }),
    }
})({
    apply({ params, dependencies }) {
        return Task.create('Create trajectory from model/topology and coordinates', async ctx => {
            const model = getMVSReferenceObject([PluginStateObject.Molecule.Model], dependencies, params.modelRef);
            const coordinates = getMVSReferenceObject([PluginStateObject.Molecule.Coordinates], dependencies, params.coordinatesRef);

            if (!model || !coordinates) {
                throw new Error('Model and coordinates are required to create a trajectory.');
            }

            const trajectory = await getTrajectory(ctx, model, coordinates.data);
            const props = { label: 'Trajectory', description: `${trajectory.frameCount} model${trajectory.frameCount === 1 ? '' : 's'}` };
            return new PluginStateObject.Molecule.Trajectory(trajectory, props);
        });
    }
});