/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */
import { Structure } from '../../../mol-model/structure/structure';
import { Task, RuntimeContext } from '../../../mol-task';
import { MembraneOrientation } from '../../../mol-model/structure/model/properties/membrane-orientation';

export function computeOPM(structure: Structure) {
    return Task.create('Parse Membrane Topology', async runtime => {
        return await calculate(runtime, structure);
    });
}

export async function calculate(runtime: RuntimeContext, structure: Structure): Promise<MembraneOrientation> {
    throw Error('impl me');
}