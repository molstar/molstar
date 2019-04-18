/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Structure } from 'mol-model/structure';
import { Task } from 'mol-task';

export namespace ComputedSecondaryStructure {
    export type Property = {} // TODO

    export function get(structure: Structure): Property | undefined {
        return structure.currentPropertyData.__ComputedSecondaryStructure__;
    }
    function set(structure: Structure, prop: Property) {
        (structure.currentPropertyData.__ComputedSecondaryStructure__ as Property) = prop;
    }

    export function createAttachTask() {
        return (structure: Structure) => Task.create('Compute Secondary Structure', async ctx => {
            if (get(structure)) return true;
            return await attachFromCifOrCompute(structure, ctx)
        });
    }

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_computed_secondary_structure',
        // TODO `cifExport` and `symbol`
    });

    export async function attachFromCifOrCompute(structure: Structure, params: {
        // TODO params
    }) {
        if (structure.customPropertyDescriptors.has(Descriptor)) return true;

        const compSecStruc = computeSecondaryStructure(structure)

        structure.customPropertyDescriptors.add(Descriptor);
        set(structure, compSecStruc);
        return true;
    }
}

// export const SecondaryStructureComputationParams = {
//     oldDefinition: PD.Boolean(true, { description: 'Whether to use the old DSSP convention for the annotation of turns and helices, causes them to be two residues shorter' }),
//     oldOrdering: PD.Boolean(true, { description: 'Alpha-helices are preferred over 3-10 helices' })
// }
// export type SecondaryStructureComputationParams = typeof SecondaryStructureComputationParams

function computeSecondaryStructure(structure: Structure): ComputedSecondaryStructure.Property {
    // TODO
    console.log('TODO calc secondary structure')
    return {}
}