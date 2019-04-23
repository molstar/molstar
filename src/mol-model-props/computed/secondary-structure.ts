/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CustomPropertyDescriptor, Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { DSSPComputationParams, computeUnitDSSP } from './secondary-structure/dssp';
import { SecondaryStructure } from 'mol-model/structure/model/properties/seconday-structure';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Unit } from 'mol-model/structure/structure';
import { idFactory } from 'mol-util/id-factory';

const nextSecondaryStructureId = idFactory()

export namespace ComputedSecondaryStructure {
    export type Property = {
        id: number
        map: Map<number, SecondaryStructure>
    }

    export function get(structure: Structure): Property | undefined {
        return structure.inheritedPropertyData.__ComputedSecondaryStructure__;
    }
    function set(structure: Structure, prop: Property) {
        (structure.inheritedPropertyData.__ComputedSecondaryStructure__ as Property) = prop;
    }

    export function createAttachTask(params: Partial<SecondaryStructureComputationProps> = {}) {
        return (structure: Structure) => Task.create('Compute Secondary Structure', async ctx => {
            if (get(structure)) return true;
            return await attachFromCifOrCompute(structure, params)
        });
    }

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_computed_secondary_structure',
        // TODO `cifExport` and `symbol`
    });

    export async function attachFromCifOrCompute(structure: Structure, params: Partial<SecondaryStructureComputationProps> = {}) {
        if (structure.customPropertyDescriptors.has(Descriptor)) return true;

        const compSecStruc = await computeSecondaryStructure(structure, params)

        structure.customPropertyDescriptors.add(Descriptor);
        set(structure, compSecStruc);
        return true;
    }
}

export const SecondaryStructureComputationParams = {
    ...DSSPComputationParams
}
export type SecondaryStructureComputationParams = typeof SecondaryStructureComputationParams
export type SecondaryStructureComputationProps = PD.Values<SecondaryStructureComputationParams>

async function computeSecondaryStructure(structure: Structure, params: Partial<SecondaryStructureComputationProps>): Promise<ComputedSecondaryStructure.Property> {
    const p = { ...PD.getDefaultValues(SecondaryStructureComputationParams), params }
    // TODO take inter-unit hbonds into account
    // TODO use Zhang-Skolnik for CA alpha only parts or for coarse parts with per-residue elements
    const map = new Map<number, SecondaryStructure>()
    for (let i = 0, il = structure.unitSymmetryGroups.length; i < il; ++i) {
        const u = structure.unitSymmetryGroups[i].units[0]
        if (Unit.isAtomic(u)) {
            const secondaryStructure = await computeUnitDSSP(u, p)
            map.set(u.invariantId, secondaryStructure)
        }
    }
    return { id: nextSecondaryStructureId(), map }
}