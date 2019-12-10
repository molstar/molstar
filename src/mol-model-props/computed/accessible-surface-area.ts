/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ShrakeRupleyComputationParams, AccessibleSurfaceArea } from './accessible-surface-area/shrake-rupley';
import { Structure, CustomPropertyDescriptor } from '../../mol-model/structure';
import { Task } from '../../mol-task';
import { idFactory } from '../../mol-util/id-factory';


const nextAccessibleSurfaceAreaId = idFactory()

export namespace ComputedAccessibleSurfaceArea {
    export type Property = {
        id: number
        asa: AccessibleSurfaceArea
    }

    export function get(structure: Structure): Property | undefined {
        return structure.inheritedPropertyData.__ComputedAccessibleSurfaceArea__;
    }
    function set(structure: Structure, prop: Property) {
        (structure.inheritedPropertyData.__ComputedAccessibleSurfaceArea__ as Property) = prop;
    }

    export function createAttachTask(params: Partial<AccessibleSurfaceAreaComputationProps> = {}) {
        return (structure: Structure) => Task.create('Compute Accessible Surface Area', async ctx => {
            if (get(structure)) return true;
            return await attachFromCifOrCompute(structure, params)
        });
    }

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_accessible_surface_area',
        // TODO `cifExport` and `symbol`
    });

    export async function attachFromCifOrCompute(structure: Structure, params: Partial<AccessibleSurfaceAreaComputationProps> = {}) {
        if (structure.customPropertyDescriptors.has(Descriptor)) return true;

        const compAccessibleSurfaceArea = await computeAccessibleSurfaceArea(structure, params)

        structure.customPropertyDescriptors.add(Descriptor);
        set(structure, compAccessibleSurfaceArea);
        return true;
    }
}

export const AccessibleSurfaceAreaComputationParams = {
    ...ShrakeRupleyComputationParams
}
export type AccessibleSurfaceAreaComputationParams = typeof AccessibleSurfaceAreaComputationParams
export type AccessibleSurfaceAreaComputationProps = PD.Values<AccessibleSurfaceAreaComputationParams>

async function computeAccessibleSurfaceArea(structure: Structure, params: Partial<AccessibleSurfaceAreaComputationProps>): Promise<ComputedAccessibleSurfaceArea.Property> {
    const p = { ...PD.getDefaultValues(AccessibleSurfaceAreaComputationParams), params };

    const accessibleSurfaceArea = await AccessibleSurfaceArea.compute(structure, p);
    return { id: nextAccessibleSurfaceAreaId(), asa: accessibleSurfaceArea }
}