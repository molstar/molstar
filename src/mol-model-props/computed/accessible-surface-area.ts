/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ShrakeRupleyComputationParams, AccessibleSurfaceArea } from './accessible-surface-area/shrake-rupley';
import { Structure, CustomPropertyDescriptor } from '../../mol-model/structure';
import { Task, RuntimeContext } from '../../mol-task';
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
        return (structure: Structure) => attachTask(structure, params)
    }

    export function attachTask(structure: Structure, params: Partial<AccessibleSurfaceAreaComputationProps> = {}) {
        return Task.create('Compute Accessible Surface Area', async ctx => {
            if (get(structure)) return;
            return await attachFromCifOrCompute(ctx, structure, params)
        });
    }

    export const Descriptor = CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_accessible_surface_area',
        // TODO `cifExport` and `symbol`
    });

    export async function attachFromCifOrCompute(ctx: RuntimeContext, structure: Structure, params: Partial<AccessibleSurfaceAreaComputationProps> = {}) {
        if (structure.customPropertyDescriptors.has(Descriptor)) return;

        const compAccessibleSurfaceArea = await computeAccessibleSurfaceArea(ctx, structure, params)

        structure.customPropertyDescriptors.add(Descriptor);
        set(structure, compAccessibleSurfaceArea);
    }
}

export const AccessibleSurfaceAreaComputationParams = {
    ...ShrakeRupleyComputationParams
}
export type AccessibleSurfaceAreaComputationParams = typeof AccessibleSurfaceAreaComputationParams
export type AccessibleSurfaceAreaComputationProps = PD.Values<AccessibleSurfaceAreaComputationParams>

async function computeAccessibleSurfaceArea(ctx: RuntimeContext, structure: Structure, params: Partial<AccessibleSurfaceAreaComputationProps>): Promise<ComputedAccessibleSurfaceArea.Property> {
    const p = { ...PD.getDefaultValues(AccessibleSurfaceAreaComputationParams), params };

    const asa = await AccessibleSurfaceArea.compute(structure, p).runInContext(ctx);
    return { id: nextAccessibleSurfaceAreaId(), asa }
}