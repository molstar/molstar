/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition'
import { ShrakeRupleyComputationParams, AccessibleSurfaceArea } from './accessible-surface-area/shrake-rupley';
import { Structure, CustomPropertyDescriptor } from '../../mol-model/structure';
import { RuntimeContext } from '../../mol-task';
import { CustomStructureProperty } from '../common/custom-property-registry';

export const AccessibleSurfaceAreaParams = {
    ...ShrakeRupleyComputationParams
}
export type AccessibleSurfaceAreaParams = typeof AccessibleSurfaceAreaParams
export type AccessibleSurfaceAreaProps = PD.Values<AccessibleSurfaceAreaParams>

export type AccessibleSurfaceAreaValue = AccessibleSurfaceArea

export const AccessibleSurfaceAreaProvider: CustomStructureProperty.Provider<AccessibleSurfaceAreaParams, AccessibleSurfaceAreaValue> = CustomStructureProperty.createProvider({
    label: 'Accessible Surface Area',
    descriptor: CustomPropertyDescriptor({
        isStatic: true,
        name: 'molstar_accessible_surface_area',
        // TODO `cifExport` and `symbol`
    }),
    defaultParams: AccessibleSurfaceAreaParams,
    getParams: (data: Structure) => AccessibleSurfaceAreaParams,
    isApplicable: (data: Structure) => true,
    compute: async (ctx: RuntimeContext, data: Structure, props: Partial<AccessibleSurfaceAreaProps>) => {
        const p = { ...PD.getDefaultValues(AccessibleSurfaceAreaParams), ...props }
        return await AccessibleSurfaceArea.compute(data, p).runInContext(ctx)
    }
})