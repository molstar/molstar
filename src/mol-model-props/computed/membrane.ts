/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Structure } from '../../mol-model/structure';
import { CustomStructureProperty } from '../common/custom-structure-property';
import { CustomProperty } from '../common/custom-property';
import { CustomPropertyDescriptor } from '../../mol-model/custom-property';
import { ANVILParams, Membrane } from './membrane/ANVIL';
import { AccessibleSurfaceAreaProvider } from './accessible-surface-area';

export const MembraneParams = {
    ...ANVILParams
};
export type MembraneParams = typeof MembraneParams
export type MembraneProps = PD.Values<MembraneParams>

export const MembraneProvider: CustomStructureProperty.Provider<MembraneParams, Membrane> = CustomStructureProperty.createProvider({
    label: 'Predicted Membrane',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_membrane',
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: MembraneParams,
    getParams: (data: Structure) => MembraneParams,
    isApplicable: (data: Structure) => true,
    // TODO needs ASA to be computed (or 'resolved' before trying computing topology) - how to achieve?
    // TODO potentially, this could behave like secondary structure info where data can be either parsed or computed
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<MembraneProps>) => {
        await AccessibleSurfaceAreaProvider.attach(ctx, data);
        const p = { ...PD.getDefaultValues(MembraneParams), ...props };
        return { value: await Membrane.compute(data, p).runInContext(ctx.runtime) };
    }
});