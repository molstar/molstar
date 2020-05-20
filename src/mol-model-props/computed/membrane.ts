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
import { ANVILParams, Membrane, ANVILProps } from './membrane/ANVIL';
import { AccessibleSurfaceAreaProvider } from './accessible-surface-area';

function getMembraneParams(data?: Structure) {
    let defaultType = 'anvil' as 'anvil' | 'opm'; // TODO flip - OPM should be default if PDB identifier is known
    return {
        type: PD.MappedStatic(defaultType, {
            'opm': PD.EmptyGroup({ label: 'OPM' }),
            'anvil': PD.Group(ANVILParams, { label: 'ANVIL' })
        }, { options: [['opm', 'OPM'], ['anvil', 'ANVIL']] })
    };
}

export const MembraneParams = getMembraneParams();
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
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<MembraneProps>) => {
        const p = { ...PD.getDefaultValues(MembraneParams), ...props };
        switch (p.type.name) {
            case 'anvil': return { value: await computeAnvil(ctx, data, p.type.params) };
            case 'opm': return { value: await computeOpm(data) };
        }
    }
});

async function computeAnvil(ctx: CustomProperty.Context, data: Structure, props: ANVILProps): Promise<Membrane> {
    await AccessibleSurfaceAreaProvider.attach(ctx, data);
    const p = { ...PD.getDefaultValues(MembraneParams), ...props };
    return await Membrane.compute(data, p).runInContext(ctx.runtime);
}

async function computeOpm(structure: Structure): Promise<Membrane> {
    throw Error('TODO impl');
}