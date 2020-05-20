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
import { ANVILParams, MembraneOrientation, ANVILProps } from './membrane-orientation/ANVIL';
import { AccessibleSurfaceAreaProvider } from './accessible-surface-area';

function getMembraneOrientationParams(data?: Structure) {
    let defaultType = 'anvil' as 'anvil' | 'opm'; // TODO flip - OPM should be default if PDB identifier is known
    return {
        type: PD.MappedStatic(defaultType, {
            'opm': PD.EmptyGroup({ label: 'OPM' }),
            'anvil': PD.Group(ANVILParams, { label: 'ANVIL' })
        }, { options: [['opm', 'OPM'], ['anvil', 'ANVIL']] })
    };
}

export const MembraneOrientationParams = getMembraneOrientationParams();
export type MembraneOrientationParams = typeof MembraneOrientationParams
export type MembraneOrientationProps = PD.Values<MembraneOrientationParams>

export const MembraneOrientationProvider: CustomStructureProperty.Provider<MembraneOrientationParams, MembraneOrientation> = CustomStructureProperty.createProvider({
    label: 'Membrane Orientation',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_computed_membrane_orientation',
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: MembraneOrientationParams,
    getParams: (data: Structure) => MembraneOrientationParams,
    isApplicable: (data: Structure) => true,
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<MembraneOrientationProps>) => {
        const p = { ...PD.getDefaultValues(MembraneOrientationParams), ...props };
        switch (p.type.name) {
            case 'anvil': return { value: await computeAnvil(ctx, data, p.type.params) };
            case 'opm': return { value: await computeOpm(data) };
        }
    }
});

async function computeAnvil(ctx: CustomProperty.Context, data: Structure, props: ANVILProps): Promise<MembraneOrientation> {
    await AccessibleSurfaceAreaProvider.attach(ctx, data);
    const p = { ...PD.getDefaultValues(MembraneOrientationParams), ...props };
    return await MembraneOrientation.compute(data, p).runInContext(ctx.runtime);
}

async function computeOpm(structure: Structure): Promise<MembraneOrientation> {
    throw Error('TODO impl');
}