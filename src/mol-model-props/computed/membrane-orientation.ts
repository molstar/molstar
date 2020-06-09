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
import { ANVILParams, ANVILProps, computeANVIL } from './membrane-orientation/ANVIL';
import { AccessibleSurfaceAreaProvider } from './accessible-surface-area';
import { MembraneOrientation } from '../../mol-model/structure/model/properties/membrane-orientation';

function getMembraneOrientationParams(data?: Structure) {
    let defaultType = 'anvil' as 'anvil' | 'db'; // TODO flip - db should be default if PDB identifier is known
    return {
        type: PD.MappedStatic(defaultType, {
            'db': PD.EmptyGroup({ label: 'DB' }),
            'anvil': PD.Group(ANVILParams, { label: 'ANVIL' })
        }, { options: [['db', 'DB'], ['anvil', 'ANVIL']] })
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
            case 'db': throw Error('TODO impl'); // TODO RCSB integration
        }
    }
});

async function computeAnvil(ctx: CustomProperty.Context, data: Structure, props: Partial<ANVILProps>): Promise<MembraneOrientation> {
    await AccessibleSurfaceAreaProvider.attach(ctx, data);
    const p = { ...PD.getDefaultValues(ANVILParams), ...props };
    return await computeANVIL(data, p).runInContext(ctx.runtime);
}