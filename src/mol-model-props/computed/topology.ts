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
import { ANVILParams, Topology } from './topology/ANVIL';

export const TopologyParams = {
    ...ANVILParams
};
export type TopologyParams = typeof TopologyParams
export type TopologyProps = PD.Values<TopologyParams>

export type TopologyValue = Map<number, Topology>

export const TopologyProvider: CustomStructureProperty.Provider<TopologyParams, Topology> = CustomStructureProperty.createProvider({
    label: 'Predicted Membrane Topology',
    descriptor: CustomPropertyDescriptor({
        name: 'molstar_topology',
        // TODO `cifExport`
    }),
    type: 'root',
    defaultParams: TopologyParams,
    getParams: (data: Structure) => TopologyParams,
    isApplicable: (data: Structure) => true, 
    // TODO needs ASA to be computed (or 'resolved' before trying computing topology) - how to achieve?
    // TODO potentially, this could behave like secondary structure info where data can be either parsed or computed
    obtain: async (ctx: CustomProperty.Context, data: Structure, props: Partial<TopologyProps>) => {
        const p = { ...PD.getDefaultValues(TopologyParams), ...props };
        return { value: await Topology.compute(data, p).runInContext(ctx.runtime) };
    }
});