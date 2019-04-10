/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from 'mol-model/structure';
import { Task, RuntimeContext } from 'mol-task';
import { getUnitConformationAndRadius } from './common';
import { PositionData, DensityData } from 'mol-math/geometry';
import { MolecularSurfaceCalculationProps, calcMolecularSurface } from 'mol-math/geometry/molecular-surface';
import { OrderedSet } from 'mol-data/int';

export function computeUnitMolecularSurface(unit: Unit, props: MolecularSurfaceCalculationProps) {
    const { position, radius } = getUnitConformationAndRadius(unit)

    return Task.create('Molecular Surface', async ctx => {
        const { indices } = position
        const n = OrderedSet.size(indices)
        const radii = new Float32Array(OrderedSet.max(indices))

        let maxRadius = 0
        for (let i = 0; i < n; ++i) {
            const j = OrderedSet.getAt(indices, i)
            const r = radius(j)
            if (maxRadius < r) maxRadius = r
            radii[j] = r + props.probeRadius

            if (i % 100000 === 0 && ctx.shouldUpdate) {
                await ctx.update({ message: 'calculating max radius', current: i, max: n })
            }
        }

        return await MolecularSurface(ctx, { ...position, radius: radii }, maxRadius, props);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number,  props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return calcMolecularSurface(ctx, position, maxRadius, props)
}