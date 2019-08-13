/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit } from '../../../../mol-model/structure';
import { Task, RuntimeContext } from '../../../../mol-task';
import { getUnitConformationAndRadius } from './common';
import { PositionData, DensityData } from '../../../../mol-math/geometry';
import { MolecularSurfaceCalculationProps, calcMolecularSurface } from '../../../../mol-math/geometry/molecular-surface';
import { OrderedSet } from '../../../../mol-data/int';

export type MolecularSurfaceProps = MolecularSurfaceCalculationProps & {
    ignoreHydrogens: boolean
}

function getPositionDataAndMaxRadius(unit: Unit, props: MolecularSurfaceProps) {
    const { probeRadius, ignoreHydrogens } = props
    const { position, radius } = getUnitConformationAndRadius(unit, ignoreHydrogens)
    const { indices } = position
    const n = OrderedSet.size(indices)
    const radii = new Float32Array(OrderedSet.end(indices))

    let maxRadius = 0
    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i)
        const r = radius(j)
        if (maxRadius < r) maxRadius = r
        radii[j] = r + probeRadius
    }

    return { position: { ...position, radius: radii }, maxRadius }
}

export function computeUnitMolecularSurface(unit: Unit, props: MolecularSurfaceProps) {
    const { position, maxRadius } = getPositionDataAndMaxRadius(unit, props)
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, maxRadius, props);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number,  props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return calcMolecularSurface(ctx, position, maxRadius, props)
}