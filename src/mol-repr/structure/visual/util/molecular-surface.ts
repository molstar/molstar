/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../../mol-model/structure';
import { Task, RuntimeContext } from '../../../../mol-task';
import { getUnitConformationAndRadius, CommonSurfaceProps } from './common';
import { PositionData, DensityData, Box3D } from '../../../../mol-math/geometry';
import { MolecularSurfaceCalculationProps, calcMolecularSurface } from '../../../../mol-math/geometry/molecular-surface';
import { OrderedSet } from '../../../../mol-data/int';

export type MolecularSurfaceProps = MolecularSurfaceCalculationProps & CommonSurfaceProps

function getPositionDataAndMaxRadius(structure: Structure, unit: Unit, props: MolecularSurfaceProps) {
    const { probeRadius } = props
    const { position, radius } = getUnitConformationAndRadius(structure, unit, props)
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

export function computeUnitMolecularSurface(structure: Structure, unit: Unit, props: MolecularSurfaceProps) {
    const { box } = unit.lookup3d.boundary
    const { position, maxRadius } = getPositionDataAndMaxRadius(structure, unit, props)
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, maxRadius, box, props);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, maxRadius: number, box: Box3D | null, props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return calcMolecularSurface(ctx, position, maxRadius, box, props)
}