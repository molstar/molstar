/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../../mol-model/structure';
import { Task, RuntimeContext } from '../../../../mol-task';
import { getUnitConformationAndRadius, CommonSurfaceProps, ensureReasonableResolution } from './common';
import { PositionData, DensityData, Box3D } from '../../../../mol-math/geometry';
import { MolecularSurfaceCalculationProps, calcMolecularSurface } from '../../../../mol-math/geometry/molecular-surface';
import { OrderedSet } from '../../../../mol-data/int';
import { Boundary } from '../../../../mol-math/geometry/boundary';

export type MolecularSurfaceProps = MolecularSurfaceCalculationProps & CommonSurfaceProps

function getPositionDataAndMaxRadius(structure: Structure, unit: Unit, props: MolecularSurfaceProps) {
    const { probeRadius } = props;
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, props);
    const { indices } = position;
    const n = OrderedSet.size(indices);
    const radii = new Float32Array(OrderedSet.end(indices));

    let maxRadius = 0;
    for (let i = 0; i < n; ++i) {
        const j = OrderedSet.getAt(indices, i);
        const r = radius(j);
        if (maxRadius < r) maxRadius = r;
        radii[j] = r + probeRadius;
    }

    return { position: { ...position, radius: radii }, boundary, maxRadius };
}

export function computeUnitMolecularSurface(structure: Structure, unit: Unit, props: MolecularSurfaceProps) {
    const { box } = unit.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, boundary, maxRadius } = getPositionDataAndMaxRadius(structure, unit, p);
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, boundary, maxRadius, box, p);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, boundary: Boundary, maxRadius: number, box: Box3D | null, props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return calcMolecularSurface(ctx, position, boundary, maxRadius, box, props);
}