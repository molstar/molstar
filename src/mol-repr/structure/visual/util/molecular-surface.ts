/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { Unit } from 'mol-model/structure';
import { Task, RuntimeContext } from 'mol-task';
import { getUnitConformationAndRadius } from './common';
import { PositionData, Box3D, DensityData } from 'mol-math/geometry';

export const MolecularSurfaceCalculationParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 10, step: 0.1 }),
    probeRadius: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
}
export const DefaultMolecularSurfaceCalculationProps = PD.getDefaultValues(MolecularSurfaceCalculationParams)
export type MolecularSurfaceCalculationProps = typeof DefaultMolecularSurfaceCalculationProps

export function computeUnitMolecularSurface(unit: Unit, props: MolecularSurfaceCalculationProps) {
    const { position, radius } = getUnitConformationAndRadius(unit)
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, unit.lookup3d.boundary.box, radius, props);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return {
        transform: Mat4,
        field: Tensor,
        idField: Tensor,
    }
}