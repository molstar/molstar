/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Unit, Structure } from '@molstar/model/model/structure';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { getUnitConformationAndRadius, type CommonSurfaceProps, ensureReasonableResolution, getStructureConformationAndRadius } from './common.js';
import type { PositionData, DensityData, Box3D } from '@molstar/core/math/geometry';
import { MolecularSurfaceCalculationParams, type MolecularSurfaceCalculationProps, calcMolecularSurface } from '@molstar/core/math/geometry/molecular-surface';
import { OrderedSet } from '@molstar/core/data/int';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export const CommonMolecularSurfaceCalculationParams = {
    ...MolecularSurfaceCalculationParams,
    resolution: { ...MolecularSurfaceCalculationParams.resolution, ...BaseGeometry.CustomQualityParamInfo },
    probePositions: { ...MolecularSurfaceCalculationParams.probePositions, ...BaseGeometry.CustomQualityParamInfo },
    floodfill: PD.Select('off', PD.arrayToOptions(['off', 'inside', 'outside']), { description: 'If and how to floodfill the molecular surface.' }),
};

export type MolecularSurfaceProps = MolecularSurfaceCalculationProps & CommonSurfaceProps

function getUnitPositionDataAndMaxRadius(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: MolecularSurfaceProps) {
    const { probeRadius } = props;
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
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

export function computeUnitMolecularSurface(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: MolecularSurfaceProps) {
    const { position, boundary, maxRadius } = getUnitPositionDataAndMaxRadius(structure, unit, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, boundary, maxRadius, boundary.box, p);
    });
}

//

function getStructurePositionDataAndMaxRadius(structure: Structure, sizeTheme: SizeTheme<any>, props: MolecularSurfaceProps) {
    const { probeRadius } = props;
    const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
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

export function computeStructureMolecularSurface(structure: Structure, sizeTheme: SizeTheme<any>, props: MolecularSurfaceProps) {
    const { position, boundary, maxRadius } = getStructurePositionDataAndMaxRadius(structure, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Molecular Surface', async ctx => {
        return await MolecularSurface(ctx, position, boundary, maxRadius, boundary.box, p);
    });
}

//

async function MolecularSurface(ctx: RuntimeContext, position: Required<PositionData>, boundary: Boundary, maxRadius: number, box: Box3D | null, props: MolecularSurfaceCalculationProps): Promise<DensityData> {
    return calcMolecularSurface(ctx, position, boundary, maxRadius, box, props);
}