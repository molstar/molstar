/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Unit, Structure } from '@molstar/model/model/structure';
import { Task, type RuntimeContext } from '@molstar/core/task';
import {
  getUnitConformationAndRadius,
  type CommonSurfaceProps,
  ensureReasonableResolution,
  getStructureConformationAndRadius,
} from './common.js';
import type { PositionData, DensityData, Box3D } from '@molstar/core/math/geometry';
import {
  MolecularSurfaceCalculationParams,
  type MolecularSurfaceCalculationProps,
  calcMolecularSurface,
} from '@molstar/core/math/geometry/molecular-surface';
import { Interval, OrderedSet } from '@molstar/core/data/int';
import type { Boundary } from '@molstar/core/math/geometry/boundary';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

export const CommonMolecularSurfaceCalculationParams = {
  ...MolecularSurfaceCalculationParams,
  resolution: { ...MolecularSurfaceCalculationParams.resolution, ...BaseGeometry.CustomQualityParamInfo },
  probePositions: { ...MolecularSurfaceCalculationParams.probePositions, ...BaseGeometry.CustomQualityParamInfo },
  floodfill: PD.Select('off', PD.arrayToOptions(['off', 'inside', 'outside']), {
    description: 'If and how to floodfill the molecular surface.',
  }),
};

export type MolecularSurfaceProps = MolecularSurfaceCalculationProps & CommonSurfaceProps;

function getUnitPositionDataAndMaxRadius(
  structure: Structure,
  unit: Unit,
  sizeTheme: SizeTheme<any>,
  props: MolecularSurfaceProps,
) {
  const { probeRadius } = props;
  const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
  const { indices, x, y, z, id } = position;
  const n = OrderedSet.size(indices);
  // compact per unit: the indices are the model's element indices, so arrays indexed by them would be sized to the
  // model, per unit (perf-audit-2 07/P6). Float64 keeps the coordinates exact whatever their array.
  const cx = new Float64Array(n);
  const cy = new Float64Array(n);
  const cz = new Float64Array(n);
  const radii = new Float32Array(n);

  let maxRadius = 0;
  for (let i = 0; i < n; ++i) {
    const j = OrderedSet.getAt(indices, i);
    const r = radius(j);
    if (maxRadius < r) maxRadius = r;
    cx[i] = x[j];
    cy[i] = y[j];
    cz[i] = z[j];
    radii[i] = r + probeRadius;
  }

  const compact = { indices: Interval.ofLength(n), x: cx, y: cy, z: cz, id, radius: radii };
  return { position: compact, boundary, maxRadius };
}

export function computeUnitMolecularSurface(
  structure: Structure,
  unit: Unit,
  sizeTheme: SizeTheme<any>,
  props: MolecularSurfaceProps,
) {
  const { position, boundary, maxRadius } = getUnitPositionDataAndMaxRadius(structure, unit, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props);
  return Task.create('Molecular Surface', async (ctx) => {
    return await MolecularSurface(ctx, position, boundary, maxRadius, boundary.box, p);
  });
}

//

function getStructurePositionDataAndMaxRadius(
  structure: Structure,
  sizeTheme: SizeTheme<any>,
  props: MolecularSurfaceProps,
) {
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

export function computeStructureMolecularSurface(
  structure: Structure,
  sizeTheme: SizeTheme<any>,
  props: MolecularSurfaceProps,
) {
  const { position, boundary, maxRadius } = getStructurePositionDataAndMaxRadius(structure, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props);
  return Task.create('Molecular Surface', async (ctx) => {
    return await MolecularSurface(ctx, position, boundary, maxRadius, boundary.box, p);
  });
}

//

async function MolecularSurface(
  ctx: RuntimeContext,
  position: Required<PositionData>,
  boundary: Boundary,
  maxRadius: number,
  box: Box3D | null,
  props: MolecularSurfaceCalculationProps,
): Promise<DensityData> {
  return calcMolecularSurface(ctx, position, boundary, maxRadius, box, props);
}
