/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsMeshParams, type UnitsVisual, UnitsMeshVisual } from '../units-visual.js';
import { type ComplexVisual, ComplexMeshParams, ComplexMeshVisual } from '../complex-visual.js';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import {
  ElementIterator,
  getElementLoci,
  eachElement,
  getSerialElementLoci,
  eachSerialElement,
} from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import {
  BlobDensityParams,
  computeUnitBlobSurface,
  computeStructureBlobSurface,
  shouldUpdateBlobGeometry,
} from './util/blob-surface.js';
import { Sphere3D } from '@molstar/core/math/geometry';
import { ValueCell } from '@molstar/core/util/value-cell';

export const BlobSurfaceMeshParams = {
  ...UnitsMeshParams,
  ...BlobDensityParams,
};
export type BlobSurfaceMeshParams = typeof BlobSurfaceMeshParams;

export const StructureBlobSurfaceMeshParams = {
  ...ComplexMeshParams,
  ...BlobDensityParams,
};
export type StructureBlobSurfaceMeshParams = typeof StructureBlobSurfaceMeshParams;

//

async function createBlobSurfaceMesh(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: PD.Values<BlobSurfaceMeshParams>,
  mesh?: Mesh,
): Promise<Mesh> {
  const { smoothness, radiusOffset } = props;
  const { transform, field, idField, radiusFactor, maxRadius } = await computeUnitBlobSurface(
    structure,
    unit,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const isoLevel = Math.exp(-smoothness) / radiusFactor;
  const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, mesh).runAsChild(
    ctx.runtime,
  );

  Mesh.transform(surface, transform);
  if (ctx.webgl && !ctx.webgl.isWebGL2) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
  const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius + extraRadius);
  surface.setBoundingSphere(sphere);

  return surface;
}

export function BlobSurfaceMeshVisual(materialId: number): UnitsVisual<BlobSurfaceMeshParams> {
  return UnitsMeshVisual<BlobSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(BlobSurfaceMeshParams),
      createGeometry: createBlobSurfaceMesh,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<BlobSurfaceMeshParams>,
        currentProps: PD.Values<BlobSurfaceMeshParams>,
      ) => {
        state.createGeometry = shouldUpdateBlobGeometry(newProps, currentProps);
      },
    },
    materialId,
  );
}

//

async function createStructureBlobSurfaceMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: PD.Values<StructureBlobSurfaceMeshParams>,
  mesh?: Mesh,
): Promise<Mesh> {
  const { smoothness, radiusOffset } = props;
  const { transform, field, idField, radiusFactor, maxRadius } = await computeStructureBlobSurface(
    structure,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const isoLevel = Math.exp(-smoothness) / radiusFactor;
  const surface = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }, mesh).runAsChild(
    ctx.runtime,
  );

  Mesh.transform(surface, transform);
  if (ctx.webgl && !ctx.webgl.isWebGL2) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const extraRadius = radiusOffset * (1 + Math.exp(-smoothness));
  const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius);
  surface.setBoundingSphere(sphere);

  return surface;
}

export function StructureBlobSurfaceMeshVisual(materialId: number): ComplexVisual<StructureBlobSurfaceMeshParams> {
  return ComplexMeshVisual<StructureBlobSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(StructureBlobSurfaceMeshParams),
      createGeometry: createStructureBlobSurfaceMesh,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<StructureBlobSurfaceMeshParams>,
        currentProps: PD.Values<StructureBlobSurfaceMeshParams>,
      ) => {
        state.createGeometry = shouldUpdateBlobGeometry(newProps, currentProps);
      },
    },
    materialId,
  );
}
