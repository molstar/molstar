/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { UnitsMeshParams, type UnitsVisual, UnitsMeshVisual } from '../units-visual.js';
import type { VisualContext } from '../../visual.js';
import type { Unit, Structure } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import {
  CommonMolecularSurfaceCalculationParams,
  computeStructureMolecularSurface,
  computeUnitMolecularSurface,
} from './util/molecular-surface.js';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import {
  ElementIterator,
  getElementLoci,
  eachElement,
  getSerialElementLoci,
  eachSerialElement,
} from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import { CommonSurfaceParams } from './util/common.js';
import { Sphere3D } from '@molstar/core/math/geometry';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { applyMeshColorSmoothing } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { ColorSmoothingParams, getColorSmoothingProps } from '@molstar/graphics/geo/geometry/base';
import { ValueCell } from '@molstar/core/util';
import { ComplexMeshVisual, type ComplexVisual } from '../complex-visual.js';
import { Tensor } from '@molstar/core/math/linear-algebra/tensor';

export const MolecularSurfaceMeshParams = {
  ...UnitsMeshParams,
  ...CommonMolecularSurfaceCalculationParams,
  ...CommonSurfaceParams,
  ...ColorSmoothingParams,
};
export type MolecularSurfaceMeshParams = typeof MolecularSurfaceMeshParams;
export type MolecularSurfaceMeshProps = PD.Values<MolecularSurfaceMeshParams>;

type MolecularSurfaceMeta = {
  resolution?: number;
  colorTexture?: Texture;
};

//

async function createMolecularSurfaceMesh(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: MolecularSurfaceMeshProps,
  mesh?: Mesh,
): Promise<Mesh> {
  const { transform, field, idField, resolution, maxRadius } = await computeUnitMolecularSurface(
    structure,
    unit,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const params = {
    isoLevel: props.probeRadius,
    scalarField:
      props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
    idField,
  };
  const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

  if (props.includeParent) {
    const iterations = Math.ceil(2 / props.resolution);
    Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
  }

  Mesh.transform(surface, transform);
  if (ctx.webgl && !ctx.webgl.isWebGL2) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, maxRadius);
  surface.setBoundingSphere(sphere);
  (surface.meta as MolecularSurfaceMeta).resolution = resolution;

  return surface;
}

export function MolecularSurfaceMeshVisual(materialId: number): UnitsVisual<MolecularSurfaceMeshParams> {
  return UnitsMeshVisual<MolecularSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(MolecularSurfaceMeshParams),
      createGeometry: createMolecularSurfaceMesh,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<MolecularSurfaceMeshParams>,
        currentProps: PD.Values<MolecularSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.probeRadius !== currentProps.probeRadius ||
          newProps.probePositions !== currentProps.probePositions ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.includeParent !== currentProps.includeParent ||
          newProps.floodfill !== currentProps.floodfill;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: MeshValues,
        geometry: Mesh,
        props: PD.Values<MolecularSurfaceMeshParams>,
        theme: Theme,
        webgl?: WebGLContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as MolecularSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp) {
          applyMeshColorSmoothing(values, csp, webgl, colorTexture);
          (geometry.meta as MolecularSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: Mesh) => {
        (geometry.meta as MolecularSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}

//

async function createStructureMolecularSurfaceMesh(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: MolecularSurfaceMeshProps,
  mesh?: Mesh,
): Promise<Mesh> {
  const { transform, field, idField, resolution, maxRadius } = await computeStructureMolecularSurface(
    structure,
    theme.size,
    props,
  ).runInContext(ctx.runtime);

  const params = {
    isoLevel: props.probeRadius,
    scalarField:
      props.floodfill !== 'off' ? Tensor.createFloodfilled(field, props.probeRadius, props.floodfill) : field,
    idField,
  };
  const surface = await computeMarchingCubesMesh(params, mesh).runAsChild(ctx.runtime);

  if (props.includeParent) {
    const iterations = Math.ceil(2 / props.resolution);
    Mesh.smoothEdges(surface, { iterations, maxNewEdgeLength: Math.sqrt(2) });
  }

  Mesh.transform(surface, transform);
  if (ctx.webgl && !ctx.webgl.isWebGL2) {
    Mesh.uniformTriangleGroup(surface);
    ValueCell.updateIfChanged(surface.varyingGroup, false);
  } else {
    ValueCell.updateIfChanged(surface.varyingGroup, true);
  }

  const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius);
  surface.setBoundingSphere(sphere);
  (surface.meta as MolecularSurfaceMeta).resolution = resolution;

  return surface;
}

export function StructureMolecularSurfaceMeshVisual(materialId: number): ComplexVisual<MolecularSurfaceMeshParams> {
  return ComplexMeshVisual<MolecularSurfaceMeshParams>(
    {
      defaultProps: PD.getDefaultValues(MolecularSurfaceMeshParams),
      createGeometry: createStructureMolecularSurfaceMesh,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<MolecularSurfaceMeshParams>,
        currentProps: PD.Values<MolecularSurfaceMeshParams>,
      ) => {
        state.createGeometry =
          newProps.resolution !== currentProps.resolution ||
          newProps.probeRadius !== currentProps.probeRadius ||
          newProps.probePositions !== currentProps.probePositions ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.includeParent !== currentProps.includeParent ||
          newProps.floodfill !== currentProps.floodfill;

        if (newProps.smoothColors.name !== currentProps.smoothColors.name) {
          state.updateColor = true;
        } else if (newProps.smoothColors.name === 'on' && currentProps.smoothColors.name === 'on') {
          if (newProps.smoothColors.params.resolutionFactor !== currentProps.smoothColors.params.resolutionFactor)
            state.updateColor = true;
          if (newProps.smoothColors.params.sampleStride !== currentProps.smoothColors.params.sampleStride)
            state.updateColor = true;
        }
      },
      processValues: (
        values: MeshValues,
        geometry: Mesh,
        props: PD.Values<MolecularSurfaceMeshParams>,
        theme: Theme,
        webgl?: WebGLContext,
      ) => {
        const { resolution, colorTexture } = geometry.meta as MolecularSurfaceMeta;
        const csp = getColorSmoothingProps(props.smoothColors, theme.color.preferSmoothing, resolution);
        if (csp) {
          applyMeshColorSmoothing(values, csp, webgl, colorTexture);
          (geometry.meta as MolecularSurfaceMeta).colorTexture = values.tColorGrid.ref.value;
        }
      },
      dispose: (geometry: Mesh) => {
        (geometry.meta as MolecularSurfaceMeta).colorTexture?.destroy();
      },
    },
    materialId,
  );
}
