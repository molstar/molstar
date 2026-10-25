/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  UnitsMeshParams,
  UnitsSpheresParams,
  type UnitsVisual,
  UnitsSpheresVisual,
  UnitsMeshVisual,
} from '../units-visual.js';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import {
  createElementSphereImpostor,
  ElementIterator,
  getElementLoci,
  eachElement,
  createElementSphereMesh,
  createStructureElementSphereImpostor,
  getSerialElementLoci,
  eachSerialElement,
  createStructureElementSphereMesh,
} from './util/element.js';
import type { VisualUpdateState } from '../../util.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { checkSphereImpostorSupport } from './util/common.js';
import {
  ComplexMeshParams,
  ComplexMeshVisual,
  ComplexSpheresParams,
  ComplexSpheresVisual,
  type ComplexVisual,
} from '../complex-visual.js';
import { chooseVisual } from '../../visual-factory.js';

export const CommonElementSphereParams = {
  sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
  detail: PD.Numeric(0, { min: 0, max: 3, step: 1 }, BaseGeometry.CustomQualityParamInfo),
  ignoreHydrogens: PD.Boolean(false),
  ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
  traceOnly: PD.Boolean(false),
  tryUseImpostor: PD.Boolean(true),
  stride: PD.Numeric(1, { min: 1, max: 100, step: 1 }),
};

//

export const ElementSphereParams = {
  ...UnitsMeshParams,
  ...UnitsSpheresParams,
  ...CommonElementSphereParams,
};
export type ElementSphereParams = typeof ElementSphereParams;

export const ElementSphereVisual = chooseVisual(
  { impostor: ElementSphereImpostorVisual, mesh: ElementSphereMeshVisual },
  (props: PD.Values<ElementSphereParams>, graphics?: GraphicsContext) =>
    props.tryUseImpostor && checkSphereImpostorSupport(graphics) ? 'impostor' : 'mesh',
);

export function ElementSphereImpostorVisual(materialId: number): UnitsVisual<ElementSphereParams> {
  return UnitsSpheresVisual<ElementSphereParams>(
    {
      defaultProps: PD.getDefaultValues(ElementSphereParams),
      createGeometry: createElementSphereImpostor,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<ElementSphereParams>,
        currentProps: PD.Values<ElementSphereParams>,
      ) => {
        state.createGeometry =
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.stride !== currentProps.stride;
      },
    },
    materialId,
  );
}

export function ElementSphereMeshVisual(materialId: number): UnitsVisual<ElementSphereParams> {
  return UnitsMeshVisual<ElementSphereParams>(
    {
      defaultProps: PD.getDefaultValues(ElementSphereParams),
      createGeometry: createElementSphereMesh,
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<ElementSphereParams>,
        currentProps: PD.Values<ElementSphereParams>,
      ) => {
        state.createGeometry =
          newProps.sizeFactor !== currentProps.sizeFactor ||
          newProps.detail !== currentProps.detail ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.stride !== currentProps.stride;
      },
    },
    materialId,
  );
}

//

export const StructureElementSphereParams = {
  ...ComplexMeshParams,
  ...ComplexSpheresParams,
  ...CommonElementSphereParams,
};
export type StructureElementSphereParams = typeof ElementSphereParams;

export const StructureElementSphereVisual = chooseVisual(
  { impostor: StructureElementSphereImpostorVisual, mesh: StructureElementSphereMeshVisual },
  (props: PD.Values<ElementSphereParams>, graphics?: GraphicsContext) =>
    props.tryUseImpostor && graphics && graphics.capabilities.fragDepth ? 'impostor' : 'mesh',
);

export function StructureElementSphereImpostorVisual(materialId: number): ComplexVisual<StructureElementSphereParams> {
  return ComplexSpheresVisual<StructureElementSphereParams>(
    {
      defaultProps: PD.getDefaultValues(StructureElementSphereParams),
      createGeometry: createStructureElementSphereImpostor,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<StructureElementSphereParams>,
        currentProps: PD.Values<StructureElementSphereParams>,
      ) => {
        state.createGeometry =
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.stride !== currentProps.stride;
      },
    },
    materialId,
  );
}

export function StructureElementSphereMeshVisual(materialId: number): ComplexVisual<StructureElementSphereParams> {
  return ComplexMeshVisual<StructureElementSphereParams>(
    {
      defaultProps: PD.getDefaultValues(StructureElementSphereParams),
      createGeometry: createStructureElementSphereMesh,
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<StructureElementSphereParams>,
        currentProps: PD.Values<StructureElementSphereParams>,
      ) => {
        state.createGeometry =
          newProps.sizeFactor !== currentProps.sizeFactor ||
          newProps.detail !== currentProps.detail ||
          newProps.ignoreHydrogens !== currentProps.ignoreHydrogens ||
          newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant ||
          newProps.traceOnly !== currentProps.traceOnly ||
          newProps.stride !== currentProps.stride;
      },
    },
    materialId,
  );
}
