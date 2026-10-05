/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import {
  BlobSurfaceMeshParams,
  BlobSurfaceMeshVisual,
  type StructureBlobSurfaceMeshParams,
  StructureBlobSurfaceMeshVisual,
} from '../visual/blob-surface-mesh.js';
import {
  BlobSurfaceWireframeParams,
  BlobSurfaceWireframeVisual,
  type StructureBlobSurfaceWireframeParams,
  StructureBlobSurfaceWireframeVisual,
} from '../visual/blob-surface-wireframe.js';
import { UnitsRepresentation } from '../units-representation.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  ComplexRepresentation,
  type StructureRepresentation,
  StructureRepresentationProvider,
  StructureRepresentationStateBuilder,
} from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const BlobSurfaceVisuals = {
  'blob-surface-mesh': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, BlobSurfaceMeshParams>,
  ) => UnitsRepresentation('Blob surface mesh', ctx, getParams, BlobSurfaceMeshVisual),
  'structure-blob-surface-mesh': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, StructureBlobSurfaceMeshParams>,
  ) => ComplexRepresentation('Structure Blob surface mesh', ctx, getParams, StructureBlobSurfaceMeshVisual),
  'blob-surface-wireframe': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, BlobSurfaceWireframeParams>,
  ) => UnitsRepresentation('Blob surface wireframe', ctx, getParams, BlobSurfaceWireframeVisual),
  'structure-blob-surface-wireframe': (
    ctx: RepresentationContext,
    getParams: RepresentationParamsGetter<Structure, StructureBlobSurfaceWireframeParams>,
  ) => ComplexRepresentation('Structure Blob surface wireframe', ctx, getParams, StructureBlobSurfaceWireframeVisual),
};

export const BlobSurfaceParams = {
  ...BlobSurfaceMeshParams,
  ...BlobSurfaceWireframeParams,
  visuals: PD.MultiSelect(['blob-surface-mesh'], PD.objectToOptions(BlobSurfaceVisuals)),
  solidInterior: PD.Boolean(true, {
    ...BaseGeometry.ShadingCategory,
    description: 'Render a solid cap where the camera near plane or a clip object cuts a closed surface',
  }),
};
export type BlobSurfaceParams = typeof BlobSurfaceParams;
export function getBlobSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
  return BlobSurfaceParams;
}

export type BlobSurfaceRepresentation = StructureRepresentation<BlobSurfaceParams>;
export function BlobSurfaceRepresentation(
  ctx: RepresentationContext,
  getParams: RepresentationParamsGetter<Structure, BlobSurfaceParams>,
): BlobSurfaceRepresentation {
  return Representation.createMulti(
    'Blob Surface',
    ctx,
    getParams,
    StructureRepresentationStateBuilder,
    BlobSurfaceVisuals as unknown as Representation.Def<Structure, BlobSurfaceParams>,
  );
}

export const BlobSurfaceRepresentationProvider = StructureRepresentationProvider({
  name: 'blob-surface',
  label: 'Blob Surface',
  description: 'Displays an approximate blobby surface useful for simplified depictions.',
  factory: BlobSurfaceRepresentation,
  getParams: getBlobSurfaceParams,
  defaultValues: PD.getDefaultValues(BlobSurfaceParams),
  defaultColorTheme: { name: 'chain-id' },
  defaultSizeTheme: { name: 'physical' },
  isApplicable: (structure: Structure) => structure.elementCount > 0,
});
