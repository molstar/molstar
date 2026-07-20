/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BlobSurfaceMeshParams, BlobSurfaceMeshVisual, StructureBlobSurfaceMeshParams, StructureBlobSurfaceMeshVisual } from '../visual/blob-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';

const BlobSurfaceVisuals = {
    'blob-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BlobSurfaceMeshParams>) => UnitsRepresentation('Blob surface mesh', ctx, getParams, BlobSurfaceMeshVisual),
    'structure-blob-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, StructureBlobSurfaceMeshParams>) => ComplexRepresentation('Structure Blob surface mesh', ctx, getParams, StructureBlobSurfaceMeshVisual),
};

export const BlobSurfaceParams = {
    ...BlobSurfaceMeshParams,
    visuals: PD.MultiSelect(['blob-surface-mesh'], PD.objectToOptions(BlobSurfaceVisuals)),
};
export type BlobSurfaceParams = typeof BlobSurfaceParams
export function getBlobSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return BlobSurfaceParams;
}

export type BlobSurfaceRepresentation = StructureRepresentation<BlobSurfaceParams>
export function BlobSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, BlobSurfaceParams>): BlobSurfaceRepresentation {
    return Representation.createMulti('Blob Surface', ctx, getParams, StructureRepresentationStateBuilder, BlobSurfaceVisuals as unknown as Representation.Def<Structure, BlobSurfaceParams>);
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
    isApplicable: (structure: Structure) => structure.elementCount > 0
});
