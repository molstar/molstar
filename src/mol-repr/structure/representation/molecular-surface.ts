/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolecularSurfaceMeshVisual, MolecularSurfaceMeshParams, StructureMolecularSurfaceMeshVisual } from '../visual/molecular-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { MolecularSurfaceWireframeParams, MolecularSurfaceWireframeVisual } from '../visual/molecular-surface-wireframe';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

const MolecularSurfaceVisuals = {
    'molecular-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceMeshParams>) => UnitsRepresentation('Molecular surface mesh', ctx, getParams, MolecularSurfaceMeshVisual),
    'structure-molecular-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceMeshParams>) => ComplexRepresentation('Structure Molecular surface mesh', ctx, getParams, StructureMolecularSurfaceMeshVisual),
    'molecular-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceWireframeParams>) => UnitsRepresentation('Molecular surface wireframe', ctx, getParams, MolecularSurfaceWireframeVisual),
};

export const MolecularSurfaceParams = {
    ...MolecularSurfaceMeshParams,
    ...MolecularSurfaceWireframeParams,
    visuals: PD.MultiSelect(['molecular-surface-mesh'], PD.objectToOptions(MolecularSurfaceVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type MolecularSurfaceParams = typeof MolecularSurfaceParams
export function getMolecularSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return MolecularSurfaceParams;
}

export type MolecularSurfaceRepresentation = StructureRepresentation<MolecularSurfaceParams>
export function MolecularSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceParams>): MolecularSurfaceRepresentation {
    return Representation.createMulti('Molecular Surface', ctx, getParams, StructureRepresentationStateBuilder, MolecularSurfaceVisuals as unknown as Representation.Def<Structure, MolecularSurfaceParams>);
}

export const MolecularSurfaceRepresentationProvider = StructureRepresentationProvider({
    name: 'molecular-surface',
    label: 'Molecular Surface',
    description: 'Displays a molecular surface.',
    factory: MolecularSurfaceRepresentation,
    getParams: getMolecularSurfaceParams,
    defaultValues: PD.getDefaultValues(MolecularSurfaceParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});