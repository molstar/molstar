/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { MolecularSurfaceMeshVisual, MolecularSurfaceMeshParams, StructureMolecularSurfaceMeshVisual } from '../visual/molecular-surface-mesh.js';
import { UnitsRepresentation } from '../units-representation.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ComplexRepresentation, type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { Representation, type RepresentationParamsGetter, type RepresentationContext } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import type { Structure } from '@molstar/model/model/structure';
import { MolecularSurfaceWireframeParams, MolecularSurfaceWireframeVisual, StructureMolecularSurfaceWireframeVisual } from '../visual/molecular-surface-wireframe.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const MolecularSurfaceVisuals = {
    'molecular-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceMeshParams>) => UnitsRepresentation('Molecular surface mesh', ctx, getParams, MolecularSurfaceMeshVisual),
    'structure-molecular-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceMeshParams>) => ComplexRepresentation('Structure Molecular surface mesh', ctx, getParams, StructureMolecularSurfaceMeshVisual),
    'molecular-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceWireframeParams>) => UnitsRepresentation('Molecular surface wireframe', ctx, getParams, MolecularSurfaceWireframeVisual),
    'structure-molecular-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceWireframeParams>) => ComplexRepresentation('Structure Molecular surface mesh', ctx, getParams, StructureMolecularSurfaceWireframeVisual),
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