/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MolecularSurfaceMeshVisual, MolecularSurfaceMeshParams } from '../visual/molecular-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { MolecularSurfaceWireframeParams, MolecularSurfaceWireframeVisual } from '../visual/molecular-surface-wireframe';

const MolecularSurfaceVisuals = {
    'molecular-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceMeshParams>) => UnitsRepresentation('Molecular surface mesh', ctx, getParams, MolecularSurfaceMeshVisual),
    'molecular-surface-wireframe': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, MolecularSurfaceWireframeParams>) => UnitsRepresentation('Molecular surface wireframe', ctx, getParams, MolecularSurfaceWireframeVisual),
};

export const MolecularSurfaceParams = {
    ...MolecularSurfaceMeshParams,
    ...MolecularSurfaceWireframeParams,
    visuals: PD.MultiSelect(['molecular-surface-mesh'], PD.objectToOptions(MolecularSurfaceVisuals)),
};
export type MolecularSurfaceParams = typeof MolecularSurfaceParams
export function getMolecularSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(MolecularSurfaceParams);
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
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0
});