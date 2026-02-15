/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BaseGeometry } from '../../../mol-geo/geometry/base';
import { Structure } from '../../../mol-model/structure';
import { Representation, RepresentationContext, RepresentationParamsGetter } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { PolyhedronMeshParams, PolyhedronMeshVisual } from '../visual/polyhedron-mesh';

const PolyhedronVisuals = {
    'polyhedron-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolyhedronMeshParams>) => ComplexRepresentation('Polyhedron mesh', ctx, getParams, PolyhedronMeshVisual),
};

export const PolyhedronParams = {
    ...PolyhedronMeshParams,
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
    visuals: PD.MultiSelect(['polyhedron-mesh'], PD.objectToOptions(PolyhedronVisuals)),
};
export type PolyhedronParams = typeof PolyhedronParams
export function getPolyhedronParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PolyhedronParams;
}

export type PolyhedronRepresentation = StructureRepresentation<PolyhedronParams>
export function PolyhedronRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolyhedronParams>): PolyhedronRepresentation {
    return Representation.createMulti('Polyhedron', ctx, getParams, StructureRepresentationStateBuilder, PolyhedronVisuals as unknown as Representation.Def<Structure, PolyhedronParams>);
}

export const PolyhedronRepresentationProvider = StructureRepresentationProvider({
    name: 'polyhedron',
    label: 'Polyhedron',
    description: 'Displays coordination polyhedra around atoms with enough bonds.',
    factory: PolyhedronRepresentation,
    getParams: getPolyhedronParams,
    defaultValues: PD.getDefaultValues(PolyhedronParams),
    defaultColorTheme: { name: 'element-symbol' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
    getData: (structure: Structure, props: PD.Values<PolyhedronParams>) => {
        return props.includeParent ? structure.asParent() : structure;
    },
    mustRecreate: (oldProps: PD.Values<PolyhedronParams>, newProps: PD.Values<PolyhedronParams>) => {
        return oldProps.includeParent !== newProps.includeParent;
    }
});
