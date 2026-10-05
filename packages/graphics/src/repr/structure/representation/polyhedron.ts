/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import type { Structure } from '@molstar/model/model/structure';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ComplexRepresentation } from '../complex-representation.js';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { CoordinationPolyhedronMeshParams, CoordinationPolyhedronMeshVisual } from '../visual/coordination-polyhedron-mesh.js';

const PolyhedronVisuals = {
    'coordination-polyhedron-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CoordinationPolyhedronMeshParams>) => ComplexRepresentation('Coordination Polyhedron mesh', ctx, getParams, CoordinationPolyhedronMeshVisual),
};

export const PolyhedronParams = {
    ...CoordinationPolyhedronMeshParams,
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
    visuals: PD.MultiSelect(['coordination-polyhedron-mesh'], PD.objectToOptions(PolyhedronVisuals)),
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
