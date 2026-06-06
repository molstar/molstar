/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Spherical Harmonic Surface representation: a molecular surface approximated by
 * a real spherical-harmonic expansion of the radial function r(theta, phi). The
 * `sphericalHarmonicL` parameter controls the level of detail.
 */

import { SphericalHarmonicSurfaceMeshVisual, SphericalHarmonicSurfaceMeshParams, StructureSphericalHarmonicSurfaceMeshVisual, ProtomerSphericalHarmonicSurfaceMeshVisual, ResidueSphericalHarmonicSurfaceMeshVisual } from '../visual/spherical-harmonic-surface-mesh';
import { UnitsRepresentation } from '../units-representation';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../../mol-repr/representation';
import { ThemeRegistryContext } from '../../../mol-theme/theme';
import { Structure } from '../../../mol-model/structure';
import { BaseGeometry } from '../../../mol-geo/geometry/base';

const SphericalHarmonicSurfaceVisuals = {
    'spherical-harmonic-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceMeshParams>) => UnitsRepresentation('Spherical harmonic surface mesh', ctx, getParams, SphericalHarmonicSurfaceMeshVisual),
    'structure-spherical-harmonic-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceMeshParams>) => ComplexRepresentation('Structure spherical harmonic surface mesh', ctx, getParams, StructureSphericalHarmonicSurfaceMeshVisual),
    'protomer-spherical-harmonic-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceMeshParams>) => ComplexRepresentation('Protomer spherical harmonic surface mesh', ctx, getParams, ProtomerSphericalHarmonicSurfaceMeshVisual),
    'residue-spherical-harmonic-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceMeshParams>) => ComplexRepresentation('Residue spherical harmonic surface mesh', ctx, getParams, ResidueSphericalHarmonicSurfaceMeshVisual),
};

export const SphericalHarmonicSurfaceParams = {
    ...SphericalHarmonicSurfaceMeshParams,
    visuals: PD.MultiSelect(['spherical-harmonic-surface-mesh'], PD.objectToOptions(SphericalHarmonicSurfaceVisuals)),
    bumpFrequency: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
};
export type SphericalHarmonicSurfaceParams = typeof SphericalHarmonicSurfaceParams
export function getSphericalHarmonicSurfaceParams(ctx: ThemeRegistryContext, structure: Structure) {
    return SphericalHarmonicSurfaceParams;
}

export type SphericalHarmonicSurfaceRepresentation = StructureRepresentation<SphericalHarmonicSurfaceParams>
export function SphericalHarmonicSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceParams>): SphericalHarmonicSurfaceRepresentation {
    return Representation.createMulti('Spherical Harmonic Surface', ctx, getParams, StructureRepresentationStateBuilder, SphericalHarmonicSurfaceVisuals as unknown as Representation.Def<Structure, SphericalHarmonicSurfaceParams>);
}

export const SphericalHarmonicSurfaceRepresentationProvider = StructureRepresentationProvider({
    name: 'spherical-harmonic-surface',
    label: 'Spherical Harmonic Surface',
    description: 'Displays a molecular surface approximated by a spherical harmonic expansion.',
    factory: SphericalHarmonicSurfaceRepresentation,
    getParams: getSphericalHarmonicSurfaceParams,
    defaultValues: PD.getDefaultValues(SphericalHarmonicSurfaceParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
});
