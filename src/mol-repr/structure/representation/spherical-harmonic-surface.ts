/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 *
 * Spherical Harmonic Surface representation: a smooth shape envelope approximated by
 * a real spherical-harmonic expansion of the radial function r(theta, phi), fit to an
 * atom cloud (each atom pushed out by its size + radiusOffset). The `sphericalHarmonicL`
 * parameter controls the level of detail. At these low degrees the expansion is low-pass,
 * so fitting atoms is equivalent to fitting the molecular surface at far lower cost
 * (Ritchie & Kemp 1999; Duncan & Olson 1993).
 *
 * A single expansion is single-valued in r about one center, so it only
 * represents star-convex shapes. Non-star-shaped inputs (elongated or
 * multi-domain) can be decomposed into several star-shaped lobes via the
 * `maxLobes` parameter, fit separately and blended into one watertight surface;
 * see the visual module for details.
 *
 * Three visual flavours: per-unit (one envelope per chain, instanced across the
 * assembly natively), structure (one merged envelope over the whole input), and
 * assembly (one envelope fit to the asymmetric unit and replicated across the
 * assembly operators).
 */

import { SphericalHarmonicSurfaceMeshVisual, SphericalHarmonicSurfaceMeshParams, StructureSphericalHarmonicSurfaceMeshVisual, AssemblySphericalHarmonicSurfaceMeshVisual } from '../visual/spherical-harmonic-surface-mesh';
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
    'assembly-spherical-harmonic-surface-mesh': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceMeshParams>) => ComplexRepresentation('Assembly spherical harmonic surface mesh', ctx, getParams, AssemblySphericalHarmonicSurfaceMeshVisual),
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
    description: 'Displays a smooth shape envelope approximated by a spherical harmonic expansion of the atom positions (star-convex per lobe; raise maxLobes to split elongated or multi-domain shapes into blended star-shaped lobes).',
    factory: SphericalHarmonicSurfaceRepresentation,
    getParams: getSphericalHarmonicSurfaceParams,
    defaultValues: PD.getDefaultValues(SphericalHarmonicSurfaceParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
});
