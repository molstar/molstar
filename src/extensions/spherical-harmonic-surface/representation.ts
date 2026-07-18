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
 * multi-domain, e.g. a threaded rRNA chain) can be decomposed into several
 * star-shaped lobes via the `lobes` parameter - either by contiguous residue
 * sequence (per chain), by spatial k-means clusters, or auto (a target lobe size that
 * sets the per-chain count for uniform blobs) - each fit separately and
 * combined as overlapping blobs (no marching cubes) or a watertight blend;
 * see the visual module for details.
 *
 * Three visual flavours, differing only in which atoms each envelope is fit to:
 * - per-unit (default): one envelope per chain, instanced across the assembly natively.
 * - structure: a single envelope over the whole input merged together - a deliberately
 *   coarse overall shape of the entire complex (e.g. a quick globular outline of a capsid
 *   or assembly as one object). Star-convex like the others, so it is a smooth bounding
 *   shape, not a per-chain surface.
 * - assembly: one envelope fit to the asymmetric unit and replicated across the assembly
 *   operators (e.g. a protomer fit once, drawn as the full symmetric assembly).
 */

import { SphericalHarmonicSurfaceMeshVisual, SphericalHarmonicSurfaceMeshParams, StructureSphericalHarmonicSurfaceMeshVisual, AssemblySphericalHarmonicSurfaceMeshVisual, LobesParam } from './mesh';
import { UnitsRepresentation } from '../../mol-repr/structure/units-representation';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ComplexRepresentation, StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../../mol-repr/structure/representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from '../../mol-repr/representation';
import { ThemeRegistryContext } from '../../mol-theme/theme';
import { Structure } from '../../mol-model/structure';
import { BaseGeometry } from '../../mol-geo/geometry/base';

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
    // size the auto target-atoms slider to the biggest chain in the loaded structure (O(units), cheap)
    let maxChainAtoms = 0;
    for (const unit of structure.units) if (unit.elements.length > maxChainAtoms) maxChainAtoms = unit.elements.length;
    return {
        ...SphericalHarmonicSurfaceParams,
        lobes: LobesParam(maxChainAtoms),
    };
}

export type SphericalHarmonicSurfaceRepresentation = StructureRepresentation<SphericalHarmonicSurfaceParams>
export function SphericalHarmonicSurfaceRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, SphericalHarmonicSurfaceParams>): SphericalHarmonicSurfaceRepresentation {
    return Representation.createMulti('Spherical Harmonic Surface', ctx, getParams, StructureRepresentationStateBuilder, SphericalHarmonicSurfaceVisuals as unknown as Representation.Def<Structure, SphericalHarmonicSurfaceParams>);
}

export const SphericalHarmonicSurfaceRepresentationProvider = StructureRepresentationProvider({
    name: 'spherical-harmonic-surface',
    label: 'Spherical Harmonic Surface',
    description: 'Displays a smooth shape envelope approximated by a spherical harmonic expansion of the atom positions (star-convex per lobe; split a chain into per-sequence, k-means, or auto (uniform-size) lobes to follow elongated or threaded shapes, as overlapping blobs or a watertight blend).',
    factory: SphericalHarmonicSurfaceRepresentation,
    getParams: getSphericalHarmonicSurfaceParams,
    defaultValues: PD.getDefaultValues(SphericalHarmonicSurfaceParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'physical' },
    isApplicable: (structure: Structure) => structure.elementCount > 0,
});
