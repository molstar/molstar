/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type Structure, Unit } from '@molstar/model/model/structure';
import { Representation, type RepresentationContext, type RepresentationParamsGetter } from '../../representation.js';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation.js';
import { UnitsRepresentation } from '../units-representation.js';
import { NucleotideBlockParams, NucleotideBlockVisual } from '../visual/nucleotide-block-mesh.js';
import { NucleotideRingParams, NucleotideRingVisual } from '../visual/nucleotide-ring-mesh.js';
import { NucleotideAtomicRingFillParams, NucleotideAtomicRingFillVisual } from '../visual/nucleotide-atomic-ring-fill.js';
import { NucleotideAtomicBondParams, NucleotideAtomicBondVisual } from '../visual/nucleotide-atomic-bond.js';
import { NucleotideAtomicElementParams, NucleotideAtomicElementVisual } from '../visual/nucleotide-atomic-element.js';
import { PolymerDirectionParams, PolymerDirectionVisual } from '../visual/polymer-direction-wedge.js';
import { PolymerGapParams, PolymerGapVisual } from '../visual/polymer-gap-cylinder.js';
import { PolymerTraceParams, PolymerTraceVisual } from '../visual/polymer-trace-mesh.js';
import { SecondaryStructureProvider } from '@molstar/model/props/computed/secondary-structure';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import { HelixOrientationProvider } from '@molstar/model/props/computed/helix-orientation';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';

const CartoonVisuals = {
    'polymer-trace': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerTraceParams>) => UnitsRepresentation('Polymer trace mesh', ctx, getParams, PolymerTraceVisual),
    'polymer-gap': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', ctx, getParams, PolymerGapVisual),
    'nucleotide-block': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideBlockParams>) => UnitsRepresentation('Nucleotide block mesh', ctx, getParams, NucleotideBlockVisual),
    'nucleotide-ring': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideRingParams>) => UnitsRepresentation('Nucleotide ring mesh', ctx, getParams, NucleotideRingVisual),
    'nucleotide-atomic-ring-fill': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideAtomicRingFillParams>) => UnitsRepresentation('Nucleotide atomic ring fill', ctx, getParams, NucleotideAtomicRingFillVisual),
    'nucleotide-atomic-bond': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideAtomicBondParams>) => UnitsRepresentation('Nucleotide atomic bond', ctx, getParams, NucleotideAtomicBondVisual),
    'nucleotide-atomic-element': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideAtomicElementParams>) => UnitsRepresentation('Nucleotide atomic element', ctx, getParams, NucleotideAtomicElementVisual),
    'direction-wedge': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerDirectionParams>) => UnitsRepresentation('Polymer direction wedge', ctx, getParams, PolymerDirectionVisual),
};

export const CartoonParams = {
    ...PolymerTraceParams,
    ...PolymerGapParams,
    ...NucleotideBlockParams,
    ...NucleotideRingParams,
    ...NucleotideAtomicBondParams,
    ...NucleotideAtomicElementParams,
    ...NucleotideAtomicRingFillParams,
    ...PolymerDirectionParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    visuals: PD.MultiSelect(['polymer-trace', 'polymer-gap', 'nucleotide-ring', 'nucleotide-atomic-ring-fill', 'nucleotide-atomic-bond', 'nucleotide-atomic-element'], PD.objectToOptions(CartoonVisuals)),
    bumpFrequency: PD.Numeric(2, { min: 0, max: 10, step: 0.1 }, BaseGeometry.ShadingCategory),
    solidInterior: PD.Boolean(false, BaseGeometry.ShadingCategory),
    density: PD.Numeric(0.1, { min: 0, max: 1, step: 0.01 }, BaseGeometry.ShadingCategory),
    colorMode: PD.Select('default', PD.arrayToOptions(['default', 'interpolate'] as const), { ...BaseGeometry.ShadingCategory, isHidden: true }),
};

export type CartoonParams = typeof CartoonParams
export function getCartoonParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(CartoonParams);
    let hasNucleotides = false;
    let hasGaps = false;
    structure.units.forEach(u => {
        if (!hasNucleotides && Unit.isAtomic(u) && u.nucleotideElements.length) hasNucleotides = true;
        if (!hasGaps && u.gapElements.length) hasGaps = true;
    });
    params.visuals.defaultValue = ['polymer-trace'];
    if (hasNucleotides) params.visuals.defaultValue.push('nucleotide-ring');
    if (hasGaps) params.visuals.defaultValue.push('polymer-gap');
    return params;
}

export type CartoonRepresentation = StructureRepresentation<CartoonParams>
export function CartoonRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CartoonParams>): CartoonRepresentation {
    return Representation.createMulti('Cartoon', ctx, getParams, StructureRepresentationStateBuilder, CartoonVisuals as unknown as Representation.Def<Structure, CartoonParams>);
}

export const CartoonRepresentationProvider = StructureRepresentationProvider({
    name: 'cartoon',
    label: 'Cartoon',
    description: 'Displays ribbons, planks, tubes smoothly following the trace atoms of polymers.',
    factory: CartoonRepresentation,
    getParams: getCartoonParams,
    defaultValues: PD.getDefaultValues(CartoonParams),
    defaultColorTheme: { name: 'chain-id' },
    defaultSizeTheme: { name: 'uniform' },
    isApplicable: (structure: Structure) => structure.polymerResidueCount > 0,
    ensureCustomProperties: {
        attach: async (ctx: CustomProperty.Context, structure: Structure) => {
            await SecondaryStructureProvider.attach(ctx, structure, void 0, true);
            for (const m of structure.models) {
                await HelixOrientationProvider.attach(ctx, m, void 0, true);
            }
        },
        detach: (data) => {
            SecondaryStructureProvider.ref(data, false);
            for (const m of data.models) {
                HelixOrientationProvider.ref(m, false);
            }
        }
    }
});
