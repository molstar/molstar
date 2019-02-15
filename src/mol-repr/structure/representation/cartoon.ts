/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerTraceVisual,  PolymerTraceParams } from '../visual/polymer-trace-mesh';
import { PolymerGapVisual, PolymerGapParams } from '../visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, NucleotideBlockParams } from '../visual/nucleotide-block-mesh';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation, StructureRepresentationProvider, StructureRepresentationStateBuilder } from '../representation';
import { Representation, RepresentationParamsGetter, RepresentationContext } from 'mol-repr/representation';
import { PolymerDirectionVisual, PolymerDirectionParams } from '../visual/polymer-direction-wedge';
import { Structure, Unit } from 'mol-model/structure';
import { ThemeRegistryContext } from 'mol-theme/theme';

const CartoonVisuals = {
    'polymer-trace': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerTraceParams>) => UnitsRepresentation('Polymer trace mesh', ctx, getParams, PolymerTraceVisual),
    'polymer-gap': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', ctx, getParams, PolymerGapVisual),
    'nucleotide-block': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, NucleotideBlockParams>) => UnitsRepresentation('Nucleotide block mesh', ctx, getParams, NucleotideBlockVisual),
    'direction-wedge': (ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, PolymerDirectionParams>) => UnitsRepresentation('Polymer direction wedge', ctx, getParams, PolymerDirectionVisual)
}
type CartoonVisualName = keyof typeof CartoonVisuals
const CartoonVisualOptions = Object.keys(CartoonVisuals).map(name => [name, name] as [CartoonVisualName, string])

export const CartoonParams = {
    ...PolymerTraceParams,
    ...PolymerGapParams,
    ...NucleotideBlockParams,
    ...PolymerDirectionParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    visuals: PD.MultiSelect<CartoonVisualName>(['polymer-trace', 'polymer-gap', 'nucleotide-block'], CartoonVisualOptions),
}
export type CartoonParams = typeof CartoonParams
export function getCartoonParams(ctx: ThemeRegistryContext, structure: Structure) {
    const params = PD.clone(CartoonParams)
    let hasNucleotides = false
    let hasGaps = false
    structure.units.forEach(u => {
        if (!hasNucleotides && Unit.isAtomic(u) && u.nucleotideElements.length) hasNucleotides = true
        if (!hasGaps && u.gapElements.length) hasGaps = true
    })
    params.visuals.defaultValue = ['polymer-trace']
    if (hasNucleotides) params.visuals.defaultValue.push('nucleotide-block')
    if (hasGaps) params.visuals.defaultValue.push('polymer-gap')
    return params
}

export type CartoonRepresentation = StructureRepresentation<CartoonParams>
export function CartoonRepresentation(ctx: RepresentationContext, getParams: RepresentationParamsGetter<Structure, CartoonParams>): CartoonRepresentation {
    return Representation.createMulti('Cartoon', ctx, getParams, StructureRepresentationStateBuilder, CartoonVisuals as unknown as Representation.Def<Structure, CartoonParams>)
}

export const CartoonRepresentationProvider: StructureRepresentationProvider<CartoonParams> = {
    label: 'Cartoon',
    description: 'Displays a ribbon smoothly following the trace atoms of polymers.',
    factory: CartoonRepresentation,
    getParams: getCartoonParams,
    defaultValues: PD.getDefaultValues(CartoonParams),
    defaultColorTheme: 'polymer-id',
    defaultSizeTheme: 'uniform'
}