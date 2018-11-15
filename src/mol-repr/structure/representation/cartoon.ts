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
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation, RepresentationParamsGetter } from 'mol-repr/representation';
import { PolymerDirectionVisual, PolymerDirectionParams } from '../visual/polymer-direction-wedge';
import { Structure } from 'mol-model/structure';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { BuiltInSizeThemeName, BuiltInSizeThemeOptions } from 'mol-theme/size';
import { BuiltInColorThemeOptions, BuiltInColorThemeName } from 'mol-theme/color';

const CartoonVisuals = {
    'polymer-trace': (getParams: RepresentationParamsGetter<Structure, PolymerTraceParams>) => UnitsRepresentation('Polymer trace mesh', getParams, PolymerTraceVisual),
    'polymer-gap': (getParams: RepresentationParamsGetter<Structure, PolymerGapParams>) => UnitsRepresentation('Polymer gap cylinder', getParams, PolymerGapVisual),
    'nucleotide-block': (getParams: RepresentationParamsGetter<Structure, NucleotideBlockParams>) => UnitsRepresentation('Nucleotide block mesh', getParams, NucleotideBlockVisual),
    'direction-wedge': (getParams: RepresentationParamsGetter<Structure, PolymerDirectionParams>) => UnitsRepresentation('Polymer direction wedge', getParams, PolymerDirectionVisual)
}
type CartoonVisualName = keyof typeof CartoonVisuals
const CartoonVisualOptions = Object.keys(CartoonVisuals).map(name => [name, name] as [CartoonVisualName, string])

export const CartoonParams = {
    ...PolymerTraceParams,
    ...PolymerGapParams,
    ...NucleotideBlockParams,
    ...PolymerDirectionParams,
    sizeFactor: PD.Numeric(0.2, { min: 0, max: 10, step: 0.01 }),
    sizeTheme: PD.Select<BuiltInSizeThemeName>('uniform', BuiltInSizeThemeOptions),
    colorTheme: PD.Select<BuiltInColorThemeName>('polymer-index', BuiltInColorThemeOptions),
    visuals: PD.MultiSelect<CartoonVisualName>(['polymer-trace', 'polymer-gap', 'nucleotide-block'], CartoonVisualOptions),
}
export type CartoonParams = typeof CartoonParams
export function getCartoonParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(CartoonParams)
}

export type CartoonRepresentation = StructureRepresentation<CartoonParams>
export function CartoonRepresentation(getParams: RepresentationParamsGetter<Structure, CartoonParams>): CartoonRepresentation {
    return Representation.createMulti('Cartoon', getParams, CartoonVisuals as unknown as Representation.Def<Structure, CartoonParams>)
}

export const CartoonRepresentationProvider: StructureRepresentationProvider<CartoonParams> = {
    label: 'Cartoon',
    description: 'Displays a ribbon smoothly following the trace atom of polymers.',
    factory: CartoonRepresentation,
    getParams: getCartoonParams,
    defaultValues: PD.getDefaultValues(CartoonParams)
}