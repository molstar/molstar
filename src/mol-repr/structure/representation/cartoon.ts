/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerTraceVisual,  PolymerTraceParams, PolymerTraceProps } from '../visual/polymer-trace-mesh';
import { PolymerGapVisual, PolymerGapParams, PolymerGapProps } from '../visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, NucleotideBlockParams, NucleotideBlockProps } from '../visual/nucleotide-block-mesh';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation } from 'mol-repr/representation';
import { PolymerDirectionVisual, PolymerDirectionParams, PolymerDirectionProps } from '../visual/polymer-direction-wedge';
import { Structure } from 'mol-model/structure';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { BuiltInSizeThemeName, BuiltInSizeThemeOptions } from 'mol-theme/size';
import { BuiltInColorThemeOptions, BuiltInColorThemeName } from 'mol-theme/color';

const CartoonVisuals = {
    'polymer-trace': (defaultProps: PolymerTraceProps) => UnitsRepresentation('Polymer trace mesh', defaultProps, PolymerTraceVisual),
    'polymer-gap': (defaultProps: PolymerGapProps) => UnitsRepresentation('Polymer gap cylinder', defaultProps, PolymerGapVisual),
    'nucleotide-block': (defaultProps: NucleotideBlockProps) => UnitsRepresentation('Nucleotide block mesh', defaultProps, NucleotideBlockVisual),
    'direction-wedge': (defaultProps: PolymerDirectionProps) => UnitsRepresentation('Polymer direction wedge', defaultProps, PolymerDirectionVisual)
}
type CartoonVisualName = keyof typeof CartoonVisuals
const CartoonVisualOptions = Object.keys(CartoonVisuals).map(name => [name, name] as [CartoonVisualName, string])

export const CartoonParams = {
    ...PolymerTraceParams,
    ...PolymerGapParams,
    ...NucleotideBlockParams,
    ...PolymerDirectionParams,
    sizeFactor: PD.Numeric('Size Factor', '', 0.2, 0, 10, 0.01),
    sizeTheme: PD.Select<BuiltInSizeThemeName>('Size Theme', '', 'uniform', BuiltInSizeThemeOptions),
    colorTheme: PD.Select<BuiltInColorThemeName>('Color Theme', '', 'polymer-index', BuiltInColorThemeOptions),
    visuals: PD.MultiSelect<CartoonVisualName>('Visuals', '', ['polymer-trace', 'polymer-gap', 'nucleotide-block'], CartoonVisualOptions),
}
export function getCartoonParams(ctx: ThemeRegistryContext, structure: Structure) {
    return CartoonParams // TODO return copy
}
export type CartoonProps = PD.DefaultValues<typeof CartoonParams>

export type CartoonRepresentation = StructureRepresentation<CartoonProps>
export function CartoonRepresentation(defaultProps: CartoonProps): CartoonRepresentation {
    return Representation.createMulti('Cartoon', defaultProps, CartoonVisuals as unknown as Representation.Def<CartoonProps>)
}

export const CartoonRepresentationProvider: StructureRepresentationProvider<typeof CartoonParams> = {
    factory: CartoonRepresentation, params: getCartoonParams
}