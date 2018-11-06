/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PolymerTraceVisual,  PolymerTraceParams } from '../visual/polymer-trace-mesh';
import { PolymerGapVisual, PolymerGapParams } from '../visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, NucleotideBlockParams } from '../visual/nucleotide-block-mesh';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';
// import { PolymerDirectionVisual, DefaultPolymerDirectionProps } from '../visual/polymer-direction-wedge';

export const CartoonParams = {
    ...PolymerTraceParams,
    ...PolymerGapParams,
    ...NucleotideBlockParams,
    // ...PolymerDirectionParams,
    sizeTheme: PD.SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: PD.NumberParam('Size Value', '', 0.6, 0, 10, 0.1),
}
export const DefaultCartoonProps = PD.paramDefaultValues(CartoonParams)
export type CartoonProps = typeof DefaultCartoonProps

export type CartoonRepresentation = StructureRepresentation<CartoonProps>

export function CartoonRepresentation(): CartoonRepresentation {
    return Representation.createMulti('Cartoon', CartoonParams, DefaultCartoonProps, [
        UnitsRepresentation('Polymer trace mesh', PolymerTraceVisual),
        UnitsRepresentation('Polymer gap cylinder', PolymerGapVisual),
        UnitsRepresentation('Nucleotide block mesh', NucleotideBlockVisual),
        // UnitsRepresentation('Polymer direction wedge', PolymerDirectionVisual)
    ] as StructureRepresentation<CartoonProps>[])
}