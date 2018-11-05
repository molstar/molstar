/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CrossLinkRestraintVisual, CrossLinkRestraintParams } from '../visual/cross-link-restraint-cylinder';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { paramDefaultValues, SelectParam, NumberParam } from 'mol-util/parameter';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const DistanceRestraintParams = {
    ...CrossLinkRestraintParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 0.25, 0, 0.05, 20),
}
export const DefaultDistanceRestraintProps = paramDefaultValues(DistanceRestraintParams)
export type DistanceRestraintProps = typeof DefaultDistanceRestraintProps

export type DistanceRestraintRepresentation = StructureRepresentation<DistanceRestraintProps>

export function DistanceRestraintRepresentation(): DistanceRestraintRepresentation {
    return Representation.createMulti('Distance restraint', DistanceRestraintParams, DefaultDistanceRestraintProps, [
        ComplexRepresentation('Cross-link restraint', CrossLinkRestraintVisual)
    ] as StructureRepresentation<DistanceRestraintProps>[])
}