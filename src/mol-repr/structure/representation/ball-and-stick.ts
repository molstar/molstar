/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { IntraUnitLinkVisual, IntraUnitLinkParams } from '../visual/intra-unit-link-cylinder';
import { InterUnitLinkVisual, InterUnitLinkParams } from '../visual/inter-unit-link-cylinder';
import { SizeThemeName, SizeThemeOptions } from 'mol-theme/size';
import { paramDefaultValues, SelectParam, NumberParam, MultiSelectParam } from 'mol-util/parameter';
import { UnitKind, UnitKindOptions } from '../visual/util/common';
import { UnitsRepresentation } from '../units-representation';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation } from '../index';
import { Representation } from 'mol-repr';

export const BallAndStickParams = {
    ...ElementSphereParams,
    ...IntraUnitLinkParams,
    ...InterUnitLinkParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 0.2, 0, 10, 0.1),
    sizeFactor: NumberParam('Size Factor', '', 1, 0, 10, 0.1),
    unitKinds: MultiSelectParam<UnitKind>('Unit Kind', '', ['atomic'], UnitKindOptions),
}
export const DefaultBallAndStickProps = paramDefaultValues(BallAndStickParams)
export type BallAndStickProps = typeof DefaultBallAndStickProps

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickProps>

export function BallAndStickRepresentation(): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', BallAndStickParams, DefaultBallAndStickProps, [
        UnitsRepresentation('Element sphere mesh', ElementSphereVisual),
        UnitsRepresentation('Intra-unit link cylinder', IntraUnitLinkVisual),
        ComplexRepresentation('Inter-unit link cylinder', InterUnitLinkVisual)
    ] as unknown as StructureRepresentation<BallAndStickProps>[]) // TODO avoid cast to unknown
}