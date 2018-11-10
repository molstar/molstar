/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams, ElementSphereProps } from '../visual/element-sphere';
import { IntraUnitLinkVisual, IntraUnitLinkParams } from '../visual/intra-unit-link-cylinder';
import { InterUnitLinkVisual, InterUnitLinkParams, InterUnitLinkProps } from '../visual/inter-unit-link-cylinder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { IntraUnitLinkProps } from '../visual/polymer-gap-cylinder';
import { BuiltInSizeThemeName, BuiltInSizeThemeOptions } from 'mol-theme/size';
import { BuiltInColorThemeName, BuiltInColorThemeOptions } from 'mol-theme/color';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

const BallAndStickVisuals = {
    'element-sphere': (defaultProps: ElementSphereProps) => UnitsRepresentation('Element sphere mesh', defaultProps, ElementSphereVisual),
    'intra-link': (defaultProps: IntraUnitLinkProps) => UnitsRepresentation('Intra-unit link cylinder', defaultProps, IntraUnitLinkVisual),
    'inter-link': (defaultProps: InterUnitLinkProps) => ComplexRepresentation('Inter-unit link cylinder', defaultProps, InterUnitLinkVisual),
}
type BallAndStickVisualName = keyof typeof BallAndStickVisuals
const BallAndStickVisualOptions = Object.keys(BallAndStickVisuals).map(name => [name, name] as [BallAndStickVisualName, string])

export const BallAndStickParams = {
    ...ElementSphereParams,
    ...IntraUnitLinkParams,
    ...InterUnitLinkParams,
    unitKinds: PD.MultiSelect<UnitKind>('Unit Kind', '', ['atomic'], UnitKindOptions),
    sizeFactor: PD.Numeric('Size Factor', '', 0.2, 0.01, 10, 0.01),
    sizeTheme: PD.Select<BuiltInSizeThemeName>('Size Theme', '', 'uniform', BuiltInSizeThemeOptions),
    colorTheme: PD.Select<BuiltInColorThemeName>('Color Theme', '', 'polymer-index', BuiltInColorThemeOptions),
    visuals: PD.MultiSelect<BallAndStickVisualName>('Visuals', '', ['element-sphere', 'intra-link', 'inter-link'], BallAndStickVisualOptions),
}
export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
    return BallAndStickParams // TODO return copy
}
export type BallAndStickProps = PD.DefaultValues<typeof BallAndStickParams>

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickProps>

export function BallAndStickRepresentation(defaultProps: BallAndStickProps): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', defaultProps, BallAndStickVisuals as unknown as Representation.Def<BallAndStickProps>)
}

export const BallAndStickRepresentationProvider: StructureRepresentationProvider<typeof BallAndStickParams> = {
    factory: BallAndStickRepresentation, params: getBallAndStickParams
}