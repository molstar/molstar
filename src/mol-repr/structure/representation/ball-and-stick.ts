/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { IntraUnitLinkVisual, IntraUnitLinkParams } from '../visual/intra-unit-link-cylinder';
import { InterUnitLinkVisual, InterUnitLinkParams } from '../visual/inter-unit-link-cylinder';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { UnitsRepresentation } from '../units-representation';
import { ComplexRepresentation } from '../complex-representation';
import { StructureRepresentation, StructureRepresentationProvider } from '../representation';
import { Representation, RepresentationParamsGetter } from 'mol-repr/representation';
import { ThemeRegistryContext } from 'mol-theme/theme';
import { Structure } from 'mol-model/structure';
import { BuiltInSizeThemeName, BuiltInSizeThemeOptions } from 'mol-theme/size';
import { BuiltInColorThemeName, BuiltInColorThemeOptions } from 'mol-theme/color';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

type ParamsGetter = RepresentationParamsGetter<Structure>
const BallAndStickVisuals = {
    'element-sphere': (getParams: ParamsGetter) => UnitsRepresentation('Element sphere mesh', getParams, ElementSphereVisual),
    'intra-link': (getParams: ParamsGetter) => UnitsRepresentation('Intra-unit link cylinder', getParams, IntraUnitLinkVisual),
    'inter-link': (getParams: ParamsGetter) => ComplexRepresentation('Inter-unit link cylinder', getParams, InterUnitLinkVisual),
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

export function BallAndStickRepresentation(getParams: ParamsGetter): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', getParams, BallAndStickVisuals as unknown as Representation.Def<Structure, BallAndStickProps>)
}

export const BallAndStickRepresentationProvider: StructureRepresentationProvider<typeof BallAndStickParams> = {
    factory: BallAndStickRepresentation, getParams: getBallAndStickParams
}