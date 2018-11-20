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
import { BuiltInColorThemeOptions, BuiltInColorThemes, ColorTheme } from 'mol-theme/color';
import { UnitKind, UnitKindOptions } from '../visual/util/common';

const BallAndStickVisuals = {
    'element-sphere': (getParams: RepresentationParamsGetter<Structure, ElementSphereParams>) => UnitsRepresentation('Element sphere mesh', getParams, ElementSphereVisual),
    'intra-link': (getParams: RepresentationParamsGetter<Structure, IntraUnitLinkParams>) => UnitsRepresentation('Intra-unit link cylinder', getParams, IntraUnitLinkVisual),
    'inter-link': (getParams: RepresentationParamsGetter<Structure, InterUnitLinkParams>) => ComplexRepresentation('Inter-unit link cylinder', getParams, InterUnitLinkVisual),
}
type BallAndStickVisualName = keyof typeof BallAndStickVisuals
const BallAndStickVisualOptions = Object.keys(BallAndStickVisuals).map(name => [name, name] as [BallAndStickVisualName, string])

export const BallAndStickParams = {
    ...ElementSphereParams,
    ...IntraUnitLinkParams,
    ...InterUnitLinkParams,
    unitKinds: PD.MultiSelect<UnitKind>(['atomic'], UnitKindOptions),
    sizeFactor: PD.Numeric(0.3, { min: 0.01, max: 10, step: 0.01 }),
    sizeAspectRatio: PD.Numeric(2/3, { min: 0.01, max: 3, step: 0.01 }),
    colorTheme: PD.Mapped('element-symbol', BuiltInColorThemeOptions, name => PD.Group((BuiltInColorThemes as { [k: string]: ColorTheme.Provider<any> })[name].getParams({}))),
    visuals: PD.MultiSelect<BallAndStickVisualName>(['element-sphere', 'intra-link', 'inter-link'], BallAndStickVisualOptions),
}
export type BallAndStickParams = typeof BallAndStickParams
export function getBallAndStickParams(ctx: ThemeRegistryContext, structure: Structure) {
    return PD.clone(BallAndStickParams)
}

export type BallAndStickRepresentation = StructureRepresentation<BallAndStickParams>
export function BallAndStickRepresentation(getParams: RepresentationParamsGetter<Structure, BallAndStickParams>): BallAndStickRepresentation {
    return Representation.createMulti('Ball & Stick', getParams, BallAndStickVisuals as unknown as Representation.Def<Structure, BallAndStickParams>)
}

export const BallAndStickRepresentationProvider: StructureRepresentationProvider<typeof BallAndStickParams> = {
    label: 'Ball & Stick',
    description: 'Displays atoms as spheres and bonds as cylinders.',
    factory: BallAndStickRepresentation,
    getParams: getBallAndStickParams,
    defaultValues: PD.getDefaultValues(BallAndStickParams)
}