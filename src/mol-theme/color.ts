/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorType } from 'mol-geo/geometry/color-data';
import { CarbohydrateSymbolColorThemeProvider } from './color/carbohydrate-symbol';
import { UniformColorThemeProvider } from './color/uniform';
import { deepEqual } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ThemeDataContext, ThemeRegistry, ThemeProvider } from './theme';
import { ChainIdColorThemeProvider } from './color/chain-id';
import { CrossLinkColorThemeProvider } from './color/cross-link';
import { ElementIndexColorThemeProvider } from './color/element-index';
import { ElementSymbolColorThemeProvider } from './color/element-symbol';
import { MoleculeTypeColorThemeProvider } from './color/molecule-type';
import { PolymerIdColorThemeProvider } from './color/polymer-id';
import { PolymerIndexColorThemeProvider } from './color/polymer-index';
import { ResidueNameColorThemeProvider } from './color/residue-name';
import { SecondaryStructureColorThemeProvider } from './color/secondary-structure';
import { SequenceIdColorThemeProvider } from './color/sequence-id';
import { ShapeGroupColorThemeProvider } from './color/shape-group';
import { UnitIndexColorThemeProvider } from './color/unit-index';
import { ScaleLegend } from 'mol-util/color/scale';
import { TableLegend } from 'mol-util/color/tables';
import { UncertaintyColorThemeProvider } from './color/uncertainty';
import { EntitySourceColorThemeProvider } from './color/entity-source';
import { IllustrativeColorThemeProvider } from './color/illustrative';
import { HydrophobicityColorThemeProvider } from './color/hydrophobicity';

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export { ColorTheme }
interface ColorTheme<P extends PD.Params> {
    readonly factory: ColorTheme.Factory<P>
    readonly granularity: ColorType
    readonly color: LocationColor
    readonly props: Readonly<PD.Values<P>>
    readonly contextHash?: number
    readonly description?: string
    readonly legend?: Readonly<ScaleLegend | TableLegend>
}
namespace ColorTheme {
    export type Props = { [k: string]: any }
    export type Factory<P extends PD.Params> = (ctx: ThemeDataContext, props: PD.Values<P>) => ColorTheme<P>
    export const EmptyFactory = () => Empty
    const EmptyColor = Color(0xCCCCCC)
    export const Empty: ColorTheme<{}> = {
        factory: EmptyFactory,
        granularity: 'uniform',
        color: () => EmptyColor,
        props: {}
    }

    export function areEqual(themeA: ColorTheme<any>, themeB: ColorTheme<any>) {
        return themeA.contextHash === themeB.contextHash && themeA.factory === themeB.factory && deepEqual(themeA.props, themeB.props)
    }

    export interface Provider<P extends PD.Params> extends ThemeProvider<ColorTheme<P>, P> { }
    export const EmptyProvider: Provider<{}> = { label: '', factory: EmptyFactory, getParams: () => ({}), defaultValues: {}, isApplicable: () => true }

    export type Registry = ThemeRegistry<ColorTheme<any>>
    export function createRegistry() {
        return new ThemeRegistry(BuiltInColorThemes as { [k: string]: Provider<any> }, EmptyProvider)
    }

    export type ParamValues<C extends ColorTheme.Provider<any>> = C extends ColorTheme.Provider<infer P> ? PD.Values<P> : never
}

export const BuiltInColorThemes = {
    'carbohydrate-symbol': CarbohydrateSymbolColorThemeProvider,
    'chain-id': ChainIdColorThemeProvider,
    'cross-link': CrossLinkColorThemeProvider,
    'element-index': ElementIndexColorThemeProvider,
    'element-symbol': ElementSymbolColorThemeProvider,
    'entity-source': EntitySourceColorThemeProvider,
    'hydrophobicity': HydrophobicityColorThemeProvider,
    'illustrative': IllustrativeColorThemeProvider,
    'molecule-type': MoleculeTypeColorThemeProvider,
    'polymer-id': PolymerIdColorThemeProvider,
    'polymer-index': PolymerIndexColorThemeProvider,
    'residue-name': ResidueNameColorThemeProvider,
    'secondary-structure': SecondaryStructureColorThemeProvider,
    'sequence-id': SequenceIdColorThemeProvider,
    'shape-group': ShapeGroupColorThemeProvider,
    'uncertainty': UncertaintyColorThemeProvider,
    'unit-index': UnitIndexColorThemeProvider,
    'uniform': UniformColorThemeProvider,
}
export type BuiltInColorThemeName = keyof typeof BuiltInColorThemes