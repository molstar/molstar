/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from 'mol-util/color';
import { Location } from 'mol-model/location';
import { ColorType } from 'mol-geo/geometry/color-data';
import { CarbohydrateSymbolColorThemeProvider } from './color/carbohydrate-symbol';
import { UniformColorTheme, UniformColorThemeProvider } from './color/uniform';
import { deepEqual } from 'mol-util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ThemeDataContext } from './theme';
import { ChainIdColorThemeProvider } from './color/chain-id';
import { CrossLinkColorThemeProvider } from './color/cross-link';
import { ElementIndexColorThemeProvider } from './color/element-index';
import { ElementSymbolColorThemeProvider } from './color/element-symbol';
import { MoleculeTypeColorThemeProvider } from './color/molecule-type';
import { PolymerIndexColorThemeProvider } from './color/polymer-index';
import { ResidueNameColorThemeProvider } from './color/residue-name';
import { SecondaryStructureColorThemeProvider } from './color/secondary-structure';
import { SequenceIdColorThemeProvider } from './color/sequence-id';
import { ShapeGroupColorThemeProvider } from './color/shape-group';
import { UnitIndexColorThemeProvider } from './color/unit-index';
import { ScaleLegend } from 'mol-util/color/scale';
import { TableLegend } from 'mol-util/color/tables';

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export type ColorThemeProps = { [k: string]: any }

export { ColorTheme }
interface ColorTheme<P extends ColorThemeProps = {}> {
    readonly granularity: ColorType
    readonly color: LocationColor
    readonly props: Readonly<P>
    readonly description?: string
    readonly legend?: Readonly<ScaleLegend | TableLegend>
}
namespace ColorTheme {
    export type Props = { [k: string]: any }
    export const Empty = UniformColorTheme({}, { value: Color(0xCCCCCC) })

    export function areEqual(themeA: ColorTheme, themeB: ColorTheme) {
        return themeA === themeB && deepEqual(themeA.props, themeB.props)
    }

    export interface Provider<P extends PD.Params> {
        readonly label: string
        readonly factory: (ctx: ThemeDataContext, props: PD.Values<P>) => ColorTheme<PD.Values<P>>
        readonly getParams: (ctx: ThemeDataContext) => P
    }
    export const EmptyProvider: Provider<{}> = { label: '', factory: () => Empty, getParams: () => ({}) }

    export class Registry {
        private _list: { name: string, provider: Provider<any> }[] = []
        private _map = new Map<string, Provider<any>>()

        get default() { return this._list[0]; }
        get types(): [string, string][] {
            return this._list.map(e => [e.name, e.provider.label] as [string, string]);
        }

        constructor() {
            Object.keys(BuiltInColorThemes).forEach(name => {
                const p = (BuiltInColorThemes as { [k: string]: Provider<any> })[name]
                this.add(name, p)
            })
        }

        add<P extends PD.Params>(name: string, provider: Provider<P>) {
            this._list.push({ name, provider })
            this._map.set(name, provider)
        }

        get<P extends PD.Params>(name: string): Provider<P> {
            return this._map.get(name) || EmptyProvider as unknown as Provider<P>
        }

        create(id: string, ctx: ThemeDataContext, props = {}) {
            const provider = this.get(id)
            return provider ? provider.factory(ctx, { ...PD.getDefaultValues(provider.getParams(ctx)), ...props }) : Empty
        }

        get list() {
            return this._list
        }
    }
}

export const BuiltInColorThemes = {
    'carbohydrate-symbol': CarbohydrateSymbolColorThemeProvider,
    'chain-id': ChainIdColorThemeProvider,
    'cross-link': CrossLinkColorThemeProvider,
    'element-index': ElementIndexColorThemeProvider,
    'element-symbol': ElementSymbolColorThemeProvider,
    'molecule-type': MoleculeTypeColorThemeProvider,
    'polymer-index': PolymerIndexColorThemeProvider,
    'residue-name': ResidueNameColorThemeProvider,
    'secondary-structure': SecondaryStructureColorThemeProvider,
    'sequence-id': SequenceIdColorThemeProvider,
    'shape-group': ShapeGroupColorThemeProvider,
    'unit-index': UnitIndexColorThemeProvider,
    'uniform': UniformColorThemeProvider,
}
export type BuiltInColorThemeName = keyof typeof BuiltInColorThemes
export const BuiltInColorThemeNames = Object.keys(BuiltInColorThemes)
export const BuiltInColorThemeOptions = BuiltInColorThemeNames.map(n => [n, n] as [BuiltInColorThemeName, string])
export const getBuiltInColorThemeParams = (name: string, ctx: ThemeDataContext = {}) => PD.Group((BuiltInColorThemes as { [k: string]: ColorTheme.Provider<any> })[name].getParams(ctx))