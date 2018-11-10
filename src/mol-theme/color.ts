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
import { ThemeDataContext } from 'mol-theme/theme';
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

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export interface ScaleLegend {
    kind: 'scale-legend'
    minLabel: string,
    maxLabel: string,
    colors: Color[]
}
export function ScaleLegend(minLabel: string, maxLabel: string, colors: Color[]): ScaleLegend {
    return { kind: 'scale-legend', minLabel, maxLabel, colors }
}

export interface TableLegend {
    kind: 'table-legend'
    table: [ string, Color ][]
}
export function TableLegend(table: [ string, Color ][]): TableLegend {
    return { kind: 'table-legend', table }
}

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
        readonly factory: (ctx: ThemeDataContext, props: PD.DefaultValues<P>) => ColorTheme<PD.DefaultValues<P>>
        readonly params: (ctx: ThemeDataContext) => P
    }

    export class Registry {
        private _list: { name: string, provider: Provider<any> }[] = []
        private _map = new Map<string, Provider<any>>()

        constructor() {
            Object.keys(BuiltInColorThemes).forEach(name => {
                const p = (BuiltInColorThemes as { [k: string]: Provider<any> })[name]
                this.add(name, p.factory, p.params)
            })
        }

        add<P extends PD.Params>(name: string, factory: Provider<P>['factory'], params: Provider<P>['params']) {
            const provider = { factory, params } as Provider<P>
            this._list.push({ name, provider })
            this._map.set(name, provider)
        }

        get(id: string) {
            return this._map.get(id)
        }

        create(id: string, ctx: ThemeDataContext, props = {}) {
            const provider = this.get(id)
            return provider ? provider.factory(ctx, { ...PD.getDefaultValues(provider.params(ctx)), ...props }) : Empty
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