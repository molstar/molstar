/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../mol-util/color';
import { Location } from '../mol-model/location';
import { ColorType } from '../mol-geo/geometry/color-data';
import { CarbohydrateSymbolColorThemeProvider } from './color/carbohydrate-symbol';
import { UniformColorThemeProvider } from './color/uniform';
import { deepEqual } from '../mol-util';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { ThemeDataContext, ThemeRegistry, ThemeProvider } from './theme';
import { ChainIdColorThemeProvider } from './color/chain-id';
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
import { ScaleLegend, TableLegend } from '../mol-util/legend';
import { UncertaintyColorThemeProvider } from './color/uncertainty';
import { EntitySourceColorThemeProvider } from './color/entity-source';
import { IllustrativeColorThemeProvider } from './color/illustrative';
import { HydrophobicityColorThemeProvider } from './color/hydrophobicity';
import { ModelIndexColorThemeProvider } from './color/model-index';
import { OccupancyColorThemeProvider } from './color/occupancy';
import { OperatorNameColorThemeProvider } from './color/operator-name';
import { OperatorHklColorThemeProvider } from './color/operator-hkl';

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export { ColorTheme };
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
    export const enum Category {
        Atom = 'Atom Property',
        Chain = 'Chain Property',
        Residue = 'Residue Property',
        Symmetry = 'Symmetry',
        Validation = 'Validation',
        Misc = 'Miscellaneous',
    }

    export type Props = { [k: string]: any }
    export type Factory<P extends PD.Params> = (ctx: ThemeDataContext, props: PD.Values<P>) => ColorTheme<P>
    export const EmptyFactory = () => Empty;
    const EmptyColor = Color(0xCCCCCC);
    export const Empty: ColorTheme<{}> = {
        factory: EmptyFactory,
        granularity: 'uniform',
        color: () => EmptyColor,
        props: {}
    };

    export function areEqual(themeA: ColorTheme<any>, themeB: ColorTheme<any>) {
        return themeA.contextHash === themeB.contextHash && themeA.factory === themeB.factory && deepEqual(themeA.props, themeB.props);
    }

    export interface Provider<P extends PD.Params = any, Id extends string = string> extends ThemeProvider<ColorTheme<P>, P, Id> { }
    export const EmptyProvider: Provider<{}> = { name: '', label: '', category: '', factory: EmptyFactory, getParams: () => ({}), defaultValues: {}, isApplicable: () => true };

    export type Registry = ThemeRegistry<ColorTheme<any>>
    export function createRegistry() {
        return new ThemeRegistry(BuiltIn as { [k: string]: Provider<any> }, EmptyProvider);
    }

    export const BuiltIn = {
        'carbohydrate-symbol': CarbohydrateSymbolColorThemeProvider,
        'chain-id': ChainIdColorThemeProvider,
        'element-index': ElementIndexColorThemeProvider,
        'element-symbol': ElementSymbolColorThemeProvider,
        'entity-source': EntitySourceColorThemeProvider,
        'hydrophobicity': HydrophobicityColorThemeProvider,
        'illustrative': IllustrativeColorThemeProvider,
        'model-index': ModelIndexColorThemeProvider,
        'molecule-type': MoleculeTypeColorThemeProvider,
        'occupancy': OccupancyColorThemeProvider,
        'operator-hkl': OperatorHklColorThemeProvider,
        'operator-name': OperatorNameColorThemeProvider,
        'polymer-id': PolymerIdColorThemeProvider,
        'polymer-index': PolymerIndexColorThemeProvider,
        'residue-name': ResidueNameColorThemeProvider,
        'secondary-structure': SecondaryStructureColorThemeProvider,
        'sequence-id': SequenceIdColorThemeProvider,
        'shape-group': ShapeGroupColorThemeProvider,
        'uncertainty': UncertaintyColorThemeProvider,
        'unit-index': UnitIndexColorThemeProvider,
        'uniform': UniformColorThemeProvider,
    };
    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type ParamValues<C extends ColorTheme.Provider<any>> = C extends ColorTheme.Provider<infer P> ? PD.Values<P> : never
    export type BuiltInParams<T extends BuiltIn> = Partial<ParamValues<_BuiltIn[T]>>
}

export function ColorThemeProvider<P extends PD.Params, Id extends string>(p: ColorTheme.Provider<P, Id>): ColorTheme.Provider<P, Id> { return p; }