/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '../mol-util/color';
import { Location } from '../mol-model/location';
import { ColorType, ColorTypeDirect, ColorTypeGrid, ColorTypeLocation } from '../mol-geo/geometry/color-data';
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
import { ResidueChargeColorThemeProvider } from './color/residue-charge';
import { SecondaryStructureColorThemeProvider } from './color/secondary-structure';
import { SequenceIdColorThemeProvider } from './color/sequence-id';
import { ShapeGroupColorThemeProvider } from './color/shape-group';
import { UnitIndexColorThemeProvider } from './color/unit-index';
import { ScaleLegend, TableLegend } from '../mol-util/legend';
import { UncertaintyColorThemeProvider } from './color/uncertainty';
import { EntitySourceColorThemeProvider } from './color/entity-source';
import { IllustrativeColorThemeProvider } from './color/illustrative';
import { HydrophobicityColorThemeProvider } from './color/hydrophobicity';
import { TrajectoryIndexColorThemeProvider } from './color/trajectory-index';
import { OccupancyColorThemeProvider } from './color/occupancy';
import { OperatorNameColorThemeProvider } from './color/operator-name';
import { OperatorHklColorThemeProvider } from './color/operator-hkl';
import { PartialChargeColorThemeProvider } from './color/partial-charge';
import { AtomIdColorThemeProvider } from './color/atom-id';
import { EntityIdColorThemeProvider } from './color/entity-id';
import type { Texture, TextureFilter } from '../mol-gl/webgl/texture';
import { VolumeValueColorThemeProvider } from './color/volume-value';
import { Vec3, Vec4 } from '../mol-math/linear-algebra';
import { ModelIndexColorThemeProvider } from './color/model-index';
import { StructureIndexColorThemeProvider } from './color/structure-index';
import { VolumeSegmentColorThemeProvider } from './color/volume-segment';
import { ExternalVolumeColorThemeProvider } from './color/external-volume';
import { ColorThemeCategory } from './color/categories';
import { CartoonColorThemeProvider } from './color/cartoon';
import { FormalChargeColorThemeProvider } from './color/formal-charge';
import { ExternalStructureColorThemeProvider } from './color/external-structure';
import { ColorListEntry } from '../mol-util/color/color';
import { getPrecision } from '../mol-util/number';
import { SortedArray } from '../mol-data/int/sorted-array';
import { normalize } from '../mol-math/interpolate';
import { VolumeInstanceColorThemeProvider } from './color/volume-instance';

export type LocationColor = (location: Location, isSecondary: boolean) => Color

export interface ColorVolume {
    colors: Texture
    dimension: Vec3
    transform: Vec4
}

export { ColorTheme };

type ColorThemeShared<P extends PD.Params, G extends ColorType> = {
    readonly factory: ColorTheme.Factory<P, G>
    readonly props: Readonly<PD.Values<P>>
    /**
     * if palette is defined, 24bit RGB color value normalized to interval [0, 1]
     * is used as index to the colors
     */
    readonly palette?: Readonly<ColorTheme.Palette>
    readonly preferSmoothing?: boolean
    readonly contextHash?: number
    readonly description?: string
    readonly legend?: Readonly<ScaleLegend | TableLegend>
}

type ColorThemeLocation<P extends PD.Params> = {
    readonly granularity: ColorTypeLocation
    readonly color: LocationColor
} & ColorThemeShared<P, ColorTypeLocation>

type ColorThemeGrid<P extends PD.Params> = {
    readonly granularity: ColorTypeGrid
    readonly grid: ColorVolume
} & ColorThemeShared<P, ColorTypeGrid>

type ColorThemeDirect<P extends PD.Params> = {
    readonly granularity: ColorTypeDirect
} & ColorThemeShared<P, ColorTypeDirect>

type ColorTheme<P extends PD.Params, G extends ColorType = ColorTypeLocation> =
    G extends ColorTypeLocation ? ColorThemeLocation<P> :
        G extends ColorTypeGrid ? ColorThemeGrid<P> :
            G extends ColorTypeDirect ? ColorThemeDirect<P> : never

namespace ColorTheme {
    export const Category = ColorThemeCategory;

    export interface Palette {
        colors: Color[],
        filter?: TextureFilter,
        domain?: [number, number],
        defaultColor?: Color,
    }

    export function Palette(list: ColorListEntry[], kind: 'set' | 'interpolate', domain?: [number, number], defaultColor?: Color): Palette {
        const colors: Color[] = [];

        const hasOffsets = list.every(c => Array.isArray(c));
        if (hasOffsets) {
            let maxPrecision = 0;
            for (const e of list) {
                if (Array.isArray(e)) {
                    const p = getPrecision(e[1]);
                    if (p > maxPrecision) maxPrecision = p;
                }
            }
            const count = Math.pow(10, maxPrecision);

            const sorted = [...list] as [Color, number][];
            sorted.sort((a, b) => a[1] - b[1]);

            const src = sorted.map(c => c[0]);
            const values = SortedArray.ofSortedArray(sorted.map(c => c[1]));

            const _off: number[] = [];
            for (let i = 0, il = values.length - 1; i < il; ++i) {
                _off.push(values[i] + (values[i + 1] - values[i]) / 2);
            }
            _off.push(values[values.length - 1]);
            const off = SortedArray.ofSortedArray(_off);

            for (let i = 0, il = Math.max(count, list.length); i < il; ++i) {
                const t = normalize(i, 0, count - 1);
                const j = SortedArray.findPredecessorIndex(off, t);
                colors[i] = src[j];
            }
        } else {
            for (const e of list) {
                if (Array.isArray(e)) colors.push(e[0]);
                else colors.push(e);
            }
        }

        return {
            colors,
            filter: kind === 'set' ? 'nearest' : 'linear',
            domain,
            defaultColor,
        };
    }

    export const PaletteScale = (1 << 24) - 2; // reserve (1 << 24) - 1 for undefiend values

    export type Props = { [k: string]: any }
    export type Factory<P extends PD.Params, G extends ColorType> = (ctx: ThemeDataContext, props: PD.Values<P>) => ColorTheme<P, G>
    export const EmptyFactory = () => Empty;
    const EmptyColor = Color(0xCCCCCC);
    export const Empty: ColorTheme<{}> = {
        factory: EmptyFactory,
        granularity: 'uniform',
        color: () => EmptyColor,
        props: {}
    };

    export function areEqual(themeA: ColorTheme<any, any>, themeB: ColorTheme<any, any>) {
        return themeA.contextHash === themeB.contextHash && themeA.factory === themeB.factory && deepEqual(themeA.props, themeB.props);
    }

    export interface Provider<P extends PD.Params = any, Id extends string = string, G extends ColorType = ColorType> extends ThemeProvider<ColorTheme<P, G>, P, Id, G> { }
    export const EmptyProvider: Provider<{}> = { name: '', label: '', category: '', factory: EmptyFactory, getParams: () => ({}), defaultValues: {}, isApplicable: () => true };

    export type Registry = ThemeRegistry<ColorTheme<any, any>>
    export function createRegistry() {
        return new ThemeRegistry(BuiltIn as { [k: string]: Provider<any, any, any> }, EmptyProvider);
    }

    export const BuiltIn = {
        'atom-id': AtomIdColorThemeProvider,
        'carbohydrate-symbol': CarbohydrateSymbolColorThemeProvider,
        'cartoon': CartoonColorThemeProvider,
        'chain-id': ChainIdColorThemeProvider,
        'element-index': ElementIndexColorThemeProvider,
        'element-symbol': ElementSymbolColorThemeProvider,
        'entity-id': EntityIdColorThemeProvider,
        'entity-source': EntitySourceColorThemeProvider,
        'external-structure': ExternalStructureColorThemeProvider,
        'external-volume': ExternalVolumeColorThemeProvider,
        'formal-charge': FormalChargeColorThemeProvider,
        'hydrophobicity': HydrophobicityColorThemeProvider,
        'illustrative': IllustrativeColorThemeProvider,
        'model-index': ModelIndexColorThemeProvider,
        'molecule-type': MoleculeTypeColorThemeProvider,
        'occupancy': OccupancyColorThemeProvider,
        'operator-hkl': OperatorHklColorThemeProvider,
        'operator-name': OperatorNameColorThemeProvider,
        'partial-charge': PartialChargeColorThemeProvider,
        'polymer-id': PolymerIdColorThemeProvider,
        'polymer-index': PolymerIndexColorThemeProvider,
        'residue-charge': ResidueChargeColorThemeProvider,
        'residue-name': ResidueNameColorThemeProvider,
        'secondary-structure': SecondaryStructureColorThemeProvider,
        'sequence-id': SequenceIdColorThemeProvider,
        'shape-group': ShapeGroupColorThemeProvider,
        'structure-index': StructureIndexColorThemeProvider,
        'trajectory-index': TrajectoryIndexColorThemeProvider,
        'uncertainty': UncertaintyColorThemeProvider,
        'unit-index': UnitIndexColorThemeProvider,
        'uniform': UniformColorThemeProvider,
        'volume-instance': VolumeInstanceColorThemeProvider,
        'volume-segment': VolumeSegmentColorThemeProvider,
        'volume-value': VolumeValueColorThemeProvider,
    };
    type _BuiltIn = typeof BuiltIn
    export type BuiltIn = keyof _BuiltIn
    export type ParamValues<C extends ColorTheme.Provider<any>> = C extends ColorTheme.Provider<infer P> ? PD.Values<P> : never
    export type BuiltInParams<T extends BuiltIn> = Partial<ParamValues<_BuiltIn[T]>>
}

export function ColorThemeProvider<P extends PD.Params, Id extends string>(p: ColorTheme.Provider<P, Id>): ColorTheme.Provider<P, Id> { return p; }