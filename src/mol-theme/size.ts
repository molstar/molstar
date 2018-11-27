/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SizeType, LocationSize } from 'mol-geo/geometry/size-data';
import { UniformSizeThemeProvider } from './size/uniform';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ThemeDataContext } from 'mol-theme/theme';
import { PhysicalSizeThemeProvider } from './size/physical';
import { deepEqual } from 'mol-util';

export { SizeTheme }
interface SizeTheme<P extends PD.Params = {}> {
    readonly factory: SizeTheme.Factory<P>
    readonly granularity: SizeType
    readonly size: LocationSize
    readonly props: Readonly<PD.Values<P>>
    readonly description?: string
}
namespace SizeTheme {
    export type Props = { [k: string]: any }
    export type Factory<P extends PD.Params> = (ctx: ThemeDataContext, props: PD.Values<P>) => SizeTheme<P>
    export const EmptyFactory = () => Empty
    export const Empty: SizeTheme<{}> = { factory: EmptyFactory, granularity: 'uniform', size: () => 1, props: {} }

    export function areEqual(themeA: SizeTheme, themeB: SizeTheme) {
        return themeA.factory === themeB.factory && deepEqual(themeA.props, themeB.props)
    }

    export interface Provider<P extends PD.Params> {
        readonly label: string
        readonly factory: Factory<P>
        readonly getParams: (ctx: ThemeDataContext) => P
    }
    export const EmptyProvider: Provider<{}> = { label: '', factory: EmptyFactory, getParams: () => ({}) }

    export class Registry {
        private _list: { name: string, provider: Provider<any> }[] = []
        private _map = new Map<string, Provider<any>>()

        get default() { return this._list[0]; }
        get types(): [string, string][] {
            return this._list.map(e => [e.name, e.provider.label] as [string, string]);
        }

        constructor() {
            Object.keys(BuiltInSizeThemes).forEach(name => {
                const p = (BuiltInSizeThemes as { [k: string]: Provider<any> })[name]
                this.add(name, p)
            })
        }

        add<P extends PD.Params>(name: string, provider: Provider<P>) {
            this._list.push({ name, provider })
            this._map.set(name, provider)
        }

        get<P extends PD.Params>(id: string) {
            return this._map.get(id) || EmptyProvider as unknown as Provider<P>
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

export const BuiltInSizeThemes = {
    'physical': PhysicalSizeThemeProvider,
    'uniform': UniformSizeThemeProvider
}