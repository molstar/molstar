/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from './color';
import { SizeTheme } from './size';
import { Structure } from '../mol-model/structure';
import { Volume } from '../mol-model/volume';
import { ParamDefinition as PD } from '../mol-util/param-definition';
import { Shape } from '../mol-model/shape';
import { CustomProperty } from '../mol-model-props/common/custom-property';
import { objectForEach } from '../mol-util/object';

export interface ThemeRegistryContext {
    colorThemeRegistry: ColorTheme.Registry
    sizeThemeRegistry: SizeTheme.Registry
}

export interface ThemeDataContext {
    [k: string]: any
    structure?: Structure
    volume?: Volume
    shape?: Shape
}

export { Theme };

interface Theme {
    color: ColorTheme<any>
    size: SizeTheme<any>
    // label: LabelTheme // TODO
}

namespace Theme {
    type Props = { [k: string]: any }

    export function create(ctx: ThemeRegistryContext, data: ThemeDataContext, props: Props, theme?: Theme) {
        theme = theme || createEmpty();

        const colorProps = props.colorTheme as PD.NamedParams;
        const sizeProps = props.sizeTheme as PD.NamedParams;

        theme.color = ctx.colorThemeRegistry.create(colorProps.name, data, colorProps.params);
        theme.size = ctx.sizeThemeRegistry.create(sizeProps.name, data, sizeProps.params);

        return theme;
    }

    export function createEmpty(): Theme {
        return { color: ColorTheme.Empty, size: SizeTheme.Empty };
    }

    export async function ensureDependencies(ctx: CustomProperty.Context, theme: ThemeRegistryContext, data: ThemeDataContext, props: Props) {
        await theme.colorThemeRegistry.get(props.colorTheme.name).ensureCustomProperties?.attach(ctx, data);
        await theme.sizeThemeRegistry.get(props.sizeTheme.name).ensureCustomProperties?.attach(ctx, data);
    }

    export function releaseDependencies(theme: ThemeRegistryContext, data: ThemeDataContext, props: Props) {
        theme.colorThemeRegistry.get(props.colorTheme.name).ensureCustomProperties?.detach(data);
        theme.sizeThemeRegistry.get(props.sizeTheme.name).ensureCustomProperties?.detach(data);
    }
}

//

export interface ThemeProvider<T extends ColorTheme<P> | SizeTheme<P>, P extends PD.Params, Id extends string = string> {
    readonly name: Id
    readonly label: string
    readonly category: string
    readonly factory: (ctx: ThemeDataContext, props: PD.Values<P>) => T
    readonly getParams: (ctx: ThemeDataContext) => P
    readonly defaultValues: PD.Values<P>
    readonly isApplicable: (ctx: ThemeDataContext) => boolean
    readonly ensureCustomProperties?: {
        attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => Promise<void>,
        detach: (data: ThemeDataContext) => void
    }
}

function getTypes(list: { name: string, provider: ThemeProvider<any, any> }[]) {
    return list.map(e => [e.name, e.provider.label, e.provider.category] as [string, string, string]);
}

export class ThemeRegistry<T extends ColorTheme<any> | SizeTheme<any>> {
    private _list: { name: string, provider: ThemeProvider<T, any> }[] = []
    private _map = new Map<string, ThemeProvider<T, any>>()
    private _name = new Map<ThemeProvider<T, any>, string>()

    get default() { return this._list[0]; }
    get list() { return this._list; }
    get types(): [string, string, string][] { return getTypes(this._list); }

    constructor(builtInThemes: { [k: string]: ThemeProvider<T, any> }, private emptyProvider: ThemeProvider<T, any>) {
        objectForEach(builtInThemes, (p, k) => {
            if (p.name !== k) throw new Error(`Fix build in themes to have matching names. ${p.name} ${k}`);
            this.add(p as any);
        });
    }

    private sort() {
        this._list.sort((a, b) => {
            if (a.provider.category === b.provider.category) {
                return a.provider.label < b.provider.label ? -1 : a.provider.label > b.provider.label ? 1 : 0;
            }
            return a.provider.category < b.provider.category ? -1 : 1;
        });
    }

    add<P extends PD.Params>(provider: ThemeProvider<T, P>) {
        if (this._map.has(provider.name)) {
            throw new Error(`${provider.name} already registered.`);
        }

        const name = provider.name;
        this._list.push({ name, provider });
        this._map.set(name, provider);
        this._name.set(provider, name);
        this.sort();
    }

    remove(provider: ThemeProvider<T, any>) {
        this._list.splice(this._list.findIndex(e => e.name === provider.name), 1);
        const p = this._map.get(provider.name);
        if (p) {
            this._map.delete(provider.name);
            this._name.delete(p);
        }
    }

    has(provider: ThemeProvider<T, any>): boolean {
        return this._map.has(provider.name);
    }

    get<P extends PD.Params>(name: string): ThemeProvider<T, P> {
        return this._map.get(name) || this.emptyProvider;
    }

    getName(provider: ThemeProvider<T, any>): string {
        if (!this._name.has(provider)) throw new Error(`'${provider.label}' is not a registered theme provider.`);
        return this._name.get(provider)!;
    }


    create(name: string, ctx: ThemeDataContext, props = {}) {
        const provider = this.get(name);
        return provider.factory(ctx, { ...PD.getDefaultValues(provider.getParams(ctx)), ...props });
    }

    getApplicableList(ctx: ThemeDataContext) {
        return this._list.filter(e => e.provider.isApplicable(ctx));
    }

    getApplicableTypes(ctx: ThemeDataContext) {
        return getTypes(this.getApplicableList(ctx));
    }
}