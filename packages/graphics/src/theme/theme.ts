/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ColorTheme } from './color.js';
import { SizeTheme } from './size.js';
import type { Structure } from '@molstar/model/model/structure';
import type { Volume } from '@molstar/model/model/volume';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Shape } from '@molstar/model/model/shape';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { ColorType } from '@molstar/graphics/geo/geometry/color-data';
import type { Location } from '@molstar/model/model/location';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';

export interface ThemeRegistryContext {
  colorThemeRegistry: ColorTheme.Registry;
  sizeThemeRegistry: SizeTheme.Registry;
}

export interface ThemeDataContext {
  [k: string]: any;
  structure?: Structure;
  volume?: Volume;
  shape?: Shape;
  particles?: ParticleList;
  /** Hint to request support for specific kinds of locations */
  locationKinds?: ReadonlyArray<Location['kind']>;
}

export { Theme };

interface Theme {
  color: ColorTheme<any, any>;
  size: SizeTheme<any>;
  // label: LabelTheme // TODO
}

namespace Theme {
  type Props = { [k: string]: any };

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

  export async function ensureDependencies(
    ctx: CustomProperty.Context,
    theme: ThemeRegistryContext,
    data: ThemeDataContext,
    props: Props,
  ) {
    await theme.colorThemeRegistry.get(props.colorTheme.name).ensureCustomProperties?.attach(ctx, data);
    await theme.sizeThemeRegistry.get(props.sizeTheme.name).ensureCustomProperties?.attach(ctx, data);
  }

  export function releaseDependencies(theme: ThemeRegistryContext, data: ThemeDataContext, props: Props) {
    theme.colorThemeRegistry.get(props.colorTheme.name).ensureCustomProperties?.detach(data);
    theme.sizeThemeRegistry.get(props.sizeTheme.name).ensureCustomProperties?.detach(data);
  }
}

//

export interface ThemeProvider<
  T extends ColorTheme<P, G> | SizeTheme<P>,
  P extends PD.Params,
  Id extends string = string,
  G extends ColorType = ColorType,
> {
  readonly name: Id;
  readonly label: string;
  readonly category: string;
  readonly factory: (ctx: ThemeDataContext, props: PD.Values<P>) => T;
  readonly getParams: (ctx: ThemeDataContext) => P;
  readonly defaultValues: PD.Values<P>;
  readonly isApplicable: (ctx: ThemeDataContext) => boolean;
  readonly ensureCustomProperties?: {
    attach: (ctx: CustomProperty.Context, data: ThemeDataContext) => Promise<void>;
    detach: (data: ThemeDataContext) => void;
  };
}

function getTypes(list: { name: string; provider: ThemeProvider<any, any> }[]) {
  return list.map((e) => [e.name, e.provider.label, e.provider.category] as [string, string, string]);
}

export class ThemeRegistry<T extends ColorTheme<any, any> | SizeTheme<any>> {
  private _list: { name: string; provider: ThemeProvider<T, any> }[] = [];
  private _map = new Map<string, ThemeProvider<T, any>>();
  private _name = new Map<ThemeProvider<T, any>, string>();
  private _count = new Map<ThemeProvider<T, any>, number>();

  /** The first registered entry in sorted order, or `undefined` when the registry is empty. */
  get default(): { name: string; provider: ThemeProvider<T, any> } | undefined {
    return this._list[0];
  }
  get list() {
    return this._list;
  }
  get types(): [string, string, string][] {
    return getTypes(this._list);
  }

  constructor(private emptyProvider: ThemeProvider<T, any>) {}

  private sort() {
    this._list.sort((a, b) => {
      if (a.provider.category === b.provider.category) {
        return a.provider.label < b.provider.label ? -1 : a.provider.label > b.provider.label ? 1 : 0;
      }
      return a.provider.category < b.provider.category ? -1 : 1;
    });
  }

  /** Returns the message `add(provider)` would throw, or `undefined` when it would succeed. Changes nothing. */
  findConflict(provider: ThemeProvider<T, any>): string | undefined {
    const existing = this._map.get(provider.name);
    if (existing && existing !== provider) {
      return `Theme '${provider.name}' is already registered with a different provider.`;
    }
    return undefined;
  }

  /** Registering the same object again increments a count; a different object under the same name throws. */
  add<P extends PD.Params>(provider: ThemeProvider<T, P>) {
    const conflict = this.findConflict(provider);
    if (conflict) throw new Error(conflict);

    const count = this._count.get(provider);
    if (count !== undefined) {
      this._count.set(provider, count + 1);
      return;
    }

    const name = provider.name;
    this._list.push({ name, provider });
    this._map.set(name, provider);
    this._name.set(provider, name);
    this._count.set(provider, 1);
    this.sort();
  }

  /** Decrements the count and removes the provider at zero. Removing an unknown provider is a no-op. */
  remove(provider: ThemeProvider<T, any>) {
    const count = this._count.get(provider);
    if (count === undefined) return;
    if (count > 1) {
      this._count.set(provider, count - 1);
      return;
    }

    this._count.delete(provider);
    this._map.delete(provider.name);
    this._name.delete(provider);
    const i = this._list.findIndex((e) => e.provider === provider);
    if (i >= 0) this._list.splice(i, 1);
  }

  has(nameOrProvider: string | ThemeProvider<T, any>): boolean {
    return typeof nameOrProvider === 'string' ? this._map.has(nameOrProvider) : this._count.has(nameOrProvider);
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
    return this._list.filter((e) => e.provider.isApplicable(ctx));
  }

  getApplicableTypes(ctx: ThemeDataContext) {
    return getTypes(this.getApplicableList(ctx));
  }

  clear() {
    this._list.length = 0;
    this._map.clear();
    this._name.clear();
    this._count.clear();
  }
}
