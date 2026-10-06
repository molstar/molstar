/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { SizeType, LocationSize } from '@molstar/graphics/geo/geometry/size-data';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type ThemeDataContext, ThemeRegistry, type ThemeProvider } from './theme.js';
import { deepEqual } from '@molstar/core/util';
import type { BuiltInSizeThemes } from './size/catalog.js';

export { SizeTheme };
interface SizeTheme<P extends PD.Params> {
  readonly factory: SizeTheme.Factory<P>;
  readonly granularity: SizeType;
  readonly size: LocationSize;
  readonly props: Readonly<PD.Values<P>>;
  readonly contextHash?: number;
  readonly description?: string;
}
namespace SizeTheme {
  export type Props = { [k: string]: any };
  export type Factory<P extends PD.Params> = (ctx: ThemeDataContext, props: PD.Values<P>) => SizeTheme<P>;
  export const EmptyFactory = () => Empty;
  export const Empty: SizeTheme<{}> = { factory: EmptyFactory, granularity: 'uniform', size: () => 1, props: {} };

  export function areEqual(themeA: SizeTheme<any>, themeB: SizeTheme<any>) {
    return (
      themeA.contextHash === themeB.contextHash &&
      themeA.factory === themeB.factory &&
      deepEqual(themeA.props, themeB.props)
    );
  }

  export interface Provider<P extends PD.Params = any, Id extends string = string>
    extends ThemeProvider<SizeTheme<P>, P, Id> {}
  export const EmptyProvider: Provider<{}> = {
    name: '',
    label: '',
    category: '',
    factory: EmptyFactory,
    getParams: () => ({}),
    defaultValues: {},
    isApplicable: () => true,
  };

  export type Registry = ThemeRegistry<SizeTheme<any>>;
  export function createRegistry() {
    const registry: Registry = new ThemeRegistry(EmptyProvider);
    return registry;
  }

  type _BuiltIn = typeof BuiltInSizeThemes;
  export type BuiltIn = keyof _BuiltIn;
  export type ParamValues<C extends SizeTheme.Provider<any>> =
    C extends SizeTheme.Provider<infer P> ? PD.Values<P> : never;
  export type BuiltInParams<T extends BuiltIn> = Partial<ParamValues<_BuiltIn[T]>>;
}
