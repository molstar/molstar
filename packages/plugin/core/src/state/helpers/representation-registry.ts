/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ThemeRegistryContext } from '@molstar/graphics/theme/theme';

/** Scopes that have their own representation and theme registries. */
export type RepresentationScope = 'structure' | 'volume' | 'particles';

/** The registry surface needed to check names and emptiness (representation or theme registry). */
interface NameRegistry {
  readonly default: unknown;
  has(name: string): boolean;
}

interface ScopeRegistries {
  registry: NameRegistry;
  themes: ThemeRegistryContext;
}

/**
 * Params of a `PD.Mapped` without any option, used for transformer and helper params when the
 * scope has no registered provider to build them from. Applying them fails with `assertRepresentationScope`.
 */
export function emptyMappedParams(): PD.Mapped<any> {
  return PD.MappedStatic('', { '': PD.EmptyGroup() }) as PD.Mapped<any>;
}

/** Values of {@link emptyMappedParams}. */
export function emptyMappedValue(): { name: string; params: any } {
  return { name: '', params: {} };
}

export function hasRepresentations(scope: ScopeRegistries) {
  return scope.registry.default !== undefined;
}

export function hasColorThemes(scope: ScopeRegistries) {
  return scope.themes.colorThemeRegistry.default !== undefined;
}

export function hasSizeThemes(scope: ScopeRegistries) {
  return scope.themes.sizeThemeRegistry.default !== undefined;
}

/** True when the scope has a representation, a color theme, and a size theme to build params from. */
export function hasRepresentationScope(scope: ScopeRegistries) {
  return hasRepresentations(scope) && hasColorThemes(scope) && hasSizeThemes(scope);
}

export function assertRepresentations(scopeName: RepresentationScope, scope: ScopeRegistries) {
  if (!hasRepresentations(scope)) throw new Error(`No ${scopeName} representations are registered in this plugin`);
}

export function assertColorThemes(scopeName: RepresentationScope, scope: ScopeRegistries) {
  if (!hasColorThemes(scope)) throw new Error(`No ${scopeName} color themes are registered in this plugin`);
}

export function assertSizeThemes(scopeName: RepresentationScope, scope: ScopeRegistries) {
  if (!hasSizeThemes(scope)) throw new Error(`No ${scopeName} size themes are registered in this plugin`);
}

/** Throws when the scope has no representation, color theme, or size theme registered. */
export function assertRepresentationScope(scopeName: RepresentationScope, scope: ScopeRegistries) {
  assertRepresentations(scopeName, scope);
  assertColorThemes(scopeName, scope);
  assertSizeThemes(scopeName, scope);
}

export type ProviderKind = 'representation' | 'color theme' | 'size theme';

/** Warns when `name` (if given) is not registered; params are then built as usual and normalization substitutes. */
export function warnIfUnregistered(
  plugin: PluginContext,
  scopeName: RepresentationScope,
  kind: ProviderKind,
  registry: NameRegistry,
  name: string | undefined,
) {
  if (name === undefined || registry.has(name)) return;
  const label = scopeName.charAt(0).toUpperCase() + scopeName.slice(1);
  plugin.log.warn(`${label} ${kind} '${name}' is not registered in this plugin; the registry default is used`);
}

interface RepresentationLike {
  defaultColorTheme: { name: string };
  defaultSizeTheme: { name: string };
}

interface RepresentationScopeRegistries extends ScopeRegistries {
  registry: NameRegistry & {
    get(name: string): RepresentationLike;
    readonly default: { provider: RepresentationLike } | undefined;
  };
}

export interface RequestedNames {
  /** Representation name or provider; the registry default when not given. */
  type?: string | RepresentationLike;
  /** Color theme name or provider; the representation's default color theme when not given. */
  color?: string | object;
  /** Size theme name or provider; the representation's default size theme when not given. */
  size?: string | object;
  /** Whether the call builds a color theme (so its name and the default color theme are checked). */
  checkColor: boolean;
  /** Whether the call builds a size theme (so its name and the default size theme are checked). */
  checkSize: boolean;
}

/**
 * Warns for each requested representation, color theme, and size theme name, and for the representation's default
 * theme names that are used, when it is not registered in the scope. Provider objects are used directly and not
 * checked. Params are built as usual afterwards; normalization substitutes the registry defaults.
 */
export function warnUnregisteredNames(
  plugin: PluginContext,
  scopeName: RepresentationScope,
  scope: RepresentationScopeRegistries,
  request: RequestedNames,
) {
  const { registry, themes } = scope;
  let repr: RepresentationLike | undefined;
  if (typeof request.type === 'string') {
    warnIfUnregistered(plugin, scopeName, 'representation', registry, request.type);
    // an unregistered name resolves to the empty provider, which has no default themes to check
    if (registry.has(request.type)) repr = registry.get(request.type);
  } else {
    repr = request.type ?? registry.default?.provider;
  }

  if (request.checkColor) {
    const color = request.color ? request.color : repr?.defaultColorTheme.name;
    if (typeof color === 'string')
      warnIfUnregistered(plugin, scopeName, 'color theme', themes.colorThemeRegistry, color);
  }
  if (request.checkSize) {
    const size = request.size ? request.size : repr?.defaultSizeTheme.name;
    if (typeof size === 'string') warnIfUnregistered(plugin, scopeName, 'size theme', themes.sizeThemeRegistry, size);
  }
}
