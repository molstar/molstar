/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { StructureRepresentationRegistry } from '@molstar/graphics/repr/structure/registry';
import { VolumeRepresentationRegistry } from '@molstar/graphics/repr/volume/registry';
import { ParticleRepresentationRegistry } from '@molstar/graphics/repr/particles/registry';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { BuiltInStructureRepresentations } from '@molstar/graphics/repr/structure/catalog';
import { BuiltInVolumeRepresentations } from '@molstar/graphics/repr/volume/catalog';
import { BuiltInParticleRepresentations } from '@molstar/graphics/repr/particles/catalog';

export type RegistryKind = 'representations' | 'color-themes' | 'size-themes';

/**
 * A plugin stand-in with real representation and theme registries for every scope, filled with every built-in
 * provider (the registries themselves start empty), and a recording `log.warn`. `empty` lists the registry kinds to
 * leave empty in every scope.
 */
export function createFakePlugin(empty: ReadonlyArray<RegistryKind> = []) {
  const warnings: string[] = [];
  const scope = <R extends { add(p: any): void }>(registry: R, representations: object) => {
    const colorThemeRegistry = ColorTheme.createRegistry();
    const sizeThemeRegistry = SizeTheme.createRegistry();
    if (!empty.includes('representations')) for (const p of Object.values(representations)) registry.add(p);
    if (!empty.includes('color-themes'))
      for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as any);
    if (!empty.includes('size-themes'))
      for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as any);
    return { registry, themes: { colorThemeRegistry, sizeThemeRegistry } };
  };

  const plugin = {
    representation: {
      structure: scope(new StructureRepresentationRegistry(), BuiltInStructureRepresentations),
      volume: scope(new VolumeRepresentationRegistry(), BuiltInVolumeRepresentations),
      particles: scope(new ParticleRepresentationRegistry(), BuiltInParticleRepresentations),
    },
    log: { warn: (msg: string) => warnings.push(msg) },
  } as unknown as PluginContext;

  return { plugin, warnings };
}
