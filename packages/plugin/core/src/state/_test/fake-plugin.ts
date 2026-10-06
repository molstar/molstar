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

export type RegistryKind = 'representations' | 'color-themes' | 'size-themes';

/**
 * A plugin stand-in with real (built-in preloaded) representation and theme registries for every scope and a
 * recording `log.warn`. `empty` lists the registry kinds to clear in every scope.
 */
export function createFakePlugin(empty: ReadonlyArray<RegistryKind> = []) {
  const warnings: string[] = [];
  const scope = <R extends { clear(): void }>(registry: R) => {
    if (empty.includes('representations')) registry.clear();
    const colorThemeRegistry = ColorTheme.createRegistry();
    const sizeThemeRegistry = SizeTheme.createRegistry();
    if (empty.includes('color-themes')) colorThemeRegistry.clear();
    if (empty.includes('size-themes')) sizeThemeRegistry.clear();
    return { registry, themes: { colorThemeRegistry, sizeThemeRegistry } };
  };

  const plugin = {
    representation: {
      structure: scope(new StructureRepresentationRegistry()),
      volume: scope(new VolumeRepresentationRegistry()),
      particles: scope(new ParticleRepresentationRegistry()),
    },
    log: { warn: (msg: string) => warnings.push(msg) },
  } as unknown as PluginContext;

  return { plugin, warnings };
}
