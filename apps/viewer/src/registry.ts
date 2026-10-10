/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultFormats, DefaultRegistry } from '@molstar/plugin/default-registry';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';

/**
 * `DefaultRegistry` followed by the entry of the user's `customFormats`, which is last so the built-in format order is
 * unchanged. As in 5.x, a custom name that equals a built-in format name overrides it: `DefaultFormats` is replaced
 * by a copy without that provider, so `dataFormats.get(name)` returns the custom one and no name is registered twice.
 */
export function createViewerRegistry(
  customFormats: [name: string, provider: DataFormatProvider.Unnamed][] | undefined,
): readonly PluginRegistryEntry[] {
  const formats = (customFormats ?? []).map(([name, provider]) => DataFormatProvider.withName(provider, name));
  const customNames = new Set(formats.map((f) => f.name));
  const overridesDefault = DefaultFormats.formats!.some((f) => customNames.has(f.name));

  const registry = overridesDefault
    ? DefaultRegistry.map((e) =>
        e === DefaultFormats ? { ...e, formats: e.formats!.filter((f) => !customNames.has(f.name)) } : e,
      )
    : DefaultRegistry;
  return [...registry, { formats }];
}
