/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';

/**
 * Development-mode check that every registered representation has its default color and size themes registered in
 * its scope. It warns through `plugin.log.warn` and never throws, registers, or resolves anything; this module must
 * not value-import catalogs.
 */
export { checkDefaultThemes };

const Scopes = [
  ['Structure', 'structure'],
  ['Volume', 'volume'],
  ['Particles', 'particles'],
] as const;

/**
 * Warns once per (scope, representation, theme kind, theme name) for the lifetime of `warned` (one set per plugin).
 * Returns the messages it logged.
 */
function checkDefaultThemes(plugin: PluginContext, warned: Set<string>): string[] {
  const messages: string[] = [];
  for (const [label, scope] of Scopes) {
    const { registry, themes } = plugin.representation[scope];
    for (const { name, provider } of registry.list) {
      const pairs = [
        ['color', provider.defaultColorTheme?.name, themes.colorThemeRegistry],
        ['size', provider.defaultSizeTheme?.name, themes.sizeThemeRegistry],
      ] as const;
      for (const [kind, themeName, themeRegistry] of pairs) {
        // Providers without default theme names (hand-written test doubles) have nothing to check.
        if (!themeName || themeRegistry.has(themeName)) continue;
        const key = `${scope}/${name}/${kind}/${themeName}`;
        if (warned.has(key)) continue;
        warned.add(key);
        const message = `${label} representation '${name}' default ${kind} theme '${themeName}' is not registered`;
        messages.push(message);
        plugin.log.warn(message);
      }
    }
  }
  return messages;
}
