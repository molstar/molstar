/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import '@molstar/plugin/state/transforms/catalog';
import { StateTransformer } from '@molstar/core/state';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { paramsToMd } from '@molstar/state-docs-cli/pd-to-md';

// Mirrors cli/state-docs/src/index.ts: the parameters of every transformer are documented with the default plugin
describe('state-docs with the default plugin', () => {
  it('documents the parameters of every transformer, listing every registered provider', async () => {
    const ctx = new PluginContext(DefaultPluginSpec());
    await ctx.init();

    const transformers = StateTransformer.getAll();
    expect(transformers.length).toBeGreaterThan(100);

    let md = '';
    for (const t of transformers) {
      if (!t.definition.params) continue;
      md += paramsToMd(t.definition.params(undefined, ctx));
    }

    for (const scope of ['structure', 'volume', 'particles'] as const) {
      const { registry, themes } = ctx.representation[scope];
      const names = [
        ...registry.list.map((p) => p.name),
        ...themes.colorThemeRegistry.list.map((p) => p.name),
        ...themes.sizeThemeRegistry.list.map((p) => p.name),
      ];
      expect(names.length).toBeGreaterThan(0);
      for (const name of names) expect(md).toContain(`**${name}**`);
    }
    ctx.dispose();
  });
});
