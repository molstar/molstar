/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { LoadVolseg, Volseg } from '@molstar/volumes-and-segmentations-extension';

const actionCount = (plugin: PluginContext) => {
  const actions = (plugin.state.data.actions as any).actions as Map<string, { count: number }>;
  return actions.get(LoadVolseg.id)?.count ?? 0;
};

describe('Volseg behavior', () => {
  const fetch = globalThis.fetch;

  beforeAll(() => {
    // the behavior asks the volume server for the list of entries when it registers
    globalThis.fetch = (async () => ({ ok: true, json: async () => ({}) })) as any;
  });

  afterAll(() => {
    globalThis.fetch = fetch;
  });

  it('registers its action through a registry entry and counts it', async () => {
    const defaults = DefaultPluginSpec();
    const plugin = new PluginContext({
      ...defaults,
      behaviors: [...defaults.behaviors, { transformer: Volseg, defaultParams: {} } as any],
    });
    await plugin.init();
    expect(actionCount(plugin)).toBe(1);

    // an app that lists the same action in its own entry
    const undo = plugin.register({ actions: [LoadVolseg] });
    expect(actionCount(plugin)).toBe(2);

    await plugin.runTask(plugin.state.behaviors.updateTree(plugin.state.behaviors.build().delete(Volseg.id)));
    expect(actionCount(plugin)).toBe(1);

    await plugin.state.updateBehavior(Volseg, () => {});
    expect(actionCount(plugin)).toBe(2);

    undo();
    expect(actionCount(plugin)).toBe(1);
    await plugin.runTask(plugin.state.behaviors.updateTree(plugin.state.behaviors.build().delete(Volseg.id)));
    expect(actionCount(plugin)).toBe(0);
    plugin.dispose();
  });
});
