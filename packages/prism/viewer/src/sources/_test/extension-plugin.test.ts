/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginConfigItem } from '@molstar/plugin/config';
import { PluginSpec } from '@molstar/plugin/spec';
import { describe, expect, it } from 'vitest';
import { headlessViewer } from '../../_test/harness.js';
import { applyExtensionPlugin } from '../extension-plugin.js';

// A lazy extension's plugin part (viewer-spec.md §2.1): applied whole or not at all, and its undo puts the config back.

const make = (name: string, category: 'custom-props' | 'misc') =>
  PluginBehavior.create({
    name,
    category,
    display: { name },
    ctor: class extends PluginBehavior.Handler {
      register() {}
    },
  });
const Props = make('prism-test-ext-props', 'custom-props');
const Quiet = make('prism-test-ext-quiet', 'misc');
const Setting = new PluginConfigItem<string>('prism-test-ext-setting', 'default');

describe('applyExtensionPlugin', () => {
  it('changes nothing when a behavior is there already', async () => {
    const v = await headlessViewer();
    const { plugin } = v;
    const undoQuiet = await applyExtensionPlugin(plugin, { behaviors: [PluginSpec.Behavior(Quiet)] });
    await expect(
      applyExtensionPlugin(plugin, {
        config: [[Setting, 'set']],
        behaviors: [PluginSpec.Behavior(Props), PluginSpec.Behavior(Quiet)],
      }),
    ).rejects.toThrow(/there already/);
    expect(plugin.state.behaviors.transforms.has(Props.id)).toBe(false);
    expect(plugin.config.get(Setting)).toBe('default');
    await undoQuiet();
    expect(plugin.state.behaviors.transforms.has(Quiet.id)).toBe(false);
  });

  it('puts its config back in the undo, unless something set it after', async () => {
    const v = await headlessViewer();
    const { plugin } = v;
    const undo = await applyExtensionPlugin(plugin, {
      config: [[Setting, 'mine']],
      behaviors: [PluginSpec.Behavior(Props)],
    });
    expect(plugin.config.get(Setting)).toBe('mine');
    expect(plugin.state.behaviors.transforms.has(Props.id)).toBe(true);
    await undo();
    expect(plugin.config.get(Setting)).toBe('default');
    expect(plugin.state.behaviors.transforms.has(Props.id)).toBe(false);

    const again = await applyExtensionPlugin(plugin, { config: [[Setting, 'mine']] });
    plugin.config.set(Setting, 'later');
    await again();
    expect(plugin.config.get(Setting)).toBe('later');
  });
});
