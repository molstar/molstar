/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { LociLabelManager, type LociLabelProvider } from '../loci-label.js';

function create() {
  const next = jest.fn();
  const addProvider = jest.fn();
  const plugin = {
    behaviors: { labels: { highlight: { next } } },
    managers: { interactivity: { lociHighlights: { addProvider } } },
  } as unknown as PluginContext;
  return { manager: new LociLabelManager(plugin) };
}

function provider(priority?: number): LociLabelProvider {
  return { label: () => 'label', priority };
}

describe('LociLabelManager', () => {
  it('registers a provider once and counts repeated registration', () => {
    const { manager } = create();
    const p = provider();
    manager.addProvider(p);
    manager.addProvider(p);
    expect(manager.providers).toEqual([p]);

    manager.removeProvider(p);
    expect(manager.providers).toEqual([p]);
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);
  });

  it('ignores removal of an unknown provider', () => {
    const { manager } = create();
    const p = provider();
    manager.removeProvider(p);
    manager.addProvider(p);
    manager.removeProvider(provider());
    expect(manager.providers).toEqual([p]);
    manager.removeProvider(p);
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);

    // a removed provider does not keep a negative count
    manager.addProvider(p);
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);
  });

  it('keeps providers sorted by priority', () => {
    const { manager } = create();
    const low = provider(1);
    const high = provider(10);
    const none = provider();
    manager.addProvider(low);
    manager.addProvider(none);
    manager.addProvider(high);
    expect(manager.providers).toEqual([high, low, none]);
  });

  it('findConflict never reports a conflict and does not count', () => {
    const { manager } = create();
    const p = provider();
    expect(manager.findConflict(p)).toBeUndefined();
    manager.addProvider(p);
    expect(manager.findConflict(p)).toBeUndefined();
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);
  });

  it('clearProviders drops providers and counts', () => {
    const { manager } = create();
    const p = provider();
    manager.addProvider(p);
    manager.addProvider(p);
    manager.clearProviders();
    expect(manager.providers).toEqual([]);

    // later removals of a cleared provider are no-ops
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);

    // and registering it again starts from a fresh count
    manager.addProvider(p);
    manager.removeProvider(p);
    expect(manager.providers).toEqual([]);
  });
});
