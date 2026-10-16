/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NavigationManager, ViewStateManager } from '@molstar/prism-ui';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { memoryStore } from '../../persist/stores.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import type { PrismStore, PrismViewerSpec } from '../../spec.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { KeepLayoutDelayMs, keepLayout, LayoutKey, layoutStore, restoreLayout } from '../layout.js';
import { mountShell } from '../shell.js';

// The layout kept in the preferences (journal.md §6, shell.md §3): written a moment after a change, restored at start
// with the host's parts left as the spec has them; another version's value is ignored; an embedded layout is not kept.

class NoResize {
  observe() {}
  unobserve() {}
  disconnect() {}
}

const cores: ViewerCore[] = [];
beforeEach(() => {
  vi.stubGlobal('ResizeObserver', NoResize);
});
afterEach(async () => {
  for (const c of cores.splice(0)) await c.dispose();
  document.body.innerHTML = '';
  vi.useRealTimers();
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

const store = () => memoryStore({ version: 'test', budget: 64 * 1024 });
const managers = (seed?: ConstructorParameters<typeof ViewStateManager>[0]) => ({
  viewState: new ViewStateManager(seed),
  navigation: new NavigationManager(),
});

describe('the kept layout', () => {
  it('is written a moment after a change, only when it changed, and once more when stopped', () => {
    vi.useFakeTimers();
    const prefs = store();
    const m = managers();
    const stop = keepLayout(m, prefs);
    expect(prefs.size).toBe(0);
    m.viewState.setDock('right', { size: 420 });
    m.viewState.setDock('right', { size: 440 });
    vi.advanceTimersByTime(KeepLayoutDelayMs - 1);
    expect(prefs.size).toBe(0);
    vi.advanceTimersByTime(1);
    expect(prefs.size).toBe(1);
    m.viewState.setTabHidden('history', true);
    stop();
    const kept = prefs.get(LayoutKey) as { view: { layout: { hiddenTabs: string[] } } };
    expect(kept.view.layout.hiddenTabs).toEqual(['history']);
  });

  it('restores the docks, the hidden tabs and the docks’ tabs, not the host’s parts, the card or the page', () => {
    const prefs = store();
    const a = managers();
    const stop = keepLayout(a, prefs);
    a.viewState.setDock('right', { size: 420, pinned: true, tabs: ['home', 'gone'] });
    a.viewState.setTabHidden('history', true);
    a.viewState.setButtons({ xr: true });
    a.viewState.setMode('desktop');
    a.navigation.setDock('top', { collapsed: true });
    a.navigation.setCard({ id: 'export' });
    stop();

    const b = managers({ compactRail: true, buttons: { xr: false } });
    expect(restoreLayout(b, prefs)).toBe(true);
    const layout = b.viewState.layout.value;
    expect(layout.docks.right).toMatchObject({ size: 420, pinned: true, tabs: ['home', 'gone'] });
    expect(layout.hiddenTabs).toEqual(['history']);
    expect(layout).toMatchObject({ mode: 'auto', compactRail: true, buttons: { xr: false } });
    expect(b.navigation.state.value.docks.top.collapsed).toBe(true);
    expect(b.navigation.state.value.card).toBeUndefined();
  });

  it('ignores a value of another version or shape', () => {
    const prefs = store();
    const m = managers();
    prefs.set(LayoutKey, { view: { version: 2, layout: { docks: {} } }, navigation: { version: 1 } });
    expect(restoreLayout(m, prefs)).toBe(false);
    prefs.set(LayoutKey, 'nonsense');
    expect(restoreLayout(m, prefs)).toBe(false);
    const throwing: PrismStore = {
      get: () => {
        throw new Error('blocked');
      },
      set: () => {},
      remove: () => {},
    };
    expect(restoreLayout(m, throwing)).toBe(false);
    expect(m.viewState.layout.value).toEqual(m.viewState.initial);
  });

  it('is kept by default, not for an embedded layout or when the spec says not to', () => {
    const preferences = store();
    const spec = (s: Partial<PrismViewerSpec>): PrismViewerSpec => ({ ...DefaultPrismViewerSpec(), ...s });
    expect(layoutStore(spec({}))).toBeUndefined();
    expect(layoutStore(spec({ persist: { key: 'k', preferences } }))).toBe(preferences);
    expect(layoutStore(spec({ persist: { key: 'k', preferences, layout: false } }))).toBeUndefined();
    expect(layoutStore(spec({ persist: { key: 'k', preferences }, ui: { layout: { mode: 'embedded' } } }))).toBe(
      undefined,
    );
    expect(layoutStore(spec({ persist: { key: 'k', preferences }, ui: { layout: { compactRail: true } } }))).toBe(
      undefined,
    );
    expect(
      layoutStore(spec({ persist: { key: 'k', preferences, layout: true }, ui: { layout: { compactRail: true } } })),
    ).toBe(preferences);
  });

  it('comes back in the next viewer on the page', async () => {
    const preferences = store();
    const mount = async () => {
      const c = await createViewerCore({
        // Home at start while nothing is loaded, as the app does.
        spec: {
          ...DefaultPrismViewerSpec(),
          options: { ui: { homeWhenEmpty: true } },
          persist: { key: 'layout-test', preferences },
        },
        createPlugin: createDefaultPlugin,
        modules: ['files', 'page'],
      });
      cores.push(c);
      vi.spyOn(c, 'mountCanvas').mockResolvedValue(false);
      const target = document.createElement('div');
      document.body.append(target);
      return { core: c, mounted: await mountShell(c, target) };
    };

    const first = await mount();
    expect(first.mounted.shell.dock('right').active.value?.id).toBe('home');
    first.mounted.shell.hideTab('home');
    first.mounted.shell.hideTab('history');
    first.mounted.shell.setSize('right', 460);
    first.mounted.dispose();

    const second = await mount();
    const shell = second.mounted.shell;
    expect(shell.dock('right').state.value.size).toBe(460);
    expect(shell.hiddenTabs.value.map((t) => t.id)).toEqual(['home', 'history']);
    expect(shell.dock('right').tabs.value.some((t) => t.id === 'history')).toBe(false);
    // Home, hidden, does not open at start.
    expect(shell.dock('right').active.value).toBeUndefined();

    // Reset layout, as the command does: the tab is back, the size the default.
    await second.core.managers.commands.open('view.reset-layout');
    expect(shell.hiddenTabs.value).toEqual([]);
    expect(shell.dock('right').state.value.size).toBe(340);
    second.mounted.dispose();
  });
});
