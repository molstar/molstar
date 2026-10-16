/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createDefaultPlugin } from '../../sources/create-plugin.js';
import type { PrismStore } from '../../spec.js';
import { createViewerCore, type ViewerCore } from '../../viewer.js';
import { mountShell } from '../shell.js';

// The page's shell over a viewer core in jsdom (no WebGL: the canvas says it could not start). What it keeps in the
// preferences is what the person chose, not what it started with; a shell that fails partway leaves nothing behind.

// jsdom has no ResizeObserver: the shell's layout watchers see no resizes.
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
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

function memoryStore(): PrismStore & { readonly values: Map<string, unknown> } {
  const values = new Map<string, unknown>();
  return {
    values,
    get: (k) => values.get(k),
    set: (k, v) => void values.set(k, v),
    remove: (k) => void values.delete(k),
  };
}

async function core(preferences?: PrismStore) {
  const c = await createViewerCore({
    spec: { ...DefaultPrismViewerSpec(), ...(preferences ? { persist: { key: 'mount-test', preferences } } : {}) },
    createPlugin: createDefaultPlugin,
    modules: ['files', 'page'],
  });
  cores.push(c);
  // No WebGL in jsdom: the canvas does not start (the shell says so).
  vi.spyOn(c, 'mountCanvas').mockResolvedValue(false);
  return c;
}

describe('the shell on a page', () => {
  it('keeps the theme only once the person changes it', async () => {
    const preferences = memoryStore();
    const c = await core(preferences);
    const target = document.createElement('div');
    document.body.append(target);
    const mounted = await mountShell(c, target);
    expect(preferences.values.has('theme')).toBe(false);
    c.managers.commands.open('view.theme');
    await new Promise((r) => setTimeout(r, 0));
    expect(preferences.values.has('theme')).toBe(true);
    mounted.dispose();
  });

  it('removes what it made when it fails partway', async () => {
    const c = await core();
    const target = document.createElement('div');
    document.body.append(target);
    vi.spyOn(c, 'attachExtensionUI').mockImplementation(() => {
      throw new Error('no port');
    });
    await expect(mountShell(c, target)).rejects.toThrow('no port');
    // Its root is gone from the shadow root, and its tabs with it (the shell retired).
    expect(target.shadowRoot?.childElementCount ?? 0).toBe(0);
  });
});
