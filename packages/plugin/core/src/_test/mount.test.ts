/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @vitest-environment jsdom

import { afterEach, describe, expect, it, vi } from 'vitest';
import { Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { PluginContext } from '@molstar/plugin/context';

// Mounting while the graphics context is made asynchronously (the backend is chosen first): `Canvas3DContext.create`
// is a promise the test settles by hand.

function deferredCreate() {
  const context = { dispose: vi.fn() } as unknown as Canvas3DContext;
  let settle!: (ok: boolean) => void;
  const create = vi.spyOn(Canvas3DContext, 'create').mockImplementation(
    () =>
      new Promise<Canvas3DContext>((resolve, reject) => {
        settle = (ok) => (ok ? resolve(context) : reject(new Error('no graphics context')));
      }),
  );
  return { create, context, settle: (ok: boolean) => settle(ok) };
}

async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [] });
  await plugin.init();
  // a failed creation also rejects `canvas3dInitialized`
  plugin.canvas3dInitialized.catch(() => {});
  return plugin;
}

describe('mounting while the graphics context is made', () => {
  afterEach(() => vi.restoreAllMocks());

  it('frees the context and builds nothing when disposed meanwhile', async () => {
    const { create, context, settle } = deferredCreate();
    const plugin = await createPlugin();
    const mounted = plugin.mountAsync(document.createElement('div'));
    await vi.waitFor(() => expect(create).toHaveBeenCalled());
    plugin.dispose();
    settle(true);
    expect(await mounted).toBe(false);
    expect(context.dispose).toHaveBeenCalledTimes(1);
    expect(plugin.canvas3dContext).toBeUndefined();
    expect(plugin.canvas3d).toBeUndefined();
  });

  it('makes one container and context for mounts that overlap', async () => {
    const { create, settle } = deferredCreate();
    const plugin = await createPlugin();
    const target = document.createElement('div');
    const first = plugin.mountAsync(target);
    // e.g. React StrictMode: mounted, unmounted (nothing to unmount yet), mounted again
    plugin.unmount();
    const second = plugin.mountAsync(target);
    await vi.waitFor(() => expect(create).toHaveBeenCalled());
    settle(false);
    expect(await Promise.all([first, second])).toEqual([false, false]);
    expect(create).toHaveBeenCalledTimes(1);
    expect(target.childElementCount).toBe(0);
    plugin.dispose();
  });
});
