/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Representation } from '@molstar/graphics/repr/representation';
import { PluginContext } from '@molstar/plugin/context';

describe('InteractivityManager.clickOwners', () => {
  it("tells the representations a tool claimed (their clicks are the tool's), until it lets go", async () => {
    const plugin = new PluginContext({ behaviors: [], registry: [] });
    await plugin.init();
    const { clickOwners } = plugin.managers.interactivity;
    const gizmo = {} as Representation.Any;
    const molecule = {} as Representation.Any;
    const release = clickOwners.claim((r) => r === gizmo);
    expect(clickOwners.owns(gizmo)).toBe(true);
    expect(clickOwners.owns(molecule)).toBe(false);
    expect(clickOwners.owns(undefined)).toBe(false);
    release();
    expect(clickOwners.owns(gizmo)).toBe(false);
    plugin.dispose();
  });
});
