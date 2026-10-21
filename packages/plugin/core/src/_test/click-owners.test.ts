/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import type { Representation } from '@molstar/graphics/repr/representation';
import { PluginContext } from '@molstar/plugin/context';

describe('PluginContext.clickOwners', () => {
  it("tells the representations a tool claimed (their clicks are the tool's), until it lets go", () => {
    const plugin = new PluginContext({ behaviors: [], registry: [] });
    const gizmo = {} as Representation.Any;
    const molecule = {} as Representation.Any;
    const release = plugin.clickOwners.claim((r) => r === gizmo);
    expect(plugin.clickOwners.owns(gizmo)).toBe(true);
    expect(plugin.clickOwners.owns(molecule)).toBe(false);
    expect(plugin.clickOwners.owns(undefined)).toBe(false);
    release();
    expect(plugin.clickOwners.owns(gizmo)).toBe(false);
    plugin.dispose();
  });
});
