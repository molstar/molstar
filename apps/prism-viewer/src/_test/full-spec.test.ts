/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import { isLazyExtension } from '@molstar/prism-viewer/spec';
import { describe, expect, it } from 'vitest';
import { FullPrismViewerSpec } from '../full-spec.js';

// The app's whole composition (app.md §2): every extension at once, headless. Ids that collide, a stub that does not
// match what it loads, a network operation that does not say what it reaches, or an extension's operation without a
// command would show here first.

describe('FullPrismViewerSpec', () => {
  it('composes headless, and every lazy extension loads as its stub says', async () => {
    const spec = FullPrismViewerSpec();
    const ids = (spec.extensions ?? []).map((e) => e.id);
    expect(new Set(ids).size).toBe(ids.length);
    const v = await createPrismViewerHeadless({ spec });
    try {
      for (const e of spec.extensions ?? []) {
        // One behind an option that is off (remote states) is not in the viewer.
        if (isLazyExtension(e) && v.extensions.status(e.id)) expect(await v.extensions.load(e.id), e.id).toBe(true);
      }
      // A network operation names what it reaches, or the agent always asks: only the generic ones may not say.
      const silent = [...v.api.definitions.value.values()]
        .filter((d) => d.approval === 'network' && !d.addresses)
        .map((d) => d.id);
      expect(silent.sort()).toEqual(['load.action', 'load.svs', 'state.action']);
      // Coverage (api.md §10), as the viewer's own (coverage.test.ts) but with the extensions: each operation has a
      // command or a tool, or says why it has neither (`surface`); reads are the agent's read tools.
      const uncommanded = [...v.api.definitions.value.values()]
        .filter((d) => !d.read && !d.surface && !v.managers.commands.get(d.id))
        .map((d) => d.id);
      expect(uncommanded.sort()).toEqual([]);
    } finally {
      await v.dispose();
    }
  });
});
