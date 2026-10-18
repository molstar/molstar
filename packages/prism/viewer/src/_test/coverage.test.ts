/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultPrismViewerSpec } from '../default-spec.js';
import { isLazyExtension } from '../spec.js';
import { headlessViewer } from './harness.js';

// Coverage (api.md §10): every operation of the default viewer, the lazy features loaded, has a command or a tool,
// or says why it has neither (`surface`). Reads (`def.read`) count as covered: they are the agent's read tools
// (agent.md §3).

describe('coverage', () => {
  it('every operation has a command, or says why not', async () => {
    // Every feature, the ones behind an option too.
    const v = await headlessViewer({ spec: { options: { ui: { showRemoteStates: true } } } });
    for (const e of DefaultPrismViewerSpec().extensions ?? []) {
      if (isLazyExtension(e)) expect(await v.extensions.load(e.id)).toBe(true);
    }
    const missing = [...v.api.definitions.value.values()]
      .filter((def) => !def.read && !def.surface && !v.managers.commands.get(def.id))
      .map((def) => def.id)
      .sort();
    expect(missing).toEqual([]);
  });
});
