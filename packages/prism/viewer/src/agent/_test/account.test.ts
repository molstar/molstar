/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer } from '../../_test/harness.js';
import { createAccount } from '../account.js';

// The scene's outline at its levels (agent.md §5): each turn the components', get_scene the full one by default.

describe('the scene outline', () => {
  it('says less at the lower levels: summary, components, full', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    const account = createAccount(v);
    const summary = account.outline(40, 'summary');
    const components = account.outline(40, 'components');
    const full = account.outline(40, 'full');
    expect(summary).toMatch(/^- structure "1CBS": "1,213 atoms"; components "Polymer", "Ligand", "Water"$/m);
    expect(summary).not.toMatch(/component "Polymer"/);
    expect(components).toMatch(/component "Polymer": "Cartoon · 1,091 atoms"/);
    expect(components).not.toMatch(/- representation|Unit cell|"Focus"/);
    expect(full).toMatch(/representation "Cartoon"/);
    expect(full).toMatch(/unit-cell "Unit cell"/);
    expect(summary.length).toBeLessThan(components.length);
    expect(components.length).toBeLessThan(full.length);
  });
});
