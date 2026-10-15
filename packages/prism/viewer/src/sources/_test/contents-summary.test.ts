/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { cbs, headlessViewer, until } from '../../_test/harness.js';

// A structure's contents in its inspector: its counts, its entities (kind, chains, atoms), each opened in the tree.

describe('the contents summary', () => {
  it('counts the structure, lists its entities, and opens one in the tree', async () => {
    const v = await headlessViewer({ files: [cbs()] });
    await until(() => v.views.scene.roots.value.length > 0);
    const contents = v.views.scene.roots.value[0].children.value.find((n) => n.kind === 'contents')!;
    const summary = v.views.scene.inspect(contents).contents!;
    expect(summary.counts.map((c) => `${c.label} ${c.value}`)).toEqual([
      'Entities 3',
      'Chains 3',
      'Residues 137',
      'Atoms 1,213',
    ]);
    expect(summary.entities.map((e) => [e.type, e.chains.join(), e.atoms])).toEqual([
      ['polymer', 'A', 1091],
      ['non-polymer', 'A', 22],
      ['water', 'A', 100],
    ]);
    summary.reveal(summary.entities[1].id);
    expect(contents.isExpanded.value).toBe(true);
    expect(v.views.scene.picked.value?.label.value).toBe('RETINOIC ACID');
  });
});
