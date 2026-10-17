/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { it } from 'vitest';
import { checkLinear, scalingAsync } from '../../../../ui/src/_test/scaling.js';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import type { PrismViewer } from '../../viewer.js';

// The representation operations over many cells are one state update (`test:scaling`, prism-ui's _test/scaling.ts;
// packages/prism/AGENTS.md, bulk operations): `representation.material` over every representation of N structures,
// `representation.set` over N representations of its type on one component. One commit each was N commits, each a walk
// of the tree and a sync of the viewer.

const N = 12;
// One component holds N representations of a type: at 50 and 200 (Mol*'s work for each is small; the commits were not).
const M = 50;

it('representation.material over N structures', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const v = await headlessViewer();
        expectOk(await v.api.load.files({ files: Array.from({ length: n }, () => cbs()) }));
        return v;
      },
      run: async (v) => {
        expectOk(await v.api.run('representation.material', { material: { metalness: 0.5, roughness: 0.2 } }));
      },
      teardown: (v) => v.dispose(),
    },
    N,
  );
  checkLinear('representation.material, N structures', s);
}, 600_000);

it('representation.set over N representations of its type on a component', async () => {
  const s = await scalingAsync<PrismViewer>(
    {
      setup: async (n) => {
        const v = await headlessViewer({ files: [cbs()] });
        const polymer = v.plugin.managers.structure.hierarchy.current.structures[0].components.find(
          (c) => c.cell.obj?.label === 'Polymer',
        )!;
        const have = polymer.representations.length;
        for (let i = have; i < n; i++) {
          await v.plugin.builders.structure.representation.addRepresentation(polymer.cell, { type: 'cartoon' });
        }
        return v;
      },
      run: async (v) => {
        expectOk(
          await v.api.run('representation.set', { component: 'Polymer', type: 'cartoon', params: { alpha: 0.5 } }),
        );
      },
      teardown: (v) => v.dispose(),
    },
    M,
  );
  checkLinear('representation.set, N of its type', s);
}, 600_000);
