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

// The viewer's bulk paths over a real plugin grow about linearly in the number of structures (`test:scaling`, prism-ui's
// _test/scaling.ts): loading them one after another, a preset applied to every one, a select of all and its clear.
// The time is Mol*'s own work and the viewer's (the scene's sync, the sources, the journal): a viewer pass over every
// structure per structure would show here as a quadratic.

interface Loaded {
  v: PrismViewer;
  n: number;
}

async function loaded(n: number): Promise<Loaded> {
  const v = await headlessViewer();
  for (let i = 0; i < n; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
  return { v, n };
}

it('loading N structures, one after another', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => ({ v: await headlessViewer(), n }),
      run: async ({ v, n }) => {
        for (let i = 0; i < n; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
      },
      teardown: ({ v }) => v.dispose(),
    },
    12,
  );
  checkLinear('load N structures', s);
}, 600_000);

it('a preset applied to N structures', async () => {
  const s = await scalingAsync(
    {
      setup: loaded,
      run: async ({ v }) => {
        expectOk(await v.api.structure.preset({ preset: 'illustrative' }));
      },
      teardown: ({ v }) => v.dispose(),
    },
    3,
  );
  checkLinear('preset over N structures', s);
}, 600_000);

it('a select of everything in N structures, and its clear', async () => {
  const s = await scalingAsync(
    {
      setup: loaded,
      run: async ({ v }) => {
        expectOk(await v.api.select.set({ target: 'polymer' }));
        expectOk(await v.api.select.clear());
      },
      teardown: ({ v }) => v.dispose(),
    },
    3,
  );
  checkLinear('select all and clear, N structures', s);
}, 600_000);
