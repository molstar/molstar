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
// _test/scaling.ts), 12 and 48: loading them (one call, and one after another), a preset over every one, a select of
// all and its clear, a component made in every one, every representation updated. The time is Mol*'s own work and the
// viewer's (the scene's sync, the sources, the journal): a viewer pass over every structure per structure (a commit, a
// sync, a name made per structure) shows here as a quadratic.

const N = 12;

interface Loaded {
  v: PrismViewer;
  n: number;
}

async function loaded(n: number): Promise<Loaded> {
  const v = await headlessViewer();
  expectOk(await v.api.load.files({ files: Array.from({ length: n }, () => cbs()) }));
  return { v, n };
}
const dispose = ({ v }: Loaded) => v.dispose();

it('loading N structures in one call', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => ({ v: await headlessViewer(), n }),
      run: async ({ v, n }) => {
        expectOk(await v.api.load.files({ files: Array.from({ length: n }, () => cbs()) }));
      },
      teardown: dispose,
    },
    N,
  );
  checkLinear('load N structures, one call', s);
}, 600_000);

it('loading N structures, one after another', async () => {
  const s = await scalingAsync(
    {
      setup: async (n) => ({ v: await headlessViewer(), n }),
      run: async ({ v, n }) => {
        for (let i = 0; i < n; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
      },
      teardown: dispose,
    },
    N,
  );
  checkLinear('load N structures, one by one', s);
}, 600_000);

it('a preset applied to N structures', async () => {
  const s = await scalingAsync(
    {
      setup: loaded,
      run: async ({ v }) => {
        expectOk(await v.api.structure.preset({ preset: 'illustrative' }));
      },
      teardown: dispose,
    },
    N,
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
      teardown: dispose,
    },
    N,
  );
  checkLinear('select all and clear, N structures', s);
}, 600_000);

it('a component made in each of N structures', async () => {
  const s = await scalingAsync(
    {
      setup: loaded,
      run: async ({ v }) => {
        expectOk(await v.api.component.add({ target: 'ligand', label: 'Ligands' } as never));
      },
      teardown: dispose,
    },
    N,
  );
  checkLinear('component in N structures', s);
}, 600_000);

it('every representation of N structures updated', async () => {
  const s = await scalingAsync(
    {
      setup: loaded,
      run: async ({ v }) => {
        expectOk(
          await v.api.representation.update({
            node: 'representation',
            params: { type: { params: { alpha: 0.5 } } },
          } as never),
        );
      },
      teardown: dispose,
    },
    N,
  );
  checkLinear('representations of N structures updated', s);
}, 600_000);
