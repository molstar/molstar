/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { it } from 'vitest';
import { checkLinear, scalingAsync } from '../../../../../ui/src/_test/scaling.js';
import { expectOk, headlessViewer } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { SnapshotsFeature } from '../index.js';

// `snapshot.move` across the list is one step (`test:scaling`, prism-ui's _test/scaling.ts): the first of N snapshots
// moved to the end. Moved by the manager's swaps of neighbors, it was N swaps, each a search of the list and a change
// told. At 250 and 1,000.

it('the first of N snapshots moved to the end', async () => {
  const s = await scalingAsync<{ v: PrismViewer; n: number; first: string }>(
    {
      setup: async (n) => {
        const v = await headlessViewer({ extensions: [SnapshotsFeature] });
        const manager = v.plugin.managers.snapshot;
        for (let i = 0; i < n; i++) {
          manager.add(
            PluginStateSnapshotManager.Entry(v.plugin.state.getSnapshot({ camera: false, image: false }), {
              name: `S${i}`,
            }),
          );
        }
        return { v, n, first: manager.state.entries.get(0)!.snapshot.id };
      },
      run: async ({ v, n, first }) => {
        expectOk(await v.api.run('snapshot.move', { entry: first, position: n }));
      },
      teardown: ({ v }) => v.dispose(),
    },
    250,
  );
  checkLinear('snapshot.move across N', s);
}, 600_000);
