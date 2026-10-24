/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { it } from 'vitest';
import { checkLinear, type Scaling, scalingAsync } from '../../../../ui/src/_test/scaling.js';
import type { Observable } from 'rxjs';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
import { createSceneSource, type SceneSource } from '../scene.js';
import { SourceScope } from '../state.js';

// The scene source's sync is incremental (`test:scaling`, prism-ui's _test/scaling.ts; sources/scene.ts `syncRows`):
// a structure whose hierarchy refs are as they were keeps its rows, so loading N structures one after another is N
// builds of the changed one, not N syncs of every one. Timed: a scene source of the test's own, between a listener
// before it and one after it on each event it syncs at (a subject tells its listeners in order), so neither Mol*'s
// load nor the viewer's other listeners (its own scene source's, which reads the hierarchies first) are in the time.
// At 12 and 48.

const N = 12;

/** The time a scene source takes to sync while `n` structures are loaded one by one. */
async function syncsWhileLoading(n: number): Promise<number> {
  const v = await headlessViewer();
  const p = v.plugin;
  const events = [
    p.managers.structure.hierarchy.behaviors.selection,
    p.managers.volume.hierarchy.behaviors.selection,
    p.managers.particles.hierarchy.behaviors.selection,
    p.managers.structure.measurement.behaviors.state,
  ] as const;
  const scope = new SourceScope();
  let ms = 0;
  let t0 = 0;
  for (const e of events) scope.add((e as Observable<unknown>).subscribe(() => void (t0 = performance.now())));
  createSceneSource(scope, p, v.api);
  for (const e of events) scope.add((e as Observable<unknown>).subscribe(() => void (ms += performance.now() - t0)));
  try {
    for (let i = 0; i < n; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
  } finally {
    scope.dispose();
    await v.dispose();
  }
  return ms;
}

async function best(n: number, reps: number): Promise<number> {
  let min = Number.POSITIVE_INFINITY;
  for (let r = 0; r < reps; r++) min = Math.min(min, await syncsWhileLoading(n));
  return min;
}

it('the scene syncs while N structures load one by one', async () => {
  await best(N, 1);
  await best(4 * N, 1);
  const a = await best(N, 3);
  const b = await best(4 * N, 3);
  const s: Scaling = { n: N, factor: 4, ms: [a, b], bytes: [1, 1], timeRatio: b / Math.max(a, 1e-3), allocRatio: 1 };
  checkLinear('scene syncs, N structures loaded one by one', s);
}, 600_000);

it('the scene syncs N structures, each with part of its contents hidden', async () => {
  // Each structure's polymer hidden in Entities: its Polymer component emptied, its own params saved. A sync reads the
  // emptied components once for every structure (`emptiedComponentsOf`), not every structure's saved params for each.
  const s = await scalingAsync(
    {
      setup: async (n) => {
        const v = await headlessViewer();
        expectOk(await v.api.load.files({ files: Array.from({ length: n }, () => cbs()) }));
        for (const st of v.plugin.managers.structure.hierarchy.current.structures) {
          expectOk(
            await v.api.run('contents.visibility', { structure: st.cell.transform.ref, entity: '1', visible: false }),
          );
        }
        return v;
      },
      run: async (v) => {
        for (let i = 0; i < 100; i++) (v.views.scene as SceneSource).sync();
      },
      teardown: (v) => v.dispose(),
    },
    N,
  );
  checkLinear('scene syncs, N structures with hidden contents', s);
}, 600_000);
