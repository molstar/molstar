/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { setCanvasModule } from '@molstar/graphics/geo/geometry/text/font-atlas';
import { beforeAll, it } from 'vitest';
import { checkLinear, scalingAsync } from '../../../../../ui/src/_test/scaling.js';
import { cbs, expectOk, headlessViewer } from '../../../_test/harness.js';
import type { PrismViewer } from '../../../viewer.js';
import { MeasurementsFeature } from '../index.js';

// Measurements added one by one grow about linearly (`test:scaling`, prism-ui's _test/scaling.ts): each operation (its
// identity check, its read back) and each scene sync read the measurements, and a measurement is described (its atoms
// labelled) once, not at every read (sources/measurements.ts `readMeasurements`). At 12 and 48.

const N = 12;

// Measurement texts are drawn into a font atlas, which Node has no canvas for: a stand-in context that draws nothing.
beforeAll(() => {
  const context = {
    measureText: () => ({ width: 8 }),
    clearRect() {},
    fillText() {},
    getImageData: (_x: number, _y: number, w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
  };
  setCanvasModule({ createCanvas: () => ({ getContext: () => context }) });
});

it('N distances measured one by one', async () => {
  const s = await scalingAsync<{ v: PrismViewer; n: number }>(
    {
      setup: async (n) => ({ v: await headlessViewer({ files: [cbs()], extensions: [MeasurementsFeature] }), n }),
      run: async ({ v, n }) => {
        for (let i = 0; i < n; i++) {
          expectOk(await v.api.run('measure.distance', { a: `A:${i + 2}@CA`, b: `A:${i + 3}@CA` }));
        }
      },
      teardown: ({ v }) => v.dispose(),
    },
    N,
  );
  checkLinear('measure N distances one by one', s);
}, 600_000);
