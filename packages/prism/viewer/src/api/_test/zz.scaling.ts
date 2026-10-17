import { expect, it } from 'vitest';
import { cbs, expectOk, headlessViewer } from '../../_test/harness.js';
it('probe', async () => {
  const out: string[] = [];
  for (const n of [12, 48]) {
    const g = globalThis as any;
    const v = await headlessViewer();
    g.__syncMs = 0;
    g.__syncN = 0;
    const t0 = performance.now();
    for (let i = 0; i < n; i++) expectOk(await v.api.load.files({ files: [cbs()] }));
    out.push(`${n}: total ${(performance.now() - t0).toFixed(0)} ms, sync ${g.__syncMs.toFixed(0)} ms in ${g.__syncN}`);
    v.dispose();
  }
  expect(out).toEqual([]);
}, 600000);
