/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { exampleFile, headlessViewer } from '../../_test/harness.js';
import { evalReport, runEval } from '../eval.js';
import { FirstScope } from './first-scope.js';
import { env, LiveModels } from './live-models.js';

// The first eval on real models (`pnpm test:live`, never in `test:unit`): each first-scope task on a headless viewer,
// checked on the scene; the report (checks, tools, tokens) is printed.

for (const m of LiveModels) {
  describe.skipIf(!env[m.key])(`the first eval on ${m.model}`, () => {
    for (const task of FirstScope) {
      it(task.title, async () => {
        const v = await headlessViewer();
        const provider = m.create();
        const r = await runEval(v, task, {
          models: {
            providerFor: (p) => (p === m.provider ? provider : undefined),
            default: { provider: m.provider, model: m.model },
          },
          file: (name) => exampleFile(name),
        });
        console.info(`${m.model}: ${evalReport(r)}`);
        expect(r.pass, evalReport(r)).toBe(true);
      });
    }
  });
}
