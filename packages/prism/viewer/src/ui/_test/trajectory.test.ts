/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ensemble } from '../../api/_test/trajectory-fixtures.js';
import { cbs, expectOk, headlessViewer, until } from '../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewerHeadless } from '../../headless.js';
import { createTrajectoryControls } from '../trajectory.js';

// The trajectory bar's model over a headless viewer: what it reads, and that it acts through the API.

describe('the trajectory controls', () => {
  it('follow the structure with several models, and seek through the API', async () => {
    const v = await headlessViewer({ files: [cbs(), ensemble()] });
    const { model, dispose } = createTrajectoryControls(v.plugin, v.api.as({ kind: 'person' }));
    try {
      expect(model.count.value).toBe(3);
      expect(model.index.value).toBe(0);
      expect(model.label.value).toMatch(/ensemble/);
      expect(model.isPlaying.value).toBe(false);

      model.seek(1);
      model.seek(2);
      await until(() => model.index.value === 2);
      // Viewing: no step of its own.
      expect(v.managers.history.steps.value.length).toBe(2);
    } finally {
      dispose();
    }
  });

  it('show nothing for structures with one model, and play and pause with a canvas', async () => {
    const plain = await headlessViewer({ files: [cbs()] });
    const none = createTrajectoryControls(plain.plugin, plain.api.as({ kind: 'person' }));
    expect(none.model.count.value).toBe(0);
    none.dispose();

    const v = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await v.api.load.files({ files: [ensemble()] }));
      const { model, dispose } = createTrajectoryControls(v.plugin, v.api.as({ kind: 'person' }));
      model.togglePlay();
      await until(() => model.isPlaying.value);
      model.togglePlay();
      await until(() => !model.isPlaying.value);
      dispose();
    } finally {
      await v.dispose();
    }
  });
});
