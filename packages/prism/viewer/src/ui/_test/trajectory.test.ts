/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
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

  it("follow a particle trajectory's frames (Simularium) as a structure's, through particles.*", async () => {
    const file = new File(
      [
        fs.readFileSync(
          path.resolve(import.meta.dirname, '../../../../../../data/examples/particles/small.simularium'),
        ),
      ],
      'small.simularium',
    );
    const v = await headlessViewer({ files: [file] });
    const { model, dispose } = createTrajectoryControls(v.plugin, v.api.as({ kind: 'person' }));
    try {
      expect(model.count.value).toBe(5);
      expect(model.unit.value).toBe('frames');
      model.seek(3);
      await until(() => model.index.value === 3);
      const list = v.plugin.managers.particles.hierarchy.current.lists[0];
      expect((list.cell.transform.params as { frameIndex: number }).frameIndex).toBe(3);
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

describe('the trajectory source', () => {
  it('listens to the plugin only while something watches, and reads once per frame', async () => {
    const v = await headlessViewer({ files: [ensemble()] });
    const { model, dispose } = createTrajectoryControls(v.plugin, v.api.as({ kind: 'person' }));
    try {
      const changed = v.plugin.state.data.events.changed as unknown as { observers: unknown[] };
      const idle = changed.observers.length;
      // Reading is not subscribing.
      expect(model.count.value).toBe(3);
      expect(model.index.value).toBe(0);
      expect(changed.observers.length).toBe(idle);

      const watching = model.index.subscribe();
      expect(changed.observers.length).toBeGreaterThan(idle);
      await v.api.trajectory.frame({ index: 2 });
      // Current at once, whether or not the frame has come.
      expect(model.index.value).toBe(2);
      watching.unsubscribe();
      expect(changed.observers.length).toBe(idle);
    } finally {
      dispose();
    }
  });

  it('has commands that are disabled, with the reason, until a structure has several models', async () => {
    const single = await headlessViewer({ files: [cbs()] });
    expect(single.managers.commands.get('trajectory.frame')!.when?.({} as never)).toBe(
      'No structure has several models.',
    );
    const many = await headlessViewer({ files: [ensemble()] });
    expect(many.managers.commands.get('trajectory.frame')!.when?.({} as never)).toBe(true);

    // Playing needs a canvas; with one, it is disabled the same way until there are models.
    const canvas = await createPrismViewerHeadless({ spec: DefaultPrismViewerSpec(), modules: ['canvas', 'files'] });
    try {
      expectOk(await canvas.api.load.files({ files: [cbs()] }));
      const play = canvas.managers.commands.get('trajectory.play')!;
      expect(play.when?.({} as never)).toBe('No structure has several models.');
      expectOk(await canvas.api.load.files({ files: [ensemble()] }, { copy: true }));
      expect(play.when?.({} as never)).toBe(true);
    } finally {
      await canvas.dispose();
    }
    // The fields are labeled ("Mode", "Seconds per pass"), not bare.
    const form = many.api.describe('trajectory.play')!;
    expect([form.params.mode.label, form.params.seconds.label]).toEqual(['Mode', 'Seconds per pass']);
  });
});
