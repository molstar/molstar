/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { ApiModuleKey } from '@molstar/prism-viewer/api/area';
import { DefaultPrismViewerSpec } from '@molstar/prism-viewer/default-spec';
import { createPrismViewerHeadless } from '@molstar/prism-viewer/headless';
import type { LazyViewerExtension, ViewerExtension } from '@molstar/prism-viewer/spec';
import type { PrismViewer } from '@molstar/prism-viewer/viewer';
import { afterEach, describe, expect, it } from 'vitest';
import { VideoExportExtension } from '../prism-viewer.js';
import { LazyVideoExport, VideoCommands, VideoOperations } from '../prism-viewer-stub.js';

// The video export in the Prism viewer (plan 5.2) on a headless viewer: no canvas, so no frame is rendered (WebGL and
// WebCodecs are the browser's); what is tested is the stub, the load, the checks before rendering and the refusals.

const viewers: PrismViewer[] = [];
afterEach(async () => {
  for (const v of viewers.splice(0)) await v.dispose();
});

async function viewer(ext: ViewerExtension | LazyViewerExtension, modules?: readonly ApiModuleKey[]) {
  const base = DefaultPrismViewerSpec();
  const v = await createPrismViewerHeadless({
    spec: { ...base, extensions: [...(base.extensions ?? []), ext] },
    ...(modules ? { modules } : {}),
  });
  viewers.push(v);
  return v;
}

describe('the video-export extension in the Prism viewer', () => {
  it('is a stub until first used: its operations and commands are there before it loads', async () => {
    const v = await viewer(LazyVideoExport);
    expect(v.extensions.status('video-export')?.value).toBe('stub');
    for (const op of Object.values(VideoOperations)) expect(v.api.describe(op.id)?.definition.title).toBe(op.title);
    for (const id of Object.values(VideoCommands)) expect(v.managers.commands.get(id)).toBeDefined();
    // Without a canvas the export cannot run, and says so before loading anything.
    const r = await v.api.run('video.export', {});
    expect(r.ok || r.code).toBe('unavailable');
    expect(v.extensions.status('video-export')?.value).toBe('stub');
    // Cancelling (nothing runs) loads it, and the real operation answers.
    const c = await v.api.run('video.cancel', {});
    expect(c.ok && c.unchanged).toBe(true);
    expect(v.extensions.status('video-export')?.value).toBe('loaded');
  });

  it('loads as its stub declares it', async () => {
    const ext = await LazyVideoExport.load();
    expect(ext).toBe(VideoExportExtension);
    const v = await viewer(LazyVideoExport);
    expect(await v.extensions.load('video-export')).toBe(true);
    for (const op of LazyVideoExport.stub.operations ?? []) {
      const def = v.api.describe(op.id)?.definition;
      expect(def?.title).toBe(op.title);
      expect(def?.undoable).toBe(op.undoable);
      expect(def?.modules).toEqual(op.modules);
    }
    for (const c of LazyVideoExport.stub.commands ?? []) expect(v.managers.commands.get(c.id)?.title).toBe(c.title);
    expect(ext.skills?.map((s) => s.id)).toEqual(LazyVideoExport.stub.skills?.map((s) => s.id));
    // Still no canvas: the loaded operation refuses too, saying why.
    const r = await v.api.video!.export({});
    expect(r.ok ? 'ok' : r.code).toBe('unavailable');
    expect(r.ok ? '' : r.message).toMatch(/canvas/);
  });

  it('checks the animation and the size before it needs the canvas', async () => {
    // A host that says it has a canvas, on a plugin that has none: the checks come first, then the refusal.
    const v = await viewer(VideoExportExtension, ['canvas', 'files', 'network']);
    const video = v.api.video!;
    const options = (v.api.describe('video.export')?.params.animation as PD.Select<string>).options;
    expect(options.map(([n]) => n)).toEqual(
      expect.arrayContaining(['built-in.animate-camera-spin', 'built-in.animate-camera-rock']),
    );
    expect((await video.export({ animation: 'no-such-animation' })).ok).toBe(false);
    // Snapshots need some.
    const snapshots = await video.export({ animation: 'built-in.animate-state-snapshots' });
    expect(snapshots.ok ? 'ok' : snapshots.code).toBe('blocked');
    const params = await video.export({ params: { noSuchParam: 1 } });
    expect(params.ok ? 'ok' : params.code).toBe('invalid');
    const size = await video.export({ size: 'custom', width: 64, height: 64 });
    expect(size.ok ? 'ok' : size.code).toBe('invalid');
    const spin = await video.export({ animation: 'built-in.animate-camera-spin', duration: 2, fps: 24 });
    expect(spin.ok ? 'ok' : spin.code).toBe('unavailable');
    // Not a step: an export changes nothing the history keeps.
    expect(v.managers.history.steps.value.length).toBe(0);
  });

  it('brings its commands and its skill', async () => {
    const v = await viewer(VideoExportExtension);
    for (const id of Object.values(VideoCommands)) expect(v.managers.commands.get(id)).toBeDefined();
    expect(VideoExportExtension.skills?.[0].tools.map((t) => t.operation)).toEqual(['video.export', 'video.cancel']);
    const r = await v.api.video!.cancel();
    expect(r.ok && r.unchanged).toBe(true);
  });
});
