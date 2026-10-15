/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { cbs, expectOk, headlessViewer, withExtensions } from '../../../_test/harness.js';
import { DefaultPrismViewerSpec } from '../../../default-spec.js';
import { createPrismViewerHeadless } from '../../../headless.js';
import type { PrismViewer } from '../../../viewer.js';
import { ExportImageFeature } from '../index.js';
import { createImageExportModel, imageValues } from '../model.js';

// The `export-image` feature (operations.md §2, `export`): its operations on a stand-in `ViewportScreenshotHelper`
// (the real one needs a canvas, WebGL and a page: what it renders is not tested here).

const Crop = { x: 0.1, y: 0.2, width: 0.5, height: 0.4 };

/** A behavior like the plugin's: a value, and `next`. */
const behavior = <T>(value: T) => ({
  value,
  next(v: T) {
    this.value = v;
  },
});

/** The helper's surface the operations use, recording what it was asked to render. */
function standIn() {
  const calls: { op: string; values: any; crop: any; cropParams: any; name?: string }[] = [];
  const behaviors = {
    values: behavior<any>({
      transparent: false,
      format: { name: 'png', params: {} },
      axes: { name: 'off', params: {} },
      resolution: { name: 'viewport', params: {} },
      illumination: { extraIterations: 1, targetIterationTimeMs: 300 },
    }),
    cropParams: behavior({ auto: true, relativePadding: 0.1 }),
    relativeCrop: behavior({ x: 0, y: 0, width: 1, height: 1 }),
  };
  const sizes: Record<string, [number, number]> = {
    viewport: [1000, 500],
    hd: [1280, 720],
    'full-hd': [1920, 1080],
  };
  const frame = () => {
    const r = behaviors.values.value.resolution;
    const [width, height] = r.name === 'custom' ? [r.params.width, r.params.height] : sizes[r.name];
    return { width, height };
  };
  const note = (op: string, name?: string) =>
    calls.push({
      op,
      values: structuredClone(behaviors.values.value),
      crop: { ...behaviors.relativeCrop.value },
      cropParams: { ...behaviors.cropParams.value },
      name,
    });
  const helper = {
    behaviors,
    get values() {
      return behaviors.values.value;
    },
    get cropParams() {
      return behaviors.cropParams.value;
    },
    get relativeCrop() {
      return behaviors.relativeCrop.value;
    },
    getSizeAndViewport() {
      const { width, height } = frame();
      const c = behaviors.relativeCrop.value;
      return {
        width,
        height,
        viewport: { x: 0, y: 0, width: Math.ceil(c.width * width), height: Math.ceil(c.height * height) },
      };
    },
    // As the real one: the preview finds the crop when it is automatic.
    async getPreview() {
      note('preview');
      if (behaviors.cropParams.value.auto) behaviors.relativeCrop.next({ ...Crop });
      const { width, height } = frame();
      return { canvas: { toDataURL: () => 'data:image/png;base64,UFJFVklFVw==' }, width, height };
    },
    async getImageDataUri() {
      note('image');
      return 'data:image/png;base64,SU1BR0U=';
    },
    async download(name?: string) {
      note('download', name);
    },
    copyToClipboard() {
      note('copy');
      return Promise.resolve();
    },
    getFilename: () => '1CBS.png',
  };
  return { helper, calls };
}

const withFeature = () => {
  const spec = DefaultPrismViewerSpec();
  return { ...spec, extensions: withExtensions(spec.extensions ?? [], [ExportImageFeature]) };
};

/** A viewer whose host has a canvas and a page, with the stand-in helper. */
async function withHelper(modules: ('canvas' | 'files' | 'page')[] = ['canvas', 'files', 'page']) {
  const v = await createPrismViewerHeadless({ spec: withFeature(), modules });
  const s = standIn();
  (v.plugin.helpers as { viewportScreenshot: unknown }).viewportScreenshot = s.helper;
  return { v, ...s };
}

let open: PrismViewer | undefined;
const track = async (p: ReturnType<typeof withHelper>) => {
  const r = await p;
  open = r.v;
  return r;
};
afterEach(async () => {
  await open?.dispose();
  open = undefined;
});

describe('export image', () => {
  it('adds the export area with its operations and commands, none a history step', async () => {
    const v = await headlessViewer({ extensions: [ExportImageFeature] });
    expect(v.api.export).toBeDefined();
    for (const id of ['export.image', 'export.image-copy', 'export.preview']) expect(v.api.describe(id)).toBeDefined();
    expect(v.api.describe('export.image')?.definition).toMatchObject({
      history: false,
      undoable: false,
      modules: ['canvas'],
    });
    expect(v.api.describe('export.image-copy')?.definition.modules).toEqual(['canvas', 'page']);
    expect(v.api.describe('export.preview')?.definition.surface?.kind).toBe('ui-only');
    for (const id of ['export.image', 'export.image-copy', 'export.show']) {
      expect(v.managers.commands.get(id), id).toBeDefined();
    }
  });

  it('says it is unavailable without a canvas', async () => {
    const v = await headlessViewer({ extensions: [ExportImageFeature], files: [cbs()] });
    expect(await v.api.export!.image()).toMatchObject({ ok: false, code: 'unavailable' });
    expect(await v.api.export!.imageCopy()).toMatchObject({ ok: false, code: 'unavailable' });
    expect(await v.api.export!.preview()).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('says it is unavailable with a canvas host but no helper (a headless plugin)', async () => {
    const v = await createPrismViewerHeadless({
      spec: withFeature(),
      modules: ['canvas', 'files'],
    });
    open = v;
    expect(await v.api.export!.image()).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('renders a data URI at the size asked, cropped, and leaves the helper as it was', async () => {
    const { v, calls, helper } = await track(withHelper());
    const steps = v.managers.history.steps.value.length;
    const r = await v.api.export!.image({ size: 'hd', transparent: true, axes: true });
    expectOk(r);
    // The automatic crop (found on a preview) is what is cut: 0.5 x 0.4 of 1280 x 720.
    expect(r.value).toEqual({ dataUri: 'data:image/png;base64,SU1BR0U=', width: 640, height: 288 });
    expect(calls.map((c) => c.op)).toEqual(['preview', 'image']);
    expect(calls[1].values).toMatchObject({
      resolution: { name: 'hd' },
      transparent: true,
      axes: { name: 'on' },
      format: { name: 'png' },
    });
    expect(calls[1].crop).toEqual(Crop);
    // Nothing of the export stays in the helper (the classic UI's), or in the history.
    expect(helper.values).toMatchObject({
      resolution: { name: 'viewport' },
      transparent: false,
      axes: { name: 'off' },
    });
    expect(helper.cropParams.auto).toBe(true);
    expect(helper.relativeCrop).toEqual({ x: 0, y: 0, width: 1, height: 1 });
    expect(v.managers.history.steps.value.length).toBe(steps);
  });

  it('takes a custom size, no crop, a box, and JPEG with a quality', async () => {
    const { v, calls } = await track(withHelper());
    const none = await v.api.export!.image({ width: 800, height: 600, crop: 'none' });
    expectOk(none);
    expect(none.value).toMatchObject({ width: 800, height: 600 });
    expect(calls.map((c) => c.op)).toEqual(['image']);
    expect(calls[0].values.resolution).toEqual({ name: 'custom', params: { width: 800, height: 600 } });

    const box = await v.api.export!.image({
      size: 'full-hd',
      format: 'jpeg',
      quality: 0.5,
      box: { x: 0, y: 0, width: 0.5, height: 0.5 },
    });
    expectOk(box);
    expect(box.value).toMatchObject({ width: 960, height: 540 });
    expect(calls[1]).toMatchObject({ op: 'image', crop: { width: 0.5 }, cropParams: { auto: false } });
    expect(calls[1].values.format).toEqual({ name: 'jpeg', params: { quality: 0.5 } });
  });

  it('fails invalid for a size or a box that cannot be', async () => {
    const { v, calls } = await track(withHelper());
    expect(await v.api.export!.image({ size: 'custom' })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.export!.image({ width: 100, height: 600 })).toMatchObject({ ok: false, code: 'invalid' });
    expect(await v.api.export!.image({ box: { x: 0.8, y: 0, width: 0.5, height: 1 } })).toMatchObject({
      ok: false,
      code: 'invalid',
    });
    expect(calls).toEqual([]);
  });

  it('downloads under a name, or the helper’s; and needs the page for it', async () => {
    const { v, calls } = await track(withHelper());
    const r = await v.api.export!.image({ download: true, name: 'mine.png', crop: 'none' });
    expectOk(r);
    expect(r.value).toMatchObject({ filename: 'mine.png' });
    expect(r.value?.dataUri).toBeUndefined();
    expect(calls[0]).toMatchObject({ op: 'download', name: 'mine.png' });
    expectOk(await v.api.export!.image({ download: true, crop: 'none' }));
    expect(calls[1]).toMatchObject({ op: 'download', name: '1CBS.png' });

    const noPage = await track(withHelper(['canvas', 'files']));
    expect(await noPage.v.api.export!.image({ download: true })).toMatchObject({ ok: false, code: 'unavailable' });
    expect(noPage.calls).toEqual([]);
  });

  it('copies to the clipboard, only where there is a page; fails cleanly when the browser cannot', async () => {
    const { v, calls, helper } = await track(withHelper());
    expectOk(await v.api.export!.imageCopy({ size: 'hd', crop: 'none' }));
    expect(calls[0]).toMatchObject({ op: 'copy', values: { resolution: { name: 'hd' } } });
    (helper as { copyToClipboard: unknown }).copyToClipboard = () => undefined;
    expect(await v.api.export!.imageCopy({ crop: 'none' })).toMatchObject({ ok: false, code: 'unavailable' });
    (helper as { copyToClipboard: unknown }).copyToClipboard = () => Promise.reject(new Error('denied'));
    expect(await v.api.export!.imageCopy({ crop: 'none' })).toMatchObject({ ok: false, code: 'failed' });

    const noPage = await track(withHelper(['canvas', 'files']));
    expect(await noPage.v.api.export!.imageCopy()).toMatchObject({ ok: false, code: 'unavailable' });
  });

  it('previews the whole frame with the crop on it', async () => {
    const { v } = await track(withHelper());
    const r = await v.api.export!.preview({ size: 'hd' });
    expectOk(r);
    expect(r.value).toEqual({ dataUri: 'data:image/png;base64,UFJFVklFVw==', width: 1280, height: 720, crop: Crop });
    const whole = await v.api.export!.preview({ size: 'hd', crop: 'none' });
    expectOk(whole);
    expect(whole.value?.crop).toEqual({ x: 0, y: 0, width: 1, height: 1 });
  });

  it('shows the Export card at the Image section', async () => {
    const v = await headlessViewer({ extensions: [ExportImageFeature] });
    await v.managers.commands.get('export.show')!.run({}, { empty: true }, new AbortController().signal);
    expect(v.managers.navigation.state.value.card).toEqual({ id: 'export', section: 'export-image.section.image' });
  });

  it('keeps the model in step with the form: its size, aspect and crop; copy and download run the operations', async () => {
    const { v, calls } = await track(withHelper());
    const posted: string[] = [];
    const model = createImageExportModel({
      api: v.api,
      sources: v.sources,
      notify: (t) => posted.push(t),
    });
    try {
      model.form.set('resolution', { name: 'hd', params: {} } as never);
      expect(model.aspect.value).toBeCloseTo(16 / 9);
      expect(model.size.value).toBe('1280 × 720 px · PNG');
      model.setCrop({ x: 0, y: 0, width: 0.5, height: 0.5 });
      expect(model.autoCrop.value).toBe(false);
      expect(model.size.value).toBe('640 × 360 px · PNG');
      model.copy();
      await vi.waitFor(() => expect(posted).toEqual(['Copied the image.']));
      expect(calls.find((c) => c.op === 'copy')).toMatchObject({ crop: { width: 0.5 } });
      model.download();
      await vi.waitFor(() => expect(posted.length).toBe(2));
      expect(posted[1]).toMatch(/^Saved 1CBS\.png/);
      expect(imageValues(model.form.values.value)).toMatchObject({ size: 'hd', width: 0, height: 0 });
    } finally {
      model.dispose();
    }
  });
});
