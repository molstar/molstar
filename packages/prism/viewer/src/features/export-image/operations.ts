/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CameraHelperParams } from '@molstar/graphics/canvas3d/helper/camera-helper';
import type { ViewportScreenshotHelper } from '@molstar/plugin/util/viewport-screenshot';
import { type ApiContext, ApiArea, ApiResult } from '../../api/area.js';

// The `export` area, image part (operations.md §2, the `export` row): an image of the canvas by Mol*'s
// `ViewportScreenshotHelper`. Exporting changes nothing in the scene, so none of these is a history step
// (`history: false`) and none is undoable. `history: false` does NOT serialize them: `api/run.ts` runs it like a read,
// directly. They all configure the one shared helper (also the classic UI's), render and put it back, so `ExportApi`
// runs them one at a time through its own queue. They give the image as bytes (a Blob and a file name); saving it is
// the page's (`ui/download.ts`), as for `session.save`. The helper exists only where there is a canvas in a page;
// elsewhere they fail `unavailable`.

/** The crop as fractions of the image (0-1), like Mol*'s `Viewport.relativeCrop`. */
export interface CropBox {
  x: number;
  y: number;
  width: number;
  height: number;
}

export type ImageSize = 'viewport' | 'hd' | 'full-hd' | 'ultra-hd' | '8k-ultra-hd' | 'custom';
export type ImageFormat = 'png' | 'jpeg' | 'webp';

/** What the image is, shared by `export.image`, `export.image-copy` and the preview. */
const ShotParams = {
  size: PD.Select<ImageSize>('viewport', [
    ['viewport', 'Viewport'],
    ['hd', 'HD (1280 x 720)'],
    ['full-hd', 'Full HD (1920 x 1080)'],
    ['ultra-hd', 'Ultra HD / 4K (3840 x 2160)'],
    ['8k-ultra-hd', '8K Ultra HD (7680 x 4320)'],
    ['custom', 'Custom'],
  ]),
  width: PD.Numeric(
    0,
    { min: 0, max: 8192, step: 1 },
    { description: 'Pixels, with height: a custom size (0: use size)' },
  ),
  height: PD.Numeric(
    0,
    { min: 0, max: 8192, step: 1 },
    { description: 'Pixels, with width: a custom size (0: use size)' },
  ),
  transparent: PD.Boolean(false, { description: 'Transparent background (PNG and WebP)' }),
  crop: PD.Select<'auto' | 'none'>(
    'auto',
    [
      ['auto', 'Auto: to the visible content'],
      ['none', 'None: the whole frame'],
    ],
    { description: 'Crop to the visible content, or keep the whole frame' },
  ),
  // `undefined` and not null: a command form shows a null value as a file picker.
  box: PD.Value<CropBox | undefined>(undefined, {
    isHidden: true,
    description: 'A crop as fractions of the frame ({ x, y, width, height }, 0-1); replaces crop',
  }),
  format: PD.Select<ImageFormat>('png', [
    ['png', 'PNG'],
    ['jpeg', 'JPEG'],
    ['webp', 'WebP'],
  ]),
  quality: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }, { description: 'JPEG and WebP quality' }),
  axes: PD.Boolean(false, { description: 'Show the orientation axes' }),
};

export const ImageParams = {
  ...ShotParams,
  name: PD.Text('', { description: 'The file name, with or without its extension (empty: from the loaded entries)' }),
};
export const CopyParams = ShotParams;
export const PreviewParams = {
  ...ShotParams,
  max: PD.Numeric(320, { min: 64, max: 1024, step: 1 }, { description: 'The preview its longest side in pixels' }),
};

export type ShotValues = PD.ValuesFor<typeof ShotParams>;

/** What `export.image` gives: the image's bytes and file name, and its size (after the crop). */
export interface ImageResult {
  blob: Blob;
  filename: string;
  width: number;
  height: number;
}

/** What `export.image-copy` gives: the size of the image copied (after the crop). */
export interface CopiedImage {
  width: number;
  height: number;
}

/** What `export.preview` gives: a small image of the whole frame, its full size, and the crop on it. */
export interface PreviewResult {
  dataUri: string;
  /** The frame in pixels, before the crop. */
  width: number;
  height: number;
  crop: CropBox;
}

const Full: CropBox = { x: 0, y: 0, width: 1, height: 1 };
const MaxSide = 8192;
const MinSide = 128;

function helperOf(ctx: ApiContext): ViewportScreenshotHelper | undefined {
  return ctx.plugin.helpers.viewportScreenshot;
}

const noCanvas = () =>
  ApiResult.fail(
    'unavailable',
    'There is no canvas to render an image from.',
    'Export an image in a viewer with a canvas.',
  );

function validBox(b: CropBox): boolean {
  return (
    [b.x, b.y, b.width, b.height].every((n) => Number.isFinite(n)) &&
    b.x >= 0 &&
    b.y >= 0 &&
    b.width > 0 &&
    b.height > 0 &&
    b.x + b.width <= 1 + 1e-6 &&
    b.y + b.height <= 1 + 1e-6
  );
}

/** Checks the values a render needs; a message when one is wrong. */
function check(v: ShotValues): string | undefined {
  const custom = v.size === 'custom' || v.width > 0 || v.height > 0;
  if (custom) {
    if (!(v.width >= MinSide && v.width <= MaxSide && v.height >= MinSide && v.height <= MaxSide)) {
      return `A custom size needs a width and a height of ${MinSide} to ${MaxSide} pixels.`;
    }
  }
  if (v.box && !validBox(v.box))
    return 'The crop box is not inside the frame: x, y, width and height are fractions of 0 to 1.';
  return undefined;
}

/**
 * Sets the helper to render `v` and returns how to put it back: the helper is the classic UI's too, so what an export
 * changed in it is restored, whether the render worked or not.
 */
function configure(helper: ViewportScreenshotHelper, v: ShotValues): () => void {
  const before = { values: helper.values, cropParams: helper.cropParams, relativeCrop: helper.relativeCrop };
  const custom = v.size === 'custom' || v.width > 0 || v.height > 0;
  const on = PD.getDefaultValues((CameraHelperParams.axes.map('on') as PD.Group<any>).params);
  helper.behaviors.values.next({
    ...before.values,
    resolution: custom
      ? { name: 'custom', params: { width: v.width, height: v.height } }
      : ({ name: v.size, params: {} } as never),
    format: v.format === 'png' ? { name: 'png', params: {} } : { name: v.format, params: { quality: v.quality } },
    transparent: v.transparent,
    axes: (v.axes ? { name: 'on', params: on } : { name: 'off', params: {} }) as never,
  });
  const auto = v.crop === 'auto' && !v.box;
  helper.behaviors.cropParams.next({ ...before.cropParams, auto });
  helper.behaviors.relativeCrop.next(v.box ? { ...v.box } : { ...Full });
  return () => {
    helper.behaviors.values.next(before.values);
    helper.behaviors.cropParams.next(before.cropParams);
    helper.behaviors.relativeCrop.next(before.relativeCrop);
  };
}

/** The crop is found on a small render (the helper does it in `getPreview`). */
async function settleCrop(ctx: ApiContext, helper: ViewportScreenshotHelper, max = 512) {
  if (!helper.cropParams.auto) return;
  await ctx.plugin.runTask(Task.create('Find the crop', (rt) => helper.getPreview(rt, max)));
}

const dims = (helper: ViewportScreenshotHelper) => {
  const { viewport } = helper.getSizeAndViewport();
  return { width: viewport.width, height: viewport.height };
};

const message = (e: unknown) => (e instanceof Error ? e.message : String(e));

/** A `data:` URI as a Blob (the helper gives the image as one; there is no blob getter). */
function dataUriToBlob(uri: string): Blob {
  const comma = uri.indexOf(',');
  const head = uri.slice(5, comma);
  const type = head.split(';')[0] || 'application/octet-stream';
  const body = uri.slice(comma + 1);
  if (!head.endsWith(';base64')) return new Blob([decodeURIComponent(body)], { type });
  const binary = atob(body);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
  return new Blob([bytes], { type });
}

const Extensions: Record<ImageFormat, string> = { png: '.png', jpeg: '.jpg', webp: '.webp' };

/** The file name: the person's, with the format's extension unless it has an image one; else the helper's. */
function fileNameOf(helper: ViewportScreenshotHelper, v: ShotValues & { name?: string }): string {
  const entered = (v.name ?? '').trim().replace(/[\\/:*?"<>|]+/g, '-');
  if (!entered) return helper.getFilename();
  return /\.(png|jpe?g|webp)$/i.test(entered) ? entered : entered + Extensions[v.format];
}

/** Runs one job after the ones before it, whatever they ended in. */
class Queue {
  #tail: Promise<unknown> = Promise.resolve();
  run<T>(job: () => Promise<T>): Promise<T> {
    const next = this.#tail.then(job, job);
    this.#tail = next.catch(() => undefined);
    return next;
  }
}

const stopped = () => ApiResult.fail('aborted', 'Stopped before it started.', 'Nothing was rendered.');

export class ExportApi extends ApiArea {
  readonly id = 'export';

  // Configure, render and restore of one operation must not interleave with another's: one queue for the area.
  readonly #queue = new Queue();

  readonly image = this.operation({
    id: 'export.image',
    title: 'Export an image',
    params: ImageParams,
    modules: ['canvas'],
    history: false,
    undoable: false,
    run: (v, ctx): Promise<ApiResult<ImageResult>> =>
      this.#queue.run(async () => {
        if (ctx.signal.aborted) return stopped();
        const helper = helperOf(ctx);
        if (!helper) return noCanvas();
        const bad = check(v);
        if (bad) return ApiResult.fail('invalid', bad);
        const restore = configure(helper, v);
        try {
          await settleCrop(ctx, helper);
          const size = dims(helper);
          const filename = fileNameOf(helper, v);
          const blob = dataUriToBlob(await helper.getImageDataUri());
          return ApiResult.ok(
            `Rendered ${filename} (${size.width} × ${size.height}).`,
            {},
            { value: { ...size, blob, filename } },
          );
        } catch (e) {
          return ApiResult.fail('failed', `Could not render the image: ${message(e)}`);
        } finally {
          restore();
        }
      }),
  });

  readonly imageCopy = this.operation({
    id: 'export.image-copy',
    title: 'Copy an image',
    params: CopyParams,
    modules: ['canvas', 'page'],
    history: false,
    undoable: false,
    run: (v, ctx): Promise<ApiResult<CopiedImage>> =>
      this.#queue.run(async () => {
        if (ctx.signal.aborted) return stopped();
        const helper = helperOf(ctx);
        if (!helper) return noCanvas();
        const bad = check(v);
        if (bad) return ApiResult.fail('invalid', bad);
        const restore = configure(helper, v);
        try {
          await settleCrop(ctx, helper);
          const size = dims(helper);
          const done = helper.copyToClipboard();
          if (!done)
            return ApiResult.fail('unavailable', 'This browser cannot copy images.', 'Download the image instead.');
          await done;
          return ApiResult.ok(
            `Copied the image (${size.width} × ${size.height}) to the clipboard.`,
            {},
            { value: size },
          );
        } catch (e) {
          return ApiResult.fail('failed', `Could not copy the image: ${message(e)}`);
        } finally {
          restore();
        }
      }),
  });

  readonly preview = this.operation({
    id: 'export.preview',
    title: 'Preview an image',
    params: PreviewParams,
    modules: ['canvas'],
    history: false,
    undoable: false,
    surface: { kind: 'ui-only', reason: "The Image section's preview; tools export the image itself" },
    run: (v, ctx): Promise<ApiResult<PreviewResult>> =>
      this.#queue.run(async () => {
        if (ctx.signal.aborted) return stopped();
        const helper = helperOf(ctx);
        if (!helper) return noCanvas();
        const bad = check(v);
        if (bad) return ApiResult.fail('invalid', bad);
        const restore = configure(helper, v);
        try {
          const shot = await ctx.plugin.runTask(Task.create('Preview', (rt) => helper.getPreview(rt, v.max)));
          if (!shot) return ApiResult.fail('failed', 'The canvas has no size to preview.');
          const { width, height } = helper.getSizeAndViewport();
          const value: PreviewResult = {
            dataUri: shot.canvas.toDataURL('image/png'),
            width,
            height,
            crop: { ...helper.relativeCrop },
          };
          return ApiResult.ok('Previewed the image.', {}, { value });
        } catch (e) {
          return ApiResult.fail('failed', `Could not preview the image: ${message(e)}`);
        } finally {
          restore();
        }
      }),
  });
}
