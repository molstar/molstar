/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import { CameraHelperParams } from '@molstar/graphics/canvas3d/helper/camera-helper';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginStateAnimation } from '@molstar/plugin/state/animation/model';
import type { ViewportScreenshotHelper } from '@molstar/plugin/util/viewport-screenshot';
import { type ApiAreaContext, ApiArea, type ApiContext, ApiResult } from '@molstar/prism-viewer/api/area';
import { OptionalParam } from '@molstar/prism-viewer/api/params';
import { holdCanvas } from '@molstar/prism-viewer/sources/canvas-presence';
import { checkValues, messageOf } from '@molstar/prism-viewer/api/run';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension } from '@molstar/prism-viewer/spec';
import { downloadBlob } from '@molstar/prism-viewer/ui/download';
import type { VideoRendererBackend } from './backend.js';
import { VideoCommands, VideoOperations, VideoSection, VideoSkill } from './prism-viewer-stub.js';

declare module '@molstar/prism-viewer/api' {
  // The `video` area: the extension's own (api.md §2, extensions add areas).
  interface ExtensionApis {
    video: VideoApi;
  }
}

// The video export in the Prism viewer (parity.md §6: Export › Animation; lazy): `video.export` renders one of
// the plugin's exportable animations frame by frame on the screenshot helper's image pass and encodes it (Mediabunny,
// WebCodecs), as v5's VideoControls did; `video.cancel` stops it. Exporting changes nothing in the scene: no history
// step, nothing to undo (as `export.image`). The file comes back as bytes; saving it is the page's (`ui/download.ts`).
// The encoder is imported on the first export, so loading the extension for its form costs little.

export type VideoSize = 'viewport' | 'hd' | 'full-hd' | 'ultra-hd' | 'custom';
export type VideoFormat = 'mp4' | 'webm';
export type VideoQuality = 'low' | 'medium' | 'high' | 'very-high';

/** What `video.export` gives: the file's bytes and name, its frame size, length and frame count. */
export interface VideoFile {
  readonly blob: Blob;
  readonly filename: string;
  readonly mimeType: string;
  readonly width: number;
  readonly height: number;
  readonly frames: number;
  readonly seconds: number;
}

const DefaultAnimation = 'built-in.animate-camera-spin';
const MinSide = 128;
const MaxSide = 8192;

/** The animations a video can be made of (v5: `isExportable`). */
export function exportableAnimations(plugin: PluginContext): PluginStateAnimation[] {
  return plugin.managers.animation.animations.filter((a) => a.isExportable);
}

export const VideoExportParams = (ctx: ApiAreaContext) => {
  const options = exportableAnimations(ctx.plugin).map((a) => [a.name, a.display.name] as [string, string]);
  const first = options.find(([n]) => n === DefaultAnimation)?.[0] ?? options[0]?.[0] ?? '';
  return {
    animation: PD.Select<string>(first, options, {
      description: 'The animation to record, by its registered name (spin, rock, snapshots, …)',
    }),
    params: PD.Value<Record<string, unknown>>(
      {},
      { description: "Some of the animation's params (animation.get lists them); the rest are its defaults" },
    ),
    duration: OptionalParam(
      PD.Numeric(
        10,
        { min: 0.1, max: 600, step: 0.1 },
        { description: "Seconds: the video's length (spin and rock play over it); absent: the animation's own" },
      ),
    ),
    fps: PD.Numeric(30, { min: 1, max: 120, step: 1 }, { description: 'Frames per second' }),
    // A preset, or custom with its sides (Prism's api.md §1: a choice with its settings is a union).
    size: PD.MappedStatic(
      'viewport',
      {
        viewport: PD.Group({}, { isFlat: true }),
        hd: PD.Group({}, { isFlat: true }),
        'full-hd': PD.Group({}, { isFlat: true }),
        'ultra-hd': PD.Group({}, { isFlat: true }),
        custom: PD.Group(
          {
            width: PD.Numeric(1920, { min: MinSide, max: MaxSide, step: 1 }, { description: 'Pixels' }),
            height: PD.Numeric(1080, { min: MinSide, max: MaxSide, step: 1 }, { description: 'Pixels' }),
          },
          { isFlat: true },
        ),
      },
      {
        options: [
          ['viewport', 'Viewport'],
          ['hd', 'HD (1280 × 720)'],
          ['full-hd', 'Full HD (1920 × 1080)'],
          ['ultra-hd', 'Ultra HD / 4K (3840 × 2160)'],
          ['custom', 'Custom'],
        ],
        description: "The video's size: a preset, or { name: 'custom', params: { width, height } }",
      },
    ),
    axes: PD.Boolean(false, { description: 'Show the orientation axes' }),
    box: PD.Value<{ x: number; y: number; width: number; height: number }>(
      { x: 0, y: 0, width: 1, height: 1 },
      {
        description:
          'The crop, as fractions of the frame (x, y from the top left, width, height); the whole frame by default',
      },
    ),
    format: PD.Select<VideoFormat>('mp4', [
      ['mp4', 'MP4'],
      ['webm', 'WebM'],
    ]),
    quality: PD.Select<VideoQuality>('high', [
      ['low', 'Low'],
      ['medium', 'Medium'],
      ['high', 'High'],
      ['very-high', 'Very high'],
    ]),
    name: PD.Text('', { description: 'The file name, with or without its extension (empty: from the loaded entries)' }),
  };
};

type ExportValues = PD.ValuesFor<ReturnType<typeof VideoExportParams>>;

const noCanvas = () =>
  ApiResult.fail(
    'unavailable',
    'There is no canvas to render a video from.',
    'Export a video in a viewer with a canvas.',
  );

/** A crop box inside the frame, at least a twentieth of it a side; the whole frame when there is none. */
export function cropOf(b: { x: number; y: number; width: number; height: number } | undefined) {
  if (!b || ![b.x, b.y, b.width, b.height].every(Number.isFinite)) return { x: 0, y: 0, width: 1, height: 1 };
  const x = Math.min(Math.max(b.x, 0), 0.95);
  const y = Math.min(Math.max(b.y, 0), 0.95);
  return { x, y, width: Math.min(Math.max(b.width, 0.05), 1 - x), height: Math.min(Math.max(b.height, 0.05), 1 - y) };
}

/** Sets the screenshot helper to the video's size (the classic UI's too) and returns how to put it back. */
function configure(helper: ViewportScreenshotHelper, v: ExportValues): () => void {
  const before = { values: helper.values, cropParams: helper.cropParams, relativeCrop: helper.relativeCrop };
  const on = PD.getDefaultValues((CameraHelperParams.axes.map('on') as PD.Group<any>).params);
  helper.behaviors.values.next({
    ...before.values,
    resolution: (v.size.name === 'custom' ? v.size : { name: v.size.name, params: {} }) as never,
    // A video has no alpha (MP4): an opaque background, whatever an image export left.
    transparent: false,
    axes: (v.axes ? { name: 'on', params: on } : { name: 'off', params: {} }) as never,
  });
  // The crop given (the Export card's view's), else the whole frame: never found per frame, which would jump.
  helper.behaviors.cropParams.next({ ...before.cropParams, auto: false });
  helper.behaviors.relativeCrop.next(cropOf(v.box));
  return () => {
    helper.behaviors.values.next(before.values);
    helper.behaviors.cropParams.next(before.cropParams);
    helper.behaviors.relativeCrop.next(before.relativeCrop);
  };
}

/** The encoder's values for a format and quality: its default codec for the container, at that quality. */
async function encoderValues(
  backend: VideoRendererBackend,
  config: { width: number; height: number; fps: number },
  format: VideoFormat,
  quality: VideoQuality,
): Promise<Record<string, unknown>> {
  const params = await backend.getParams(config, { container: format });
  const defaults = PD.getDefaultValues(params) as Record<string, any>;
  if (defaults.container !== format) throw new Error(`This browser cannot encode ${format.toUpperCase()} video.`);
  return {
    container: format,
    codec: { name: defaults.codec.name, params: { rate: { name: quality, params: {} } } },
  };
}

function fileNameOf(helper: ViewportScreenshotHelper, anim: PluginStateAnimation, v: ExportValues, ext: string) {
  const entered = v.name.trim().replace(/[\\/:*?"<>|]+/g, '-');
  if (entered) return new RegExp(`\\.${ext}$`, 'i').test(entered) ? entered : `${entered}.${ext}`;
  const kind = anim.display.name
    .toLowerCase()
    .replace(/\s/g, '-')
    .replace(/[^a-z0-9_-]/g, '');
  return `${helper.getFilename('')}_${kind}.${ext}`;
}

export class VideoApi extends ApiArea {
  readonly id = 'video';

  /** Stops the export in progress, when there is one. */
  #cancel: (() => void) | undefined;

  readonly export = this.operation({
    ...VideoOperations.export,
    params: VideoExportParams,
    // A file for the person; the scene is as it was.
    history: false,
    run: (v, ctx) => this.#export(v, ctx),
  });

  readonly cancel = this.operation({
    ...VideoOperations.cancel,
    params: {},
    history: false,
    run: async () => {
      if (!this.#cancel) return ApiResult.already('No video is being exported.');
      this.#cancel();
      return ApiResult.ok('Cancelling the video export.');
    },
  });

  async #export(v: ExportValues, ctx: ApiContext): Promise<ApiResult<VideoFile>> {
    if (this.#cancel) {
      return ApiResult.fail('blocked', 'A video is being exported.', 'Wait for it, or cancel it (video.cancel).');
    }
    const plugin = ctx.plugin;
    const all = exportableAnimations(plugin);
    const anim = all.find((a) => a.name === v.animation);
    if (!anim) {
      return ApiResult.fail(
        'invalid',
        `There is no animation '${v.animation}' to export.`,
        `Animations: ${all.map((a) => a.name).join(', ') || 'none'}.`,
      );
    }
    const can = anim.canApply?.(plugin);
    if (can && !can.canApply) {
      return ApiResult.fail('blocked', `${anim.display.name} cannot play now: ${can.reason ?? 'nothing to animate'}.`);
    }
    const checked = checkValues(anim.params(plugin) as PD.Params, v.params ?? {});
    if (ApiResult.isFailure(checked)) return checked;
    // A duration sets the animation's own (spin, rock) and the video's length (any animation with no end of its own).
    const duration = (v as { duration?: number }).duration;
    const customMs = duration !== undefined ? Math.round(duration * 1000) : undefined;
    if (customMs !== undefined && 'durationInMs' in checked.values) checked.values.durationInMs = customMs;
    const lengthMs = customMs ?? fixedDuration(plugin, anim, checked.values);
    if (!lengthMs) {
      return ApiResult.fail(
        'invalid',
        `${anim.display.name} has no length of its own.`,
        'Give a duration in seconds (snapshots: add some first).',
      );
    }
    const sides = v.size.name === 'custom' ? (v.size.params as { width?: number; height?: number }) : undefined;
    if (
      sides &&
      !(sides.width! >= MinSide && sides.width! <= MaxSide && sides.height! >= MinSide && sides.height! <= MaxSide)
    ) {
      return ApiResult.fail('invalid', `A custom size needs a width and a height of ${MinSide} to ${MaxSide} pixels.`);
    }
    // The values are checked first: they do not need the canvas.
    const helper = plugin.helpers.viewportScreenshot;
    if (!plugin.canvas3d || !helper) return noCanvas();
    if (ctx.signal.aborted) return ApiResult.fail('aborted', 'Stopped before it started.', 'Nothing was rendered.');

    // Rendered on purpose: the canvas is held for the whole export (its context restored if it was released; its loop
    // left to the encoder, not paused or started by the canvas' presence: Prism's README decision 29).
    const unhold = await holdCanvas(plugin);
    const restore = configure(helper, v);
    let aborted = false;
    const stop = () => {
      aborted = true;
      if (task) plugin.managers.task.requestAbort(task.id, 'Cancelled');
    };
    let task: Task<{ data: Uint8Array<ArrayBuffer>; mimeType: string; extension: string }> | undefined;
    this.#cancel = stop;
    ctx.signal.addEventListener('abort', stop);
    try {
      const [{ encodeVideoAnimation }, { MediabunnyBackend }] = await Promise.all([
        import('./encoder.js'),
        import('./backends/mediabunny.js'),
      ]);
      const { width, height, viewport } = helper.getSizeAndViewport();
      // Even sides: H.264 and most players need them.
      const frame = {
        ...viewport,
        width: viewport.width - (viewport.width % 2),
        height: viewport.height - (viewport.height % 2),
      };
      let backendParams: Record<string, unknown>;
      try {
        backendParams = await encoderValues(MediabunnyBackend, { ...frame, fps: v.fps }, v.format, v.quality);
      } catch (e) {
        return ApiResult.fail('unavailable', messageOf(e), 'Try the other format, or a browser with WebCodecs.');
      }
      if (aborted) return ApiResult.fail('aborted', 'The video export was cancelled.', 'Nothing was saved.');
      const t = Task.create('Export a video', (rt) =>
        encodeVideoAnimation(plugin, rt, {
          animation: { definition: anim, params: checked.values, ...(customMs ? { customDurationMs: customMs } : {}) },
          width,
          height,
          viewport: frame,
          fps: v.fps,
          backend: MediabunnyBackend,
          backendParams,
          pass: helper.imagePass,
        }),
      );
      task = t;
      const movie = await plugin.runTask(t, { useOverlay: true });
      const filename = fileNameOf(helper, anim, v, movie.extension);
      const seconds = lengthMs / 1000;
      const frames = Math.ceil(seconds * v.fps);
      return ApiResult.ok(
        `Rendered ${filename} (${frame.width} × ${frame.height}, ${seconds.toFixed(1)} s).`,
        {},
        {
          value: {
            blob: new Blob([movie.data], { type: movie.mimeType }),
            filename,
            mimeType: movie.mimeType,
            width: frame.width,
            height: frame.height,
            frames,
            seconds,
          },
        },
      );
    } catch (e) {
      if (aborted || Task.isAbort(e)) {
        return ApiResult.fail('aborted', 'The video export was cancelled.', 'Nothing was saved.');
      }
      return ApiResult.fail('failed', `Could not export the video: ${messageOf(e)}`);
    } finally {
      ctx.signal.removeEventListener('abort', stop);
      this.#cancel = undefined;
      restore();
      unhold();
    }
  }
}

/** The animation's length in milliseconds when it has a fixed one (0: none, or it never ends). */
function fixedDuration(plugin: PluginContext, anim: PluginStateAnimation, params: Record<string, unknown>): number {
  const d = anim.getDuration?.(params, plugin);
  return d?.kind === 'fixed' ? d.durationMs : 0;
}

export const VideoExportExtension: ViewerExtension = {
  id: 'video-export',
  title: 'Video export',
  api: { video: (ctx) => new VideoApi(ctx) },
  commands: (api) => {
    const x = api.video;
    if (!x) return [];
    return [
      {
        ...operationCommand(api, x.export, {
          category: 'Export',
          icon: 'clapperboard',
          title: 'Export a video…',
          alias: ['video', 'movie', 'mp4', 'webm', 'record animation', 'export animation'],
          action: 'Render',
        }),
        // The operation gives the file's bytes; the command saves it through the page.
        async run(args) {
          const r = await api.video!.export(args as never);
          if (!r.ok) throw new Error(r.hint ? `${r.message} ${r.hint}` : r.message);
          if (r.value) downloadBlob(r.value.blob, r.value.filename);
          return { message: `Saved ${r.value?.filename}.`, ops: `video.export(${JSON.stringify(args)})` };
        },
      },
      operationCommand(api, x.cancel, { category: 'Export', icon: 'x', alias: ['stop video', 'cancel render'] }),
      {
        id: VideoCommands.open,
        title: 'Show the Animation export',
        category: 'Export',
        icon: 'clapperboard',
        alias: ['video export', 'animation export'],
        // Opens the card at the section: a view, not a change.
        async run() {
          // Animation is a choice of the Image section's (`export.picture`).
          const r = await api.view.open({
            card: 'core.overlay.export',
            section: `export-image.section.image/${VideoSection}`,
          });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Export card.', ops: 'view.open({ card: "export" })' };
        },
      },
    ];
  },
  skills: [VideoSkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).videoExportUI(ctx),
};
