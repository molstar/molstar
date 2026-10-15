/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { MessageKind } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { observableState } from '@molstar/prism-ui/reactive/observable-state';
import { combineLatest, debounceTime, type Subscription, tap } from 'rxjs';
import type { ViewerApi } from '../../api/area.js';
import type { ImageExportModel } from '../../model/image-export.js';
import type { ViewerSources } from '../../sources/viewer-sources.js';
import type { ExportFrameValues } from '../../spec.js';
import type { CropBox, ImageSize } from './operations.js';

// The view model of the Image section: the form, the preview and the crop, over `export.preview`, `export.image` and
// `export.image-copy` (the people's API; no plugin). The preview is a small render, asked for after the settings or the
// canvas stop changing, and only while the section shows (`setActive`: a folded section holds no subscription and
// renders nothing; AGENTS.md, hidden UI). The operations serialize on the one shared helper themselves (a queue in
// `ExportApi`), so copy and download simply wait behind a running preview. DOM-free.

/** The form: the screenshot settings, as Mol*'s classic export shows them (the operations' params, grouped). */
/** The Image section's choice: `image`, or the id of an `export.picture` view (an animation). */
const isImage = (v: { mode?: string }) => !v.mode || v.mode === 'image';

export const FormParams = {
  // What the section exports: set by its switch, not a row of the form.
  mode: PD.Text('image', { isHidden: true }),
  resolution: PD.MappedStatic(
    'viewport',
    {
      viewport: PD.Group({}, { isFlat: true }),
      hd: PD.Group({}, { isFlat: true }),
      'full-hd': PD.Group({}, { isFlat: true }),
      'ultra-hd': PD.Group({}, { isFlat: true }),
      '8k-ultra-hd': PD.Group({}, { isFlat: true }),
      custom: PD.Group(
        {
          width: PD.Numeric(1920, { min: 128, max: 8192, step: 1 }),
          height: PD.Numeric(1080, { min: 128, max: 8192, step: 1 }),
        },
        { isFlat: true },
      ),
    },
    {
      options: [
        ['viewport', 'Viewport'],
        ['hd', 'HD (1280 × 720)'],
        ['full-hd', 'Full HD (1920 × 1080)'],
        ['ultra-hd', 'Ultra HD (3840 × 2160)'],
        ['8k-ultra-hd', '8K (7680 × 4320)'],
        ['custom', 'Custom'],
      ],
    },
  ),
  format: PD.MappedStatic(
    'png',
    {
      png: PD.Group({}, { isFlat: true }),
      jpeg: PD.Group({ quality: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }) }, { isFlat: true }),
      webp: PD.Group({ quality: PD.Numeric(0.9, { min: 0, max: 1, step: 0.01 }) }, { isFlat: true }),
    },
    {
      options: [
        ['png', 'PNG'],
        ['jpeg', 'JPEG'],
        ['webp', 'WebP'],
      ],
      hideIf: (v: { mode?: string }) => !isImage(v),
    },
  ),
  transparent: PD.Boolean(false, { hideIf: (v: { mode?: string }) => !isImage(v) }),
  axes: PD.Boolean(false, { label: 'Axes' }),
};

type FormValues = PD.Values<typeof FormParams>;

const Presets: Record<Exclude<ImageSize, 'viewport' | 'custom'>, [number, number]> = {
  hd: [1280, 720],
  'full-hd': [1920, 1080],
  'ultra-hd': [3840, 2160],
  '8k-ultra-hd': [7680, 4320],
};

const Full: CropBox = { x: 0, y: 0, width: 1, height: 1 };

/** The form's values as an image operation's. */
export function imageValues(v: FormValues): Record<string, unknown> {
  const custom = v.resolution.name === 'custom';
  return {
    size: v.resolution.name,
    width: custom ? (v.resolution.params as { width: number }).width : 0,
    height: custom ? (v.resolution.params as { height: number }).height : 0,
    format: v.format.name,
    quality: v.format.name === 'png' ? 0.9 : (v.format.params as { quality: number }).quality,
    transparent: isImage(v) && v.transparent,
    axes: v.axes,
  };
}

/** The frame the form asks for, when it says (the viewport's is the canvas's). */
/** The frame an `export.picture` view renders (`ExportFrame`): the viewport, or its pixels. */
export function frameValues(v: FormValues): ExportFrameValues {
  const r = v.resolution;
  if (r.name === 'viewport') return { size: 'viewport', width: 0, height: 0, axes: v.axes };
  const { width, height } = frameOf(v, undefined)!;
  return { size: 'custom', width, height, axes: v.axes };
}

function frameOf(v: FormValues, viewport: { width: number; height: number } | undefined) {
  const r = v.resolution;
  if (r.name === 'viewport') return viewport;
  if (r.name === 'custom') return r.params as { width: number; height: number };
  const [width, height] = Presets[r.name as keyof typeof Presets];
  return { width, height };
}

export interface ImageExportDeps {
  readonly api: ViewerApi;
  readonly sources: ViewerSources;
  notify(text: string, kind?: MessageKind): void;
  /** Hands the image to the page to save (`ui/download.ts`). */
  save(blob: Blob, filename: string): void;
}

export interface ImageExport extends ImageExportModel {
  /** Starts (true) or stops (false) the preview's subscriptions and work: the view's activity. Starts active. */
  setActive(active: boolean): void;
  dispose(): void;
}

const PreviewDelayMs = 250;
const FallbackAspect = 16 / 9;

export function createImageExportModel(deps: ImageExportDeps): ImageExport {
  const form = new ParamsFormModel(FormParams);
  const autoCrop = observableState(true);
  const manual = observableState<CropBox>({ ...Full });
  // The crop the content has, from the last preview.
  const found = observableState<CropBox>({ ...Full });
  const crop = observableState<CropBox>({ ...Full });
  const preview = observableState<string | undefined>(undefined);
  const viewport = observableState<{ width: number; height: number } | undefined>(undefined);
  const aspect = observableState(FallbackAspect);
  const size = observableState('');
  const busy = observableState(false);
  // Whether the preview is of the settings now: false from a change until its preview is in.
  const settled = observableState(false);

  let disposed = false;
  let token = 0;
  // The preview asked for last: a newer one stops an older that has not started (a running render cannot be stopped).
  let asked: AbortController | undefined;
  // Previews come from the script origin: a failed one is not worth a message of the person's.
  const quiet = deps.api.as({ kind: 'script' });

  const compute = () => {
    const v = form.values.value;
    const frame = frameOf(v, viewport.value);
    // An animation renders the whole frame: a crop found for one image does not fit a moving camera.
    const still = isImage(v);
    const c = !still ? Full : autoCrop.value ? found.value : manual.value;
    crop.set(c);
    aspect.set(frame ? frame.width / frame.height : FallbackAspect);
    const px = frame ? `${Math.round(frame.width * c.width)} × ${Math.round(frame.height * c.height)} px` : '';
    size.set(frame && still ? `${px} · ${v.format.name.toUpperCase()}` : px);
  };

  const refresh = async () => {
    if (disposed || !deps.sources.canvas.value || !quiet.export) return;
    const mine = ++token;
    asked?.abort();
    const mineAsked = (asked = new AbortController());
    const r = await quiet.export.preview(
      { ...imageValues(form.values.value), crop: autoCrop.value && isImage(form.values.value) ? 'auto' : 'none' },
      { signal: mineAsked.signal },
    );
    if (disposed || mine !== token || !r.ok || !r.value) return;
    const { dataUri, width, height, crop: box } = r.value;
    preview.set(`center / 100% 100% no-repeat url("${dataUri}")`);
    viewport.set({ width, height });
    if (autoCrop.value && isImage(form.values.value)) found.set(box);
    settled.set(true);
  };

  // Turning auto-crop off starts a manual crop from the whole frame (not from the content's box, which would look as if
  // nothing changed).
  let wasAuto = autoCrop.value;
  let subs: Subscription[] = [];
  const start = () => {
    if (disposed || subs.length) return;
    subs = [
      autoCrop.asReadonly.subscribe((auto) => {
        if (!auto && wasAuto) manual.set({ ...Full });
        wasAuto = auto;
      }),
      combineLatest([
        form.values,
        autoCrop.asReadonly,
        manual.asReadonly,
        found.asReadonly,
        viewport.asReadonly,
      ]).subscribe(compute),
      combineLatest([form.values, autoCrop.asReadonly, deps.sources.canvas, deps.sources.structures])
        .pipe(
          // The preview no longer matches from the change on, not from when it is asked for.
          tap(() => settled.set(false)),
          debounceTime(PreviewDelayMs),
        )
        .subscribe(() => void refresh()),
    ];
  };
  const stop = () => {
    token++;
    asked?.abort();
    asked = undefined;
    for (const s of subs.splice(0)) s.unsubscribe();
  };

  /** Runs an operation with the form's values as the person (a failure says itself); `done` handles the result. */
  const run = async <R>(op: 'image' | 'imageCopy', done: (value: R | undefined, did: string) => void) => {
    const fn = deps.api.export?.[op] as ((v: never) => Promise<{ ok: boolean; did: string; value?: R }>) | undefined;
    if (!fn || busy.value) return;
    busy.set(true);
    try {
      const box = autoCrop.value ? { crop: 'auto' } : { crop: 'none', box: { ...manual.value } };
      const r = await fn({ ...imageValues(form.values.value), ...box } as never);
      if (r.ok) done(r.value, r.did);
    } finally {
      busy.set(false);
    }
  };

  compute();
  start();

  return {
    // `preview` is a CSS background value: the image, under the crop box.
    preview: preview.asReadonly,
    aspect: aspect.asReadonly,
    autoCrop,
    settled: settled.asReadonly,
    crop: crop.asReadonly,
    setCrop: (c) => {
      autoCrop.set(false);
      manual.set(c);
    },
    size: size.asReadonly,
    form,
    busy: busy.asReadonly,
    copy: () => void run('imageCopy', () => deps.notify('Copied the image.', 'success')),
    download: () =>
      void run<{ blob: Blob; filename: string; width: number; height: number }>('image', (image) => {
        if (!image) return;
        deps.save(image.blob, image.filename);
        deps.notify(`Saved ${image.filename} (${image.width} × ${image.height}).`, 'success');
      }),
    setActive: (active) => (active ? start() : stop()),
    dispose() {
      disposed = true;
      stop();
    },
  };
}
