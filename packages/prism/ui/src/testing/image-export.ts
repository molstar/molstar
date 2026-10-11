/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { combineLatest, distinctUntilChanged, map, type Observable } from 'rxjs';
import { ParamsFormModel } from '../forms/params-form.js';
import type { CropBox, ImageExportModel } from '../model/image-export.js';
import { observableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

/** A state following a stream (the mocks live as long as the story). */
function readonlyState<T>(source: Observable<T>): ReadonlyState<T> {
  const state = observableState<T>(undefined as T);
  source.subscribe((v) => state.set(v));
  return state.asReadonly;
}

// Mol*'s ViewportScreenshotHelper params.
const Params = {
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
    },
  ),
  transparent: PD.Boolean(false),
  axes: PD.MappedStatic(
    'off',
    {
      on: PD.Group(
        {
          location: PD.Select('bottom-left', [
            ['bottom-left', 'Bottom left'],
            ['bottom-right', 'Bottom right'],
            ['top-left', 'Top left'],
            ['top-right', 'Top right'],
          ]),
          scale: PD.Numeric(0.33, { min: 0.1, max: 2, step: 0.01 }),
          labels: PD.Boolean(true),
        },
        { isFlat: true },
      ),
      off: PD.Group({}, { isFlat: true }),
    },
    { cycle: true, label: 'Axes' },
  ),
  illumination: PD.Group(
    {
      extraIterations: PD.Numeric(1, { min: 0, max: 5, step: 1 }, { label: 'Extra iterations' }),
      targetIterationTimeMs: PD.Numeric(300, { min: 100, max: 3000, step: 10 }, { label: 'Iteration time (ms)' }),
    },
    { label: 'Illumination' },
  ),
};

const Sizes: Record<string, [number, number]> = {
  viewport: [1366, 812],
  hd: [1280, 720],
  'full-hd': [1920, 1080],
  'ultra-hd': [3840, 2160],
  '8k-ultra-hd': [7680, 4320],
};

const AutoCrop: CropBox = { x: 0.17, y: 0.1, width: 0.66, height: 0.8 };

/** Image export for the stories: the params, a drawn preview, the crop. */
export function mockImageExport(post: (text: string, kind?: 'info' | 'success') => void): ImageExportModel {
  const form = new ParamsFormModel(Params);
  const autoCrop = observableState(true);
  const manual = observableState<CropBox>({ x: 0.05, y: 0.05, width: 0.9, height: 0.9 });
  const busy = observableState(false);
  const dims = form.values.pipe(
    map((v: any): [number, number] =>
      v.resolution.name === 'custom'
        ? [v.resolution.params.width, v.resolution.params.height]
        : Sizes[v.resolution.name],
    ),
  );
  const crop = combineLatest([autoCrop.asReadonly, manual.asReadonly]).pipe(map(([a, m]) => (a ? AutoCrop : m)));
  const background = form.values.pipe(
    map((v: any) =>
      v.transparent
        ? 'repeating-conic-gradient(#d8d8de 0 25%, #f2f2f5 0 50%) 0 0 / 12px 12px'
        : 'radial-gradient(ellipse at 50% 45%, #fff 0, #eef0f4 70%)',
    ),
  );
  const molecule =
    'radial-gradient(circle at 34% 62%, #6aa6d8 0 6%, transparent 6.5%), radial-gradient(circle at 46% 54%, #7a6fd6 0 7%, transparent 7.5%), radial-gradient(circle at 58% 40%, #b46ad0 0 8%, transparent 8.5%), radial-gradient(circle at 68% 52%, #d06a8c 0 6%, transparent 6.5%)';
  return {
    preview: readonlyState(background.pipe(map((b) => `${molecule}, ${b}`))),
    aspect: readonlyState(
      dims.pipe(
        map(([w, h]) => w / h),
        distinctUntilChanged(),
      ),
    ),
    autoCrop,
    crop: readonlyState(crop),
    setCrop: (c) => {
      autoCrop.set(false);
      manual.set(c);
    },
    size: readonlyState(
      combineLatest([dims, crop, form.values]).pipe(
        map(([[w, h], c, v]: any) => {
          const cw = Math.round(w * c.width);
          const ch = Math.round(h * c.height);
          return `${cw} × ${ch} px · ${String(v.format.name).toUpperCase()}`;
        }),
      ),
    ),
    form,
    busy: busy.asReadonly,
    copy: () => post('Copied the image', 'success'),
    download: () => {
      busy.set(true);
      setTimeout(() => {
        busy.set(false);
        post(`Saved molstar-image.${(form.values.value as any).format.name}`, 'success');
      }, 700);
    },
  };
}
