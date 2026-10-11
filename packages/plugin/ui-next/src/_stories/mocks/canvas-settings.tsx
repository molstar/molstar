/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ParamsForm } from '../../components/params/params-form.js';
import { domView } from '../../dom/view.js';
import { ParamsFormModel } from '../../forms/params-form.js';
import type { ShellModel } from '../../shell/shell.js';
import { SectionList } from '../../views/shell/section-list.js';

// Canvas3D's params, grouped by what people reach for: the background, the lighting and effects, the camera, clipping
// and fog. Shaped like Mol*'s (renderer, postprocessing, camera, cameraFog, cameraClipping).
const Background = {
  color: PD.Color(0xffffff as any, { isEssential: true }),
  transparent: PD.Boolean(false, { label: 'Transparent', isEssential: true }),
};

const Lighting = {
  style: PD.Select(
    'default',
    [
      ['default', 'Default'],
      ['plastic', 'Plastic'],
      ['flat', 'Flat'],
      ['celshaded', 'Cel-shaded'],
      ['illustrative', 'Illustrative'],
    ],
    PD.Essential,
  ),
  // Mol*'s postprocessing pattern: on (with its settings) or off.
  occlusion: PD.MappedStatic(
    'on',
    {
      on: PD.Group({
        samples: PD.Numeric(32, { min: 1, max: 256, step: 1 }),
        multiScale: PD.Boolean(false, { label: 'Multi scale' }),
        radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }),
        bias: PD.Numeric(0.8, { min: 0, max: 3, step: 0.1 }),
        blurKernelSize: PD.Numeric(15, { min: 1, max: 25, step: 2 }, { label: 'Blur kernel size' }),
        blurDepthBias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }, { label: 'Blur depth bias' }),
      }),
      off: PD.EmptyGroup(),
    },
    { label: 'Ambient occlusion', isEssential: true },
  ),
  shadow: PD.MappedStatic(
    'off',
    {
      on: PD.Group({
        steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }),
        maxDistance: PD.Numeric(3, { min: 0, max: 256, step: 1 }, { label: 'Max distance' }),
        tolerance: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
      }),
      off: PD.EmptyGroup(),
    },
    { label: 'Shadows', isEssential: true },
  ),
  outline: PD.MappedStatic(
    'off',
    {
      on: PD.Group({
        scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
        threshold: PD.Numeric(0.33, { min: 0.01, max: 1, step: 0.01 }),
        color: PD.Color(0x000000 as any),
        includeTransparent: PD.Boolean(true, { label: 'Include transparent' }),
      }),
      off: PD.EmptyGroup(),
    },
    { label: 'Outline', isEssential: true },
  ),
  exposure: PD.Numeric(1, { min: 0, max: 3, step: 0.01 }),
};

const Camera = {
  mode: PD.Select(
    'perspective',
    [
      ['perspective', 'Perspective'],
      ['orthographic', 'Orthographic'],
    ],
    PD.Essential,
  ),
  fov: PD.Numeric(45, { min: 10, max: 130, step: 1 }, { label: 'Field of view', isEssential: true }),
  stereo: PD.Boolean(false, { label: 'Stereo' }),
  helper: PD.Boolean(false, { label: 'Axes helper' }),
};

const ClippingFog = {
  fog: PD.Numeric(50, { min: 0, max: 100, step: 1 }, { label: 'Fog intensity', isEssential: true }),
  radius: PD.Numeric(100, { min: 0, max: 99, step: 1 }, { label: 'Clip radius', isEssential: true }),
  far: PD.Boolean(true, { label: 'Clip far plane' }),
};

/** The Canvas settings card (canvas stack ⚙︎): sections like the Export panel's, so extensions can add their own. */
export function registerMockCanvasSettings(shell: ShellModel) {
  shell.registerOverlays([
    {
      id: 'canvas-settings',
      title: 'Canvas',
      render: domView(({ ctx }) => <SectionList shell={ctx.shell} slot="canvas-settings" focus={ctx.focus} />),
    },
  ]);
  const section = (id: string, title: string, icon: any, params: PD.Params, order: number) => {
    const form = new ParamsFormModel(params);
    return { id, slot: 'canvas-settings', title, icon, order, render: domView(() => <ParamsForm form={form} />) };
  };
  shell.registerSections([
    section('background', 'Background', 'palette', Background, 0),
    section('lighting', 'Lighting & effects', 'lightbulb', Lighting, 1),
    section('camera', 'Camera', 'camera', Camera, 2),
    section('clipping', 'Clipping & fog', 'scan', ClippingFog, 3),
  ]);
}
