/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { onCleanup } from '../../dom/index.js';
import { ParamsFormModel } from '../params-form.js';
import { type SettingsGroup, SettingsView } from '../settings-view.js';

// Settings over generic param groups: a search over all of them, each group an `Expandable` (styles §8).

const meta: Meta = { title: 'Views/Settings' };
export default meta;

const groups = (): SettingsGroup[] => [
  {
    id: 'appearance',
    title: 'Appearance',
    icon: 'palette',
    form: new ParamsFormModel({
      background: PD.Color(0xfcfbf9 as any, { label: 'Background' }),
      transparent: PD.Boolean(false, { label: 'Transparent background' }),
      outline: PD.Boolean(true, { label: 'Outline', description: 'A line around every object' }),
      occlusion: PD.Numeric(0.6, { min: 0, max: 1, step: 0.01 }, { label: 'Ambient occlusion' }),
    }),
  },
  {
    id: 'camera',
    title: 'Camera',
    icon: 'camera',
    form: new ParamsFormModel({
      mode: PD.Select('perspective', [
        ['perspective', 'Perspective'],
        ['orthographic', 'Orthographic'],
      ]),
      fov: PD.Numeric(45, { min: 10, max: 90, step: 1 }, { label: 'Field of view' }),
      clip: PD.Interval([1, 100], { min: 1, max: 100, step: 1 }, { label: 'Clipping' }),
    }),
  },
  {
    id: 'interaction',
    title: 'Interaction',
    icon: 'mouse-pointer-2',
    form: new ParamsFormModel({
      highlight: PD.Boolean(true, { label: 'Highlight on hover' }),
      granularity: PD.Select('residue', [
        ['atom', 'Atom'],
        ['residue', 'Residue'],
        ['chain', 'Chain'],
      ]),
      zoomSpeed: PD.Numeric(1, { min: 0.1, max: 4, step: 0.1 }, { label: 'Zoom speed' }),
    }),
  },
  {
    id: 'performance',
    title: 'Performance',
    icon: 'gauge',
    form: new ParamsFormModel({
      quality: PD.Select('auto', [
        ['auto', 'Auto'],
        ['low', 'Low'],
        ['high', 'High'],
      ]),
      pixelScale: PD.Numeric(1, { min: 0.5, max: 2, step: 0.25 }, { label: 'Pixel scale' }),
    }),
  },
];

/** Groups of settings, the first open; type in the search: matching groups open, the rest hide. */
export const Settings: StoryObj = {
  render: () =>
    story(() => (
      <div
        ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
        class="msui-surface"
        style="width: 320px; margin: var(--msui-space-4); padding: 6px;"
      >
        <SettingsView groups={groups()} />
      </div>
    )),
};
