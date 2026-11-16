/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import type { Child } from '../../../dom/index.js';
import { onCleanup } from '../../../dom/index.js';
import { Note, page, story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import type { IconName } from '../../icon/icon.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { Toolbar, type ToolbarFace, type ToolbarGroup, type ToolbarMenuRow } from '../toolbar.js';

const meta: Meta = { title: 'Components/Toolbar' };
export default meta;

/** A story's frame: tooltips on, and what the last button or row did. */
function Frame(props: { children: (ran: (what: string) => void) => Child; note: string; area?: string }) {
  const ran = state('nothing yet');
  return (
    <div
      style={page}
      ref={(el: HTMLElement) => {
        onCleanup(installTooltips(el));
      }}
    >
      <Note>{props.note}</Note>
      <div style={props.area ?? ''}>{props.children((what) => ran.set(what))}</div>
      <Note>ran: {ran}</Note>
    </div>
  );
}

/** A single choice's rows: each pressed while it is the value. */
function choice(
  value: State<string>,
  group: string,
  options: readonly (readonly [id: string, label: string, icon?: IconName])[],
  heading?: string,
): ToolbarMenuRow[] {
  return options.map(([id, label, icon], i) => ({
    id,
    label,
    ...(icon ? { icon } : {}),
    group,
    ...(heading && i === 0 ? { heading } : {}),
    choice: true,
    pressed: derive(value, (v) => v === id),
    run: () => value.set(id),
  }));
}

/** The viewer's bar: camera and picking, view, lighting and style; the picking level is a menu of two groups. */
function viewerGroups(ran: (what: string) => void): ToolbarGroup[] {
  const spin = state(false);
  const level = state('residue');
  const instances = state('off');
  const lighting = state('studio');
  const ghost = state(false);
  const levelIcon = derive(level, (l): ToolbarFace => {
    const icon = `pick-${l}` as IconName;
    return { icon, title: `Pick: ${l}` };
  });
  return [
    {
      id: 'interact',
      label: 'Interaction',
      items: [
        {
          id: 'rotate',
          icon: 'rotate-3d',
          title: 'Spin',
          pressed: spin,
          run: () => {
            spin.set(!spin.value);
            ran(`spin ${spin.value ? 'on' : 'off'}`);
          },
          menu: [
            { id: 'slow', label: 'Slow', hint: '0.5×', group: 'speed', run: () => ran('spin slow') },
            { id: 'normal', label: 'Normal', hint: '1×', group: 'speed', run: () => ran('spin normal') },
            { id: 'fast', label: 'Fast', hint: '2×', group: 'speed', run: () => ran('spin fast') },
            { id: 'rock', label: 'Rock instead', icon: 'orbit', group: 'mode', run: () => ran('rock') },
          ],
        },
        {
          id: 'select',
          icon: 'mouse-pointer-2',
          title: 'Picking level',
          face: levelIcon,
          menu: [
            ...choice(level, 'level', [
              ['atom', 'Atom', 'pick-atom'],
              ['residue', 'Residue', 'pick-residue'],
              ['chain', 'Chain', 'pick-chain'],
              ['entity', 'Entity', 'pick-entity'],
              ['model', 'Model', 'pick-model'],
            ]),
            ...choice(
              instances,
              'instances',
              [
                ['off', 'Only the picked one'],
                ['same-residue', 'Same residue name'],
                ['same-entity', 'Same entity'],
              ],
              'All instances',
            ),
          ],
        },
      ],
    },
    {
      id: 'camera',
      label: 'Camera',
      items: [
        { id: 'focus', icon: 'focus', title: 'Focus on the selection', run: () => ran('focus') },
        { id: 'expand', icon: 'maximize-2', title: 'Expand', run: () => ran('expand') },
      ],
    },
    {
      id: 'look',
      label: 'Look',
      items: [
        {
          id: 'lighting',
          icon: 'lightbulb',
          title: 'Lighting',
          menu: choice(lighting, 'preset', [
            ['flat', 'Flat'],
            ['studio', 'Studio'],
            ['illustrative', 'Illustrative'],
            ['occlusion', 'Ambient occlusion'],
          ]),
        },
        {
          id: 'style',
          icon: 'sparkles',
          title: 'More style',
          menu: [
            {
              id: 'ghost',
              label: 'Ghost the rest',
              icon: 'eye-off',
              group: 'view',
              pressed: ghost,
              run: () => ghost.set(!ghost.value),
            },
            {
              id: 'outline',
              label: 'Outline',
              icon: 'paintbrush',
              group: 'view',
              items: [
                { id: 'thin', label: 'Thin', run: () => ran('outline thin') },
                { id: 'thick', label: 'Thick', run: () => ran('outline thick') },
                { id: 'none', label: 'None', run: () => ran('outline none') },
              ],
            },
            { id: 'export', label: 'Export image', icon: 'image', hint: '⌘E', group: 'out', run: () => ran('export') },
            { id: 'xr', label: 'Enter XR', icon: 'glasses', group: 'out', disabled: 'No XR device', run: () => {} },
          ],
        },
      ],
    },
  ];
}

/**
 * The viewer's bar: blocks of square 32px buttons with 16px icons, the gaps between blocks the grouping. Buttons with a
 * menu have a corner mark: hover opens the menu (it stays while the pointer is over it), so does ↓ from the keyboard.
 * Picking level: a choice of five (the chosen one in the accent) and a second choice headed "All instances"; its icon
 * follows the level. Spin is a toggle with a menu: a click runs it, hover shows its speeds. More style has a check, a
 * submenu (→), a hint and a disabled row (its reason in the tooltip).
 */
export const Horizontal: StoryObj = {
  render: () =>
    story(() => (
      <Frame note="Tab to the bar, ← → between buttons, ↓ opens a menu; or hover a button with a corner mark.">
        {(ran) => <Toolbar groups={viewerGroups(ran)} aria-label="Viewer" />}
      </Frame>
    )),
};

/** Over a canvas, at its bottom edge: raised blocks, menus open above (flipped below where there is no room). */
export const Above: StoryObj = {
  render: () =>
    story(() => (
      <Frame
        note="The bar at the bottom of a canvas; its menus open upward."
        area="display: flex; align-items: flex-end; justify-content: flex-end; height: 360px; padding: var(--msui-space-3); border-radius: var(--msui-radius-lg); background: var(--msui-surface-sunken); box-shadow: inset 0 0 0 1px var(--msui-line);"
      >
        {(ran) => <Toolbar groups={viewerGroups(ran)} menuSide="above" elevated aria-label="Viewer" />}
      </Frame>
    )),
};

/**
 * One surface with separators (`grouping="separators"`): the same groups on one block, a hairline between them, in a
 * row and in a column; the arrows and menus as with blocks.
 */
export const Separators: StoryObj = {
  render: () =>
    story(() => (
      <Frame note="One block for the bar, a hairline between its groups.">
        {(ran) => (
          <div style="display: flex; gap: var(--msui-space-5); align-items: flex-start;">
            <Toolbar groups={viewerGroups(ran)} grouping="separators" aria-label="Viewer" />
            <Toolbar
              groups={viewerGroups(ran)}
              grouping="separators"
              orientation="vertical"
              aria-label="Viewer column"
            />
          </div>
        )}
      </Frame>
    )),
};

/** A column: ↓ ↑ move, → opens a menu (to the right, flipped at the window edge), ← in it closes it. */
export const Vertical: StoryObj = {
  render: () =>
    story(() => (
      <Frame note="A vertical bar: menus open beside it, the first row level with the button.">
        {(ran) => <Toolbar groups={viewerGroups(ran)} orientation="vertical" aria-label="Viewer" />}
      </Frame>
    )),
};

/**
 * States: toggles (pressed), captions (beside the icon in a row, under it in a column), a button disabled with a reason
 * (it keeps its place for the arrows; the reason is in its tooltip), one that hides, and a face that changes (the
 * theme: icon and name).
 */
export const States: StoryObj = {
  render: () =>
    story(() => (
      <Frame note="Toggle the axes and the grid; the theme button changes its face; Measure is disabled until Select is on.">
        {(ran) => {
          const axes = state(true);
          const grid = state(false);
          const select = state(false);
          const extra = state(true);
          const theme = state<'light' | 'dark'>('light');
          const face = derive(
            theme,
            (t): ToolbarFace =>
              t === 'light' ? { icon: 'sun', title: 'Theme: light' } : { icon: 'moon', title: 'Theme: dark' },
          );
          const measureOff = derive(select, (s) => (s ? false : 'Turn on Select first'));
          const groups = (captioned: boolean): ToolbarGroup[] => [
            {
              id: 'toggles',
              items: [
                {
                  id: 'axes',
                  icon: 'axis-3d',
                  title: 'Axes',
                  pressed: axes,
                  ...(captioned ? { caption: 'Axes' } : {}),
                  run: () => axes.set(!axes.value),
                },
                {
                  id: 'grid',
                  icon: 'box',
                  title: 'Box',
                  pressed: grid,
                  ...(captioned ? { caption: 'Box' } : {}),
                  run: () => grid.set(!grid.value),
                },
              ],
            },
            {
              id: 'tools',
              items: [
                {
                  id: 'select',
                  icon: 'mouse-pointer-click',
                  title: 'Select',
                  pressed: select,
                  ...(captioned ? { caption: 'Select' } : {}),
                  run: () => select.set(!select.value),
                },
                {
                  id: 'measure',
                  icon: 'ruler',
                  title: 'Measure',
                  disabled: measureOff,
                  ...(captioned ? { caption: 'Measure' } : {}),
                  run: () => ran('measure'),
                },
                {
                  id: 'extra',
                  icon: 'scissors',
                  title: 'Clip (hides with the theme dark)',
                  visible: extra,
                  run: () => ran('clip'),
                },
              ],
            },
            {
              id: 'app',
              items: [
                {
                  id: 'theme',
                  icon: 'sun',
                  title: 'Theme',
                  face,
                  ...(captioned ? { caption: '⌘T' } : {}),
                  run: () => {
                    theme.set(theme.value === 'light' ? 'dark' : 'light');
                    extra.set(theme.value === 'light');
                  },
                },
              ],
            },
          ];
          return (
            <div style="display: flex; gap: var(--msui-space-6); align-items: flex-start;">
              <div style="display: grid; gap: var(--msui-space-4);">
                <Toolbar groups={groups(false)} aria-label="Plain" />
                <Toolbar groups={groups(true)} aria-label="Captioned" />
              </div>
              <Toolbar groups={groups(true)} orientation="vertical" menuSide="left" aria-label="Captioned column" />
            </div>
          );
        }}
      </Frame>
    )),
};
