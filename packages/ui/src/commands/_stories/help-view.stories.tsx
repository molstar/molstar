/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/support.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { onCleanup } from '../../dom/index.js';
import { HelpView, type InputBindingGroup } from '../help-view.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition } from '../types.js';

// Help over a demo command set (commands §6.6): an app's mouse and touch bindings, the keys everywhere, the commands'
// keys, a conflict; "Change keys" records, clears and resets (remaps kept for the story's life).

const meta: Meta = { title: 'Views/Help' };
export default meta;

const command = (
  id: string,
  title: string,
  category: string,
  shortcut?: CommandDefinition<object>['shortcut'],
): CommandDefinition<object> => ({ id, title, category, ...(shortcut ? { shortcut } : {}), run: () => ({}) });

const bindings: readonly InputBindingGroup[] = [
  {
    title: 'Camera',
    bindings: [
      { action: 'Rotate', input: ['Drag'], or: ['One finger'] },
      { action: 'Move', input: ['Right drag'], or: ['Two fingers'] },
      { action: 'Zoom', input: ['Wheel'], or: ['Pinch'] },
      { action: 'Focus', input: ['Click'] },
    ],
  },
  {
    title: 'Selection',
    bindings: [
      { action: 'Select', input: ['Shift', 'Click'] },
      { action: 'Add a range', input: ['Shift', 'Alt', 'Click'] },
    ],
  },
];

function manager(): CommandManager<object> {
  const m = new CommandManager<object>({
    context: new State<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  m.register([
    command('menu.open', 'Command menu', 'App', { key: 'k', mod: true }),
    command('help.open', 'Help', 'App', { key: '?' }),
    command('view.reset', 'Reset the camera', 'View', { key: 'r', alt: true }),
    command('view.focus', 'Focus on the selection', 'View', { key: 'f' }),
    command('step.next', 'Next step', 'Playback', { key: 'ArrowRight' }),
    command('step.previous', 'Previous step', 'Playback', { key: 'ArrowLeft' }),
    // Two on one key: a conflict, said first.
    command('look.cartoon', 'Show as cartoon', 'Look', { key: 'c', alt: true }),
    command('look.color', 'Color by chain', 'Look', { key: 'c', alt: true }),
    command('export.image', 'Save an image', 'Export'),
    command('export.session', 'Save the session', 'Export'),
  ]);
  return m;
}

/** Help in a view's width: bindings, the keys everywhere, the commands with keys; Change keys edits them. */
export const Help: StoryObj = {
  render: () =>
    story(() => (
      <div
        ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
        class="msui-surface"
        style="width: 340px; margin: var(--msui-space-4); padding: var(--msui-space-1) 6px 6px;"
      >
        <HelpView
          commands={manager()}
          bindings={bindings}
          global={[
            { action: 'Command menu', input: ['⌘', 'K'] },
            { action: 'Help', input: ['?'] },
            { action: 'Close a menu or a card', input: ['Esc'] },
          ]}
        />
      </div>
    )),
};
