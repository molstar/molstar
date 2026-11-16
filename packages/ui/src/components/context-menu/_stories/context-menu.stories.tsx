/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { onCleanup } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { type ContextMenuItem, openContextMenu } from '../context-menu.js';

const meta: Meta = { title: 'Components/Context menu' };
export default meta;

const area =
  'display: grid; place-items: center; height: 240px; margin: var(--msui-space-4); border-radius: var(--msui-radius-lg); background: var(--msui-surface); color: var(--msui-text-muted);';

function Area(props: { items: (ran: (what: string) => void) => readonly ContextMenuItem[]; hint: string }) {
  const ran = state('nothing yet');
  let root!: HTMLElement;
  return (
    <div
      ref={(el: HTMLElement) => {
        root = el;
        onCleanup(installTooltips(el));
      }}
    >
      <div
        style={area}
        onContextMenu={(e: MouseEvent) => {
          e.preventDefault();
          openContextMenu(
            root,
            { x: e.clientX, y: e.clientY },
            props.items((what) => ran.set(what)),
          );
        }}
      >
        {props.hint}
      </div>
      <p style="margin: 0 var(--msui-space-4); color: var(--msui-text-muted);">ran: {ran}</p>
    </div>
  );
}

/**
 * A right-click opens it at the pointer (flipped at the window's edges). Rows are list rows: an icon slot on the inset
 * (kept empty for an item without an icon), the label after it, the key ending on the inner padding. A disabled item
 * shows why in its tooltip; the arrows skip it.
 */
export const WithIcons: StoryObj = {
  render: () =>
    story(() => (
      <Area
        hint="Right-click here"
        items={(ran) => [
          { id: 'focus', label: 'Focus', icon: 'scan', shortcut: 'F', run: () => ran('Focus') },
          { id: 'hide', label: 'Hide', icon: 'eye-off', shortcut: 'H', run: () => ran('Hide') },
          { id: 'snapshot', label: 'Take a snapshot', icon: 'camera', run: () => ran('Snapshot') },
          { id: 'measure', label: 'Measure', run: () => ran('Measure') },
          { id: 'remove', label: 'Remove', icon: 'trash-2', disabled: 'Nothing selected', run: () => ran('Remove') },
        ]}
      />
    )),
};

/** Without icons: no slot, the labels on the inset. */
export const Plain: StoryObj = {
  render: () =>
    story(() => (
      <Area
        hint="Right-click here"
        items={(ran) => [
          { id: 'copy', label: 'Copy', shortcut: '⌘C', run: () => ran('Copy') },
          { id: 'paste', label: 'Paste', shortcut: '⌘V', run: () => ran('Paste') },
          { id: 'select-all', label: 'Select all', shortcut: '⌘A', run: () => ran('Select all') },
        ]}
      />
    )),
};
