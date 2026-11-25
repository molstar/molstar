/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { LucideIcons } from '../_generated/lucide.js';
import { Icon, type IconName } from '../icon.js';
import { MolecularIcons } from '../molecular.js';

const meta: Meta = { title: 'Components/Icon' };
export default meta;

const grid =
  'display: grid; grid-template-columns: repeat(auto-fill, minmax(112px, 1fr)); gap: var(--msui-space-1); padding: var(--msui-space-3);';
const cell =
  'display: flex; flex-direction: column; align-items: center; gap: var(--msui-space-1); padding: var(--msui-space-2); border: 1px solid var(--msui-line); border-radius: var(--msui-radius-md); color: var(--msui-text); background: var(--msui-surface);';
const sizes = 'display: flex; align-items: center; gap: var(--msui-space-2);';
const caption =
  'font-size: var(--msui-text-xs); color: var(--msui-text-muted); text-align: center; word-break: break-all;';
const heading =
  'margin: 0; padding: var(--msui-space-3) var(--msui-space-3) 0; color: var(--msui-text); font-size: var(--msui-text-sm);';

function Gallery(props: { names: readonly IconName[] }) {
  return (
    <div style={grid}>
      {props.names.map((name) => (
        <div style={cell} title={name}>
          <div style={sizes}>
            <Icon name={name} size={16} />
            <Icon name={name} size={24} />
          </div>
          <span style={caption}>{name}</span>
        </div>
      ))}
    </div>
  );
}

export const All: StoryObj = {
  render: () =>
    story(() => (
      <div>
        <h3 style={heading}>Molecular ({Object.keys(MolecularIcons).length})</h3>
        <Gallery names={Object.keys(MolecularIcons) as IconName[]} />
        <h3 style={heading}>Lucide ({Object.keys(LucideIcons).length})</h3>
        <Gallery names={Object.keys(LucideIcons) as IconName[]} />
      </div>
    )),
};

export const Labelled: StoryObj = {
  render: () =>
    story(() => (
      <div style={`${sizes} padding: var(--msui-space-3); color: var(--msui-accent);`}>
        <Icon name="eye" label="Visible" />
        <Icon name="cartoon" size={32} label="Cartoon representation" />
        <span style="color: var(--msui-danger)">
          <Icon name="triangle-alert" size={20} />
        </span>
      </div>
    )),
};
