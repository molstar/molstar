/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../_stories/story-root.js';
import { LucideIcons } from './_generated/lucide.js';
import { Icon, type IconName } from './icon.js';
import { MolecularIcons } from './molecular.js';

const meta: Meta = { title: 'Components/Icon' };
export default meta;

const grid =
  'display: grid; grid-template-columns: repeat(auto-fill, minmax(112px, 1fr)); gap: var(--msp-space-2); padding: var(--msp-space-4);';
const cell =
  'display: flex; flex-direction: column; align-items: center; gap: var(--msp-space-2); padding: var(--msp-space-3); border: 1px solid var(--msp-color-border); border-radius: var(--msp-radius-md); color: var(--msp-color-text); background: var(--msp-color-panel);';
const sizes = 'display: flex; align-items: center; gap: var(--msp-space-3);';
const caption =
  'font-size: var(--msp-font-size-xs); color: var(--msp-color-text-muted); text-align: center; word-break: break-all;';
const heading =
  'margin: 0; padding: var(--msp-space-4) var(--msp-space-4) 0; color: var(--msp-color-text); font-size: var(--msp-font-size-md);';

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
      <div style={`${sizes} padding: var(--msp-space-4); color: var(--msp-color-accent);`}>
        <Icon name="eye" label="Visible" />
        <Icon name="cartoon" size={32} label="Cartoon representation" />
        <span style="color: var(--msp-color-danger)">
          <Icon name="triangle-alert" size={20} />
        </span>
      </div>
    )),
};
