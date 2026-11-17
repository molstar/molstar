/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { Button } from '../../components/button/button.js';
import { Icon } from '../../components/icon/icon.js';
import { Heading, Note, page, story } from '../support.js';

const meta: Meta = { title: 'Foundations/Icons' };
export default meta;

const Sizes = [
  ['--msui-icon-sm', 'sm', 'Commands: in 20px controls, menus and option lists; with small text'],
  ['--msui-icon', 'sm', 'Objects: in list and tree rows, headers, default buttons'],
  ['--msui-icon-lg', 'lg', 'With large text; empty states'],
] as const;

/**
 * Icon sizing (styles §4.1): three sizes on the grid, each paired with its text; an icon takes its size from where it
 * sits (`--msui-icon-size`, set by a button, a slot, a header), so a container and its icons never drift apart.
 */
export const Sizing: StoryObj = {
  render: () =>
    story(() => (
      <div style={page}>
        <Heading>Sizes and their text</Heading>
        <div style="display: grid; gap: var(--msui-space-3);">
          {Sizes.map(([token, text, use]) => (
            <div style="display: grid; grid-template-columns: 160px 280px minmax(0, 1fr); column-gap: var(--msui-space-4); align-items: center;">
              <code style="font-size: var(--msui-text-xs);">{token}</code>
              <span
                style={{
                  display: 'inline-flex',
                  alignItems: 'center',
                  gap: 'var(--msui-space-1)',
                  fontSize: `var(--msui-text-${text})`,
                  lineHeight: `var(--msui-leading-${text})`,
                  '--msui-icon-size': `var(${token})`,
                }}
              >
                <Icon name="layers" />
                Structure components
              </span>
              <span style="color: var(--msui-text-muted)">{use}</span>
            </div>
          ))}
        </div>
        <Heading>From the context</Heading>
        <Note>
          The same icon, no size given: a row's slot, a header, a menu item and each button size set it. An icon-only
          button is square at its height.
        </Note>
        <div style="display: grid; grid-template-columns: 280px minmax(0, 1fr); gap: var(--msui-space-5); align-items: start;">
          <div class="msui-surface">
            <div class="msui-header">
              <span class="msui-slot">
                <Icon name="layers" />
              </span>
              Header
            </div>
            <div class="msui-row">
              <span class="msui-slot">
                <Icon name="layers" />
              </span>
              <span>A row</span>
              <span />
            </div>
          </div>
          <div style="display: grid; gap: var(--msui-space-3);">
            {(['md', 'sm', 'xs'] as const).map((size) => (
              <div style="display: flex; gap: var(--msui-space-2); align-items: center;">
                <code style="width: 32px; font-size: var(--msui-text-xs);">{size}</code>
                <Button size={size} icon={<Icon name="layers" />}>
                  With icon
                </Button>
                <Button size={size} icon={<Icon name="layers" />} title="Icon only" />
                <Button size={size} variant="ghost" icon={<Icon name="layers" />} title="Icon only, ghost" />
                <Button size={size}>Text only</Button>
              </div>
            ))}
          </div>
        </div>
      </div>
    )),
};
