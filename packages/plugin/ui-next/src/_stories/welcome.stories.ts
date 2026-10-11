/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';

const meta: Meta = {
  title: 'ui-next/Welcome',
  render: (_args, context) => {
    const root = document.createElement('div');
    root.style.cssText = 'font-family: sans-serif; padding: 16px;';

    const title = document.createElement('h1');
    title.textContent = 'Mol* ui-next';
    root.appendChild(title);

    const theme = document.createElement('p');
    theme.textContent = `Theme: ${context.globals.theme ?? 'light'}`;
    root.appendChild(theme);

    return root;
  },
};

export default meta;

export const Default: StoryObj = {};
