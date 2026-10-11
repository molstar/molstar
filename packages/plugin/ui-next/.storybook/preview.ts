/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Preview } from '@storybook/web-components-vite';

const preview: Preview = {
  parameters: {
    layout: 'fullscreen',
  },
  initialGlobals: {
    theme: 'light',
  },
  globalTypes: {
    theme: {
      description: 'Color theme (sets `data-theme` on the document root)',
      toolbar: {
        title: 'Theme',
        icon: 'paintbrush',
        items: [
          { value: 'light', title: 'Light' },
          { value: 'dark', title: 'Dark' },
        ],
        dynamicTitle: true,
      },
    },
  },
  decorators: [
    (story, context) => {
      document.documentElement.dataset.theme = context.globals.theme ?? 'light';
      return story();
    },
  ],
};

export default preview;
