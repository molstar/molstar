/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Preview } from '@storybook/web-components-vite';
// The decorators of the registered packages' preview modules (packages.ts).
import packageDecorators from 'virtual:story-package-previews';
import { ThemedDocsContainer } from './docs-container.js';

const preview: Preview = {
  parameters: {
    layout: 'fullscreen',
    docs: { container: ThemedDocsContainer },
    options: {
      // One group per registered package group (packages.ts). Prism: Overview first, then the real viewer from the
      // whole to its parts and data, then prism-ui's components over stand-in models, the runtime last.
      storySort: {
        order: [
          'Prism',
          [
            'Overview',
            ['Introduction', 'Foundations'],
            'Viewer',
            [
              'Shell',
              'Home',
              'Structure panel',
              'Volumes',
              'Playback',
              'Assistant',
              'Data',
              'MVS stories',
              'With extensions',
            ],
            'Components',
            [
              'Button',
              'Icon',
              'Logo',
              'Text input',
              'Select',
              'Choices',
              'Params form',
              'Legend',
              'Markdown',
              'Tooltip',
              'Scroll area',
              'Drag to reorder',
              'Action footer',
              'Playback',
            ],
            'Internals',
            ['Runtime', 'Runtime performance', 'State tree performance', 'React and Preact'],
          ],
          // Seq*: the docs pages first (in reading order), then the cases, then the topics.
          'Seq*',
          [
            'Docs',
            [
              'Introduction',
              'Concepts',
              'Principles',
              'Features',
              'Format',
              'Query language',
              'Streaming',
              'Color schemes',
              'Extending',
            ],
            'Cases',
            ['Viewer', 'Format', 'Mol*', 'Nightingale', 'Genome browser', 'Default UIs'],
            'Views and layouts',
            'Themes',
            'Interaction',
            'Query and search',
            'Streaming',
            'File formats',
            'Performance',
          ],
        ],
      },
    },
  },
  initialGlobals: {
    theme: 'light',
  },
  globalTypes: {
    theme: {
      description: 'UI theme (sets `data-msp-theme` on the document root)',
      toolbar: {
        title: 'Theme',
        icon: 'paintbrush',
        items: [
          { value: 'light', title: 'Light' },
          { value: 'dark', title: 'Dark' },
          { value: 'blue', title: 'Blue' },
          { value: 'auto', title: 'Auto (system)' },
        ],
        dynamicTitle: true,
      },
    },
  },
  decorators: [
    // First, so innermost: a package's decorator may render the story itself (React) instead of calling it, which would
    // skip the decorators inside it, so the theme sits outside.
    ...packageDecorators,
    (story, context) => {
      document.documentElement.dataset.mspTheme = context.globals.theme ?? 'light';
      document.body.style.background = 'var(--msp-color-bg)';
      return story();
    },
  ],
};

export default preview;
