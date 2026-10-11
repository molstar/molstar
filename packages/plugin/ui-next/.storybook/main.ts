/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fileURLToPath } from 'node:url';
import type { StorybookConfig } from '@storybook/web-components-vite';
import { defaultClientConditions, mergeConfig } from 'vite';

const monorepoRoot = fileURLToPath(new URL('../../../..', import.meta.url));

const config: StorybookConfig = {
  stories: ['../src/**/*.mdx', '../src/**/*.stories.@(ts|tsx)'],
  addons: ['@storybook/addon-docs', '@storybook/addon-a11y'],
  framework: '@storybook/web-components-vite',
  viteFinal: (viteConfig) =>
    mergeConfig(viteConfig, {
      resolve: {
        // Resolve workspace `@molstar/*` packages to source instead of built `lib/`.
        conditions: ['molstar-src', ...defaultClientConditions],
        dedupe: ['rxjs'],
      },
      server: {
        fs: { allow: [monorepoRoot] },
      },
    }),
};

export default config;
