/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fileURLToPath } from 'node:url';
import type { StorybookConfig } from '@storybook/web-components-vite';
import { defaultClientConditions, mergeConfig } from 'vite';
import { cssModulesPlugin, cssPackage } from '../scripts/css.mjs';

const monorepoRoot = fileURLToPath(new URL('../../../..', import.meta.url));
// Packages built on prism-ui show their stories here too (the agent module: Shell › Assistant).
const agentRoot = fileURLToPath(new URL('../../agent', import.meta.url));

const config: StorybookConfig = {
  stories: ['../src/**/*.mdx', '../src/**/*.stories.@(ts|tsx)', '../../agent/src/**/*.stories.@(ts|tsx)'],
  addons: ['@storybook/addon-docs', '@storybook/addon-a11y'],
  framework: '@storybook/web-components-vite',
  viteFinal: (viteConfig) =>
    mergeConfig(viteConfig, {
      // Keeps the generated `*.css.ts` modules current while Storybook runs.
      plugins: [cssModulesPlugin([cssPackage(), cssPackage(agentRoot, 'agent')])],
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
