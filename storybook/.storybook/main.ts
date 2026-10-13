/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { StorybookConfig } from '@storybook/web-components-vite';
import { defaultClientConditions, mergeConfig } from 'vite';
import { monorepoRoot, packageRoot, StoryPackages } from '../packages.js';

// The stories of every registered package (../packages.ts), under its group; by convention they sit in `_stories/`
// folders.
const config: StorybookConfig = {
  stories: StoryPackages.map((pkg) => ({
    directory: `../../${pkg.dir}/src`,
    files: '**/*.@(mdx|stories.ts|stories.tsx)',
    titlePrefix: pkg.group,
  })),
  addons: ['@storybook/addon-docs', '@storybook/addon-a11y'],
  framework: '@storybook/web-components-vite',
  viteFinal: (viteConfig) =>
    mergeConfig(viteConfig, {
      plugins: StoryPackages.flatMap((pkg) => pkg.plugins?.(packageRoot(pkg)) ?? []),
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
