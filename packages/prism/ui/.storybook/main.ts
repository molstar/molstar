/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fileURLToPath } from 'node:url';
import type { StorybookConfig } from '@storybook/web-components-vite';
import { defaultClientConditions, mergeConfig, type Plugin } from 'vite';
import { cssModulesPlugin, cssPackage } from '../scripts/css.mjs';

const monorepoRoot = fileURLToPath(new URL('../../../..', import.meta.url));
// Packages built on prism-ui show their stories here too (the viewer mode: Shell › Viewer mode, Views; the assistant: Shell › Assistant; story mode: Shell › Story mode).
const assistantRoot = fileURLToPath(new URL('../../assistant', import.meta.url));
const storiesRoot = fileURLToPath(new URL('../../stories', import.meta.url));
const viewerRoot = fileURLToPath(new URL('../../viewer', import.meta.url));

// Internals › React and Preact runs one adapter on both: the modules in src/react/_stories/preact/, and the adapter
// modules they import, resolve `react` to `preact/compat` (as a Preact app's bundler does for the whole app). The adapter
// modules are loaded a second time for them (`?preact`), so the React story keeps real React.
const preactFolder = fileURLToPath(new URL('../src/react/_stories/preact/', import.meta.url));
const adapterFolder = fileURLToPath(new URL('../src/react/', import.meta.url));
const PreactCopy = '?preact';
const ReactToPreact: Record<string, string> = {
  react: 'preact/compat',
  'react-dom': 'preact/compat',
  'react-dom/client': 'preact/compat/client',
  'react/jsx-runtime': 'preact/jsx-runtime',
};
const BareReact = /(\bfrom\s*|\bimport\s*\(?\s*)(['"])(react|react-dom|react-dom\/client|react\/jsx-runtime)\2/g;

function preactFolderPlugin(): Plugin {
  return {
    name: 'msp-preact-folder',
    enforce: 'pre',
    // The adapter modules this folder imports get a copy of their own (`?preact`).
    async resolveId(source, importer, options) {
      if (!importer || !source.startsWith('.')) return null;
      const file = importer.split('?')[0];
      if (!importer.endsWith(PreactCopy) && !file.startsWith(preactFolder)) return null;
      const resolved = await this.resolve(source, file, { ...options, skipSelf: true });
      if (!resolved || resolved.external) return resolved;
      const id = resolved.id.split('?')[0];
      const adapter = id.startsWith(adapterFolder) && !id.startsWith(preactFolder) && !id.includes('/_stories/');
      return adapter ? { ...resolved, id: `${id}${PreactCopy}` } : resolved;
    },
    // In those copies, `react` is `preact/compat` (rewritten in the source: Vite resolves bare imports to its
    // pre-bundled dependencies before plugins see them).
    transform(code, id) {
      if (!id.endsWith(PreactCopy)) return null;
      return {
        code: code.replace(BareReact, (_, lead, q, name) => `${lead}${q}${ReactToPreact[name]}${q}`),
        map: null,
      };
    },
  };
}

const config: StorybookConfig = {
  stories: [
    '../src/**/*.mdx',
    '../src/**/*.stories.@(ts|tsx)',
    '../../assistant/src/**/*.stories.@(ts|tsx)',
    '../../stories/src/**/*.stories.@(ts|tsx)',
    '../../viewer/src/**/*.stories.@(ts|tsx)',
  ],
  addons: ['@storybook/addon-docs', '@storybook/addon-a11y'],
  framework: '@storybook/web-components-vite',
  viteFinal: (viteConfig) =>
    mergeConfig(viteConfig, {
      // Keeps the generated `*.css.ts` modules current while Storybook runs.
      plugins: [
        cssModulesPlugin([
          cssPackage(),
          cssPackage(assistantRoot, 'assistant'),
          cssPackage(storiesRoot, 'stories'),
          cssPackage(viewerRoot, 'viewer'),
        ]),
        preactFolderPlugin(),
      ],
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
