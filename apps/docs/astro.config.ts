/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { unified } from '@astrojs/markdown-remark';
import starlight from '@astrojs/starlight';
import { defineConfig, passthroughImageService } from 'astro/config';
import fs from 'node:fs';
import rehypeKatex from 'rehype-katex';
import remarkMath from 'remark-math';
import { DEFAULT_VERSION } from './src/versions';

const base = '/docs';

// The pages the previous (MkDocs) site served at /docs/<page>/ now live under /docs/v5/. They keep working through
// redirects that always lead to v5, whatever the default version is.
function v5Redirects(): Record<string, string> {
  const dir = new URL('./src/content/docs/v5/', import.meta.url);
  const redirects: Record<string, string> = {};
  for (const file of fs.readdirSync(dir, { recursive: true, encoding: 'utf8' })) {
    if (!file.endsWith('.md')) continue;
    const page = file
      .replace(/\\/g, '/')
      .replace(/\.md$/, '')
      .replace(/(^|\/)index$/, '');
    if (page) redirects[`/${page}`] = `${base}/v5/${page}/`;
  }
  return redirects;
}

export default defineConfig({
  site: 'https://molstar.org',
  base,
  redirects: {
    ...v5Redirects(),
    '/': `${base}/${DEFAULT_VERSION}/`,
  },
  // The few images are served as they are: no native image library (sharp) in the workspace.
  image: { service: passthroughImageService() },
  // The equations need remark and rehype plugins, which the default (Sätteri) processor does not run.
  markdown: {
    processor: unified({ remarkPlugins: [remarkMath], rehypePlugins: [rehypeKatex] }),
  },
  integrations: [
    starlight({
      title: 'Mol* Documentation',
      description: 'Documentation of Mol*, a library and set of apps for molecular structures and volumes.',
      favicon: '/favicon.ico',
      social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/molstar/molstar' }],
      customCss: ['katex/dist/katex.min.css'],
      components: { ThemeSelect: './src/components/ThemeSelect.astro' },
      routeMiddleware: './src/route-middleware.ts',
      // One group per version; the route middleware shows the group of the page's version.
      sidebar: [
        {
          label: 'v6',
          items: [{ slug: 'v6', label: 'Introduction' }],
        },
        {
          label: 'v5',
          items: [
            { slug: 'v5', label: 'Installation' },
            {
              label: 'Plugin',
              items: [
                { slug: 'v5/plugin/instance', label: 'Creating Instance' },
                { slug: 'v5/plugin/examples', label: 'Examples' },
                { slug: 'v5/plugin/custom-library', label: 'Custom Library' },
                { slug: 'v5/plugin/selections', label: 'Selections' },
                { slug: 'v5/plugin/superposition', label: 'Superposition' },
                { slug: 'v5/plugin/viewer-state', label: 'Viewer State' },
                { slug: 'v5/plugin/data-state', label: 'Data State' },
                { slug: 'v5/plugin/file-formats', label: 'File Formats' },
                { slug: 'v5/plugin/cif-schemas', label: 'CIF Schemas' },
                {
                  label: 'Managers',
                  items: [{ slug: 'v5/plugin/managers/markdown-extensions', label: 'Markdown Extensions' }],
                },
                {
                  label: 'State Transforms',
                  items: [
                    { slug: 'v5/plugin/transforms/custom-trajectory', label: 'Custom Trajectory' },
                    { slug: 'v5/plugin/transforms/custom-conformation', label: 'Custom Conformation' },
                  ],
                },
              ],
            },
            {
              label: 'Data Access Tools',
              items: [
                { slug: 'v5/data-access-tools/model-server' },
                {
                  label: 'Volume Server',
                  items: [
                    { slug: 'v5/data-access-tools/volume-server', label: 'Overview' },
                    { slug: 'v5/data-access-tools/volume-server/examples', label: 'Examples' },
                    { slug: 'v5/data-access-tools/volume-server/how-it-works', label: 'How it Works' },
                    { slug: 'v5/data-access-tools/volume-server/response-data-format', label: 'Data Format' },
                  ],
                },
                { slug: 'v5/data-access-tools/plugin-state-server' },
                { slug: 'v5/data-access-tools/convert-to-bcif' },
                { slug: 'v5/data-access-tools/create-ccd-table' },
                { slug: 'v5/data-access-tools/extract-ccd-ions' },
              ],
            },
            {
              label: 'Extensions',
              items: [
                { slug: 'v5/extensions/mvs', label: 'MolViewSpec' },
                { slug: 'v5/extensions/struct-conn', label: 'wwPDB StructConn' },
                { slug: 'v5/extensions/tunnels', label: 'Tunnels' },
                { slug: 'v5/extensions/interactions', label: 'Interactions' },
              ],
            },
            {
              label: 'Misc',
              items: [
                { slug: 'v5/misc/interesting-pdb-entries', label: 'Interesting PDB entries' },
                { slug: 'v5/misc/exporting-components', label: 'Exporting component data' },
              ],
            },
          ],
        },
      ],
    }),
  ],
});
