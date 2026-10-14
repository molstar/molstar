/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';
import { fileURLToPath } from 'node:url';
import type { PluginOption } from 'vite';
// The Prism CSS generator, run as a Vite plugin (the Prism packages' `styles` scripts run the same file).
import { cssModulesPlugin, cssPackage } from '../packages/prism/ui/scripts/css.mjs';
import { preactFolderPlugin } from './.storybook/preact-folder.js';
import { nightingaleSideEffects, seqstarFixtureFs } from '../packages/sequences/storybook/vite.js';

export const monorepoRoot = fileURLToPath(new URL('..', import.meta.url));

/** A package whose stories Storybook shows. */
export interface StoryPackage {
  /** The package folder, relative to the repository root. Every `*.stories.ts(x)` and `*.mdx` under its `src/` is shown. */
  dir: string;
  /**
   * The sidebar group its stories go under (Storybook's `titlePrefix`), e.g. `Prism` for `Prism/Shell/…`; none when the
   * story titles already start with their group.
   */
  group?: string;
  /** Vite plugins its stories need (given the package's absolute path). */
  plugins?: (root: string) => PluginOption[];
  /**
   * A module of the package, relative to its folder, that exports `decorators`, the ones its stories need (and which
   * may import the global styles they rely on); it joins the preview. E.g. a decorator that mounts a story the web
   * components renderer cannot show.
   */
  preview?: string;
}

/** The packages Storybook collects stories from; register a package here to show its stories. */
export const StoryPackages: StoryPackage[] = [
  {
    dir: 'packages/prism/ui',
    group: 'Prism',
    plugins: (root) => [cssModulesPlugin([cssPackage(root)]), preactFolderPlugin(root)],
  },
  // Built on prism-ui (the viewer mode: Prism › Shell › Viewer mode, Views; story mode: Prism › Shell › Story mode; the
  // assistant: Prism › Shell › Assistant). Their CSS class names carry the package prefix.
  {
    dir: 'packages/prism/viewer',
    group: 'Prism',
    plugins: (root) => [cssModulesPlugin([cssPackage(root, 'viewer')])],
  },
  {
    dir: 'packages/prism/stories',
    group: 'Prism',
    plugins: (root) => [cssModulesPlugin([cssPackage(root, 'stories')])],
  },
  {
    dir: 'packages/prism/assistant',
    group: 'Prism',
    plugins: (root) => [cssModulesPlugin([cssPackage(root, 'assistant')])],
  },
  // Seq*, the sequence viewer: React stories and MDX pages (`Seq*/…` titles). Its preview module mounts the stories with
  // `react-dom/client` and brings Prism's tokens; its Vite plugins serve the test fixtures the stories compile cases
  // from, and keep the Nightingale elements in a build.
  {
    dir: 'packages/sequences/storybook',
    plugins: () => [seqstarFixtureFs(), nightingaleSideEffects()],
    preview: 'src/preview.ts',
  },
];

export function packageRoot(pkg: StoryPackage) {
  return path.join(monorepoRoot, pkg.dir);
}
