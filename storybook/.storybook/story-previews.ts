/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';
import type { Plugin } from 'vite';
import { packageRoot, type StoryPackage } from '../packages.js';

const ModuleId = 'virtual:story-package-previews';
const ResolvedId = `\0${ModuleId}`;

/**
 * `virtual:story-package-previews`: the `decorators` of the packages' preview modules (`StoryPackage.preview`), in
 * registration order, for preview.ts to apply. A package that needs none has no entry.
 */
export function storyPreviewsPlugin(packages: StoryPackage[]): Plugin {
  const previews = packages.flatMap((pkg) => (pkg.preview ? [path.join(packageRoot(pkg), pkg.preview)] : []));
  return {
    name: 'msp-story-package-previews',
    resolveId: (id) => (id === ModuleId ? ResolvedId : null),
    load(id) {
      if (id !== ResolvedId) return null;
      const imports = previews.map((file, i) => `import { decorators as d${i} } from ${JSON.stringify(file)};`);
      return `${imports.join('\n')}\nexport default [${previews.map((_, i) => `...d${i}`).join(', ')}];\n`;
    },
  };
}
