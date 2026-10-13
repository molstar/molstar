/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import path from 'node:path';
import type { Plugin } from 'vite';

// Internals › React and Preact runs one adapter on both: the modules in prism-ui's src/react/_stories/preact/, and the
// adapter modules they import, resolve `react` to `preact/compat` (as a Preact app's bundler does for the whole app).
// The adapter modules are loaded a second time for them (`?preact`), so the React story keeps real React.
const PreactCopy = '?preact';
const ReactToPreact: Record<string, string> = {
  react: 'preact/compat',
  'react-dom': 'preact/compat',
  'react-dom/client': 'preact/compat/client',
  'react/jsx-runtime': 'preact/jsx-runtime',
};
const BareReact = /(\bfrom\s*|\bimport\s*\(?\s*)(['"])(react|react-dom|react-dom\/client|react\/jsx-runtime)\2/g;

/** The plugin for prism-ui (`uiRoot`, its absolute path). */
export function preactFolderPlugin(uiRoot: string): Plugin {
  const adapterFolder = path.join(uiRoot, 'src/react') + path.sep;
  const preactFolder = path.join(adapterFolder, '_stories/preact') + path.sep;
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
