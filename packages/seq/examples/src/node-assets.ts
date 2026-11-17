/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFile } from 'node:fs/promises';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { AssetRefused, checkLocation, type ReaderAssetsPort } from '@molstar/seqx-core/context/index';

// The assets port's Node adapter for the examples (and a second adapter of the port, for its suite): relative locations
// are read from the package folder (`fixtures/p53/tp53-gene.seq`), nothing else: no network, no absolute path, nothing
// above the folder. A location it refuses is the same refusal the browser's policy gives, which a source's opening turns
// into an issue.

/** The package folder: where `fixtures/` is, from the sources and from the built library alike. */
export const PackageFolder = fileURLToPath(new URL('..', import.meta.url));

/** The adapter: a port to give a context (`ports: { assets: () => nodeFixtureAssets() }`). */
export function nodeFixtureAssets(root: string = PackageFolder): ReaderAssetsPort {
  const base = pathToFileURL(root.endsWith('/') ? root : `${root}/`).href;
  const policy = { schemes: [], relative: true, confine: true, base } as const;
  return {
    reader: {
      async readText(uri, signal) {
        const checked = checkLocation(uri, policy);
        if ('refused' in checked) throw new AssetRefused(uri, checked.refused);
        return readFile(fileURLToPath(checked.url), { encoding: 'utf8', signal });
      },
    },
  };
}
