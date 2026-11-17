/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFile, realpath } from 'node:fs/promises';
import { sep } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
import { AssetRefused, checkLocation, type ReaderAssetsPort } from '@molstar/seqx-core/context';

// The assets port's Node adapter for the examples (and a second adapter of the port, for its suite): a location is read
// relative to the package folder (`fixtures/p53/tp53-gene.seq`, the way a spec names it) and only when it is a file
// under the package's `fixtures/` folder: no network, no absolute path, nothing above the folder, and no link that
// leads out of it (the real path of the file must be under the real path of the folder). A refusal is the same one the
// browser's policy gives, which a source's opening turns into an issue.

/** The package folder: where `fixtures/` is, from the sources and from the built library alike. */
export const PackageFolder = fileURLToPath(new URL('..', import.meta.url));

/** The adapter: a port to give a context (`ports: { assets: () => nodeFixtureAssets() }`). */
export function nodeFixtureAssets(packageFolder: string = PackageFolder): ReaderAssetsPort {
  const folder = packageFolder.endsWith('/') ? packageFolder : `${packageFolder}/`;
  const base = pathToFileURL(folder).href;
  const fixtures = pathToFileURL(`${folder}fixtures/`).pathname;
  const policy = { schemes: [], relative: true, confine: true, base } as const;
  return {
    reader: {
      async readText(uri, signal) {
        const checked = checkLocation(uri, policy);
        if ('refused' in checked) throw new AssetRefused(uri, checked.refused);
        const url = new URL(checked.url);
        if (!url.pathname.startsWith(fixtures)) throw new AssetRefused(uri, 'outside the fixtures folder');
        let file: string;
        let root: string;
        try {
          file = await realpath(fileURLToPath(url));
          root = await realpath(fileURLToPath(`file://${fixtures}`));
        } catch {
          throw new AssetRefused(uri, 'not a file of the fixtures folder');
        }
        // the real path: a link in the folder may lead anywhere
        if (!file.startsWith(root.endsWith(sep) ? root : `${root}${sep}`))
          throw new AssetRefused(uri, 'outside the fixtures folder');
        return readFile(file, { encoding: 'utf8', signal });
      },
    },
  };
}
