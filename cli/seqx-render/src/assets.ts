/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { realpath } from 'node:fs/promises';
import { sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { AssetRefused, checkLocation, headlessAssetPolicy, policyAssets } from '@molstar/seqx-core/context';
import type { AssetReader } from '@molstar/seqx-core/data';

// The CLI's assets under the headless policy, with two more rules the policy cannot hold on a location alone: a file
// is read only when its real path (links resolved) is under the real path of the root, so a link in the spec's folder
// that leads out of it is refused; and each refusal is kept with its reason, for the CLI to print (the source's issue
// does not quote a location or its cause).

export interface CliAssets {
  readonly reader: AssetReader;
  /** The refusals so far, each as the policy words it (`The location is not read: <uri> (<reason>)`). */
  readonly refusals: ReadonlySet<string>;
}

/** Reads under `root` and from the `https:` `hosts` named, through `manager` (a `DefaultAssetManager`). */
export function cliAssets(root: string, hosts: readonly string[], manager = new DefaultAssetManager()): CliAssets {
  const policy = headlessAssetPolicy({ root, hosts });
  const inner = policyAssets(manager, policy);
  const refusals = new Set<string>();
  let rootReal: Promise<string> | undefined;

  async function guard(uri: string) {
    const checked = checkLocation(uri, policy);
    // a refused location is the inner reader's to refuse, with its reason
    if ('refused' in checked || !checked.url.startsWith('file:')) return;
    rootReal ??= realpath(root);
    let real: string;
    try {
      real = await realpath(fileURLToPath(checked.url));
    } catch {
      return; // not there: the read says so
    }
    const base = await rootReal;
    if (real !== base && !real.startsWith(base.endsWith(sep) ? base : base + sep))
      throw new AssetRefused(uri, 'a link that leads outside the root');
  }

  async function guarded<T>(uri: string, read: () => Promise<T>): Promise<T> {
    try {
      await guard(uri);
      return await read();
    } catch (e) {
      if (e instanceof AssetRefused) refusals.add(e.message);
      throw e;
    }
  }

  return {
    refusals,
    reader: {
      readText: (uri, signal) => guarded(uri, () => inner.readText(uri, signal)),
      readBytes: (uri, signal) => guarded(uri, () => inner.readBytes(uri, signal)),
    },
  };
}
