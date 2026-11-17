/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset, type AssetManager } from '@molstar/core/util/assets';
import type { AssetReader } from '@molstar/seqx-core/data';

// The assets port's policy wrapper (_spec/README.md "Context", Ports): the application's `AssetManager` behind the
// rules of what may be read (schemes, hosts, relative locations under a base), how large one asset may be, how long
// a read may take, and no credentials, adapted to the source protocol's `AssetReader`. A read the policy refuses
// rejects with `AssetRefused`; a source's opening turns that into an issue on the source, never a throw to the caller.

export interface AssetPolicy {
  /** Absolute locations may use these schemes (`https:`); none by default. */
  readonly schemes?: readonly string[];
  /** Absolute locations may reach these hosts only; any host of an allowed scheme when absent. */
  readonly hosts?: readonly string[];
  /** Relative locations are read, resolved against `base` (true by default when there is a base). */
  readonly relative?: boolean;
  /** The location relative ones resolve against: a page, a folder (`https://x/y/`). Relative reads need it. */
  readonly base?: string;
  /** A relative location may not leave the folder of `base` (`../`): true by default. */
  readonly confine?: boolean;
  /** The largest asset, in characters of text (a manager reads the whole body before it is checked). */
  readonly maxBytes?: number;
  /** A read that takes longer is abandoned, ms. */
  readonly timeoutMs?: number;
}

/** The browser's: `https:` and relative locations under the page, 256 MB an asset, 30 s a read. */
export const BrowserAssetPolicy: AssetPolicy = {
  schemes: ['https:'],
  relative: true,
  confine: true,
  maxBytes: 256 * 1024 * 1024,
  timeoutMs: 30_000,
};

/** The headless default: nothing may be read (a caller names what may). */
export const HeadlessAssetPolicy: AssetPolicy = { schemes: [], relative: false };

/** A location the policy refuses, or a read it abandoned. */
export class AssetRefused extends Error {
  constructor(
    readonly uri: string,
    why: string,
  ) {
    super(`The location is not read: ${uri} (${why})`);
    this.name = 'AssetRefused';
  }
}

const HasScheme = /^[a-z][a-z0-9+.-]*:/i;

/** The URL a location resolves to under the policy, or why it is refused. */
export function checkLocation(uri: string, policy: AssetPolicy): { url: string } | { refused: string } {
  const relative = !HasScheme.test(uri) && !uri.startsWith('//');
  if (relative) {
    if (policy.relative === false || policy.base === undefined) return { refused: 'relative locations are not read' };
    let url: URL;
    try {
      url = new URL(uri, policy.base);
    } catch {
      return { refused: 'not a location' };
    }
    if (policy.confine !== false && !url.href.startsWith(new URL('./', policy.base).href))
      return { refused: 'outside the base folder' };
    return { url: url.href };
  }
  let url: URL;
  try {
    url = new URL(uri, policy.base);
  } catch {
    return { refused: 'not a location' };
  }
  if (!(policy.schemes ?? []).includes(url.protocol)) return { refused: `the scheme ${url.protocol} is not allowed` };
  if (policy.hosts && !policy.hosts.includes(url.host)) return { refused: `the host ${url.host} is not allowed` };
  return { url: url.href };
}

/** The manager's text reads under the policy, as the source protocol's reader. Never sends credentials or headers. */
export function policyAssets(manager: AssetManager, policy: AssetPolicy): AssetReader {
  return {
    async readText(uri, signal) {
      const checked = checkLocation(uri, policy);
      if ('refused' in checked) throw new AssetRefused(uri, checked.refused);
      if (signal.aborted) throw signal.reason ?? new Error('aborted');
      const read = manager.resolve(Asset.Url(checked.url), 'string', { store: false }).run();
      let timer: ReturnType<typeof setTimeout> | undefined;
      let onAbort: (() => void) | undefined;
      const limits = new Promise<never>((_, reject) => {
        onAbort = () => reject(signal.reason ?? new Error('aborted'));
        signal.addEventListener('abort', onAbort, { once: true });
        if (policy.timeoutMs !== undefined)
          timer = setTimeout(() => reject(new AssetRefused(uri, 'the read took too long')), policy.timeoutMs);
      });
      let wrapper: Awaited<typeof read>;
      try {
        wrapper = await Promise.race([read, limits]);
      } catch (e) {
        // a read that lost the race still finishes: its result is let go
        read.then(
          (w) => w.dispose(),
          () => {},
        );
        throw e;
      } finally {
        if (timer !== undefined) clearTimeout(timer);
        if (onAbort) signal.removeEventListener('abort', onAbort);
      }
      try {
        const text = wrapper.data as string;
        if (policy.maxBytes !== undefined && text.length > policy.maxBytes)
          throw new AssetRefused(uri, 'larger than the limit');
        return text;
      } finally {
        wrapper.dispose();
      }
    },
  };
}
