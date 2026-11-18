/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset, type AssetManager } from '@molstar/core/util/assets';
import { SizeLimitExceeded } from '@molstar/core/util/data-source';
import type { AssetReader } from '@molstar/seqx-core/data';

// The assets port's policy wrapper (_spec/README.md "Context", Ports): the application's `AssetManager` behind the
// rules of what may be read (schemes, hosts, relative locations under a base), how large one asset may be and how long
// a read may take, adapted to the source protocol's `AssetReader`. A read the policy refuses rejects with
// `AssetRefused`; a source's opening turns that into an issue on the source, never a throw to the caller.
//
// What the policy enforces is a decision on the location, made on the URL it resolves to (never on the string it came
// in as), and four choices it passes to the manager on every read (`AssetManager.ResolveOptions`, core's
// `util/_spec/assets.md` §8; a manager honors each or refuses the read): `credentials` (`'omit'` unless the policy says
// otherwise), `redirect: 'error'` (a redirect is not followed: the location checked is the location read), the read's
// signal (the caller's, and the timeout's) and `maxBytes`, enforced while the body streams. The policy also checks the
// size once the body is read, and abandons a read the manager did not end in time, for a manager that falls short.

export interface AssetPolicy {
  /** Absolute locations may use these schemes (`https:`); none by default. */
  readonly schemes?: readonly string[];
  /** Absolute locations may reach these hosts only; any host of an allowed scheme when absent. */
  readonly hosts?: readonly string[];
  /** Relative locations are read, resolved against `base` (true by default when there is a base). */
  readonly relative?: boolean;
  /** The location relative ones resolve against: a page, a folder (`https://x/y/`). Relative reads need it. */
  readonly base?: string;
  /** A relative location may not leave the folder of `base` (`../`, an absolute path): true by default. */
  readonly confine?: boolean;
  /** The largest asset, in bytes: enforced by the manager while the body streams, and checked once it is read. */
  readonly maxBytes?: number;
  /** A read that takes longer is abandoned, ms: the manager is told by the read's signal, which stops the download. */
  readonly timeoutMs?: number;
  /**
   * Whether cookies and HTTP authentication go with a read (`fetch`'s `credentials`): `'omit'` by default, so a
   * composition that needs a login (a lab's server) says `'include'` or `'same-origin'`.
   */
  readonly credentials?: 'omit' | 'same-origin' | 'include';
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

/** What the URL parser drops before it reads a location: leading and trailing C0 controls and spaces, tabs, newlines. */
const cleaned = (uri: string) => uri.replace(/^[\u0000-\u0020]+|[\u0000-\u0020]+$/g, '').replace(/[\t\n\r]/g, '');
const HasScheme = /^[a-z][a-z0-9+.-]*:/i;

/**
 * The URL a location resolves to under the policy, or why it is refused. A location is classified by what it resolves
 * to: one that has a scheme is absolute (its scheme and host are the policy's to allow); one that has none is relative
 * and must resolve to the base's own scheme and host (so `//host`, `\\host` and `/\host`, which the URL parser reads as
 * another host, are refused), and with `confine` stay under the base's folder (`..`, `%2e%2e` and a rooted path
 * resolve out of it; an encoded separator is refused).
 */
export function checkLocation(uri: string, policy: AssetPolicy): { url: string } | { refused: string } {
  const raw = cleaned(uri);
  const relative = !HasScheme.test(raw);
  let url: URL;
  try {
    url = new URL(raw, relative ? policy.base : undefined);
  } catch {
    return { refused: relative && policy.base === undefined ? 'relative locations are not read' : 'not a location' };
  }
  if (relative) {
    if (policy.relative === false || policy.base === undefined) return { refused: 'relative locations are not read' };
    const base = new URL(policy.base);
    if (url.protocol !== base.protocol || url.host !== base.host) return { refused: 'another host than the base' };
    // a server may decode an encoded separator into the traversal the parser did not see
    if (/%2f|%5c/i.test(url.pathname)) return { refused: 'an encoded path separator' };
    if (policy.confine !== false && !url.pathname.startsWith(new URL('./', base).pathname))
      return { refused: 'outside the base folder' };
    return { url: url.href };
  }
  if (!(policy.schemes ?? []).includes(url.protocol)) return { refused: `the scheme ${url.protocol} is not allowed` };
  if (policy.hosts && !policy.hosts.includes(url.host)) return { refused: `the host ${url.host} is not allowed` };
  return { url: url.href };
}

/** The UTF-8 length of a string, without making the bytes. */
export function utf8Length(s: string): number {
  let n = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i);
    if (c < 0x80) n += 1;
    else if (c < 0x800) n += 2;
    else if (c >= 0xd800 && c <= 0xdbff && i + 1 < s.length && (s.charCodeAt(i + 1) & 0xfc00) === 0xdc00) {
      n += 4;
      i++;
    } else n += 3;
  }
  return n;
}

/** The manager's text reads under the policy, as the source protocol's reader. It sets no header of its own. */
export function policyAssets(manager: AssetManager, policy: AssetPolicy): AssetReader {
  return {
    async readText(uri, signal) {
      const checked = checkLocation(uri, policy);
      if ('refused' in checked) throw new AssetRefused(uri, checked.refused);
      if (signal.aborted) throw signal.reason ?? new Error('aborted');
      // the manager's signal: the caller's or the timeout's, whichever comes first
      const stop = new AbortController();
      const read = manager
        .resolve(Asset.Url(checked.url), 'string', {
          store: false,
          credentials: policy.credentials ?? 'omit',
          redirect: 'error',
          signal: stop.signal,
          ...(policy.maxBytes !== undefined ? { maxBytes: policy.maxBytes } : {}),
        })
        .run();
      let timer: ReturnType<typeof setTimeout> | undefined;
      let onAbort: (() => void) | undefined;
      const limits = new Promise<never>((_, reject) => {
        onAbort = () => {
          stop.abort(signal.reason);
          reject(signal.reason ?? new Error('aborted'));
        };
        signal.addEventListener('abort', onAbort, { once: true });
        if (policy.timeoutMs !== undefined)
          timer = setTimeout(() => {
            const refused = new AssetRefused(uri, 'the read took too long');
            stop.abort(refused);
            reject(refused);
          }, policy.timeoutMs);
      });
      let wrapper: Awaited<typeof read>;
      try {
        wrapper = await Promise.race([read, limits]);
      } catch (e) {
        // a read that lost the race still finishes (a manager that ignores the signal): its result is let go
        read.then(
          (w) => w.dispose(),
          () => {},
        );
        if (e instanceof SizeLimitExceeded) throw new AssetRefused(uri, 'larger than the limit');
        throw e;
      } finally {
        if (timer !== undefined) clearTimeout(timer);
        if (onAbort) signal.removeEventListener('abort', onAbort);
      }
      try {
        const text = wrapper.data as string;
        if (policy.maxBytes !== undefined && utf8Length(text) > policy.maxBytes)
          throw new AssetRefused(uri, 'larger than the limit');
        return text;
      } finally {
        wrapper.dispose();
      }
    },
  };
}
