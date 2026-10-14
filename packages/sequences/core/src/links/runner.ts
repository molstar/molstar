/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The runner (seqspec:libraries#links): a renderer's link component calls it on pointer enter, leave and click. It runs
// the invocations of the keys it knows through the function it is given (`ctx.commands.run(inv, "host:links")`), and
// remembers the links entered so one removed while hovered, or a disposed runner, leaves nothing lit.

import type { Invocation } from "@molstar/seq-core/format";
import { type ParsedLink, parseLink } from "./codec.js";
import { DefaultLinkKeys, type LinkKey } from "./keys.js";

export interface LinkRunner {
  /** False: no key of the href has an entry (not a command, or only other viewers' keys). */
  run(href: string, event: "enter" | "leave" | "click", link: string): boolean;
  /** The link left the page: its pending leave runs. */
  release(link: string): void;
  /** Releases every entered link. */
  dispose(): void;
}

export function createLinkRunner(
  run: (invocation: Invocation) => unknown,
  opts: { keys?: readonly LinkKey[]; cacheLinks?: number } = {},
): LinkRunner {
  const keys = new Map((opts.keys ?? DefaultLinkKeys).map((k) => [k.key, k]));
  const limit = opts.cacheLinks ?? 256;
  /** Parsed hrefs, oldest first (a Map keeps insertion order). */
  const cache = new Map<string, readonly ParsedLink[]>();
  const entered = new Map<string, string>();
  const parsed = (href: string) => {
    let p = cache.get(href);
    if (p) return p;
    p = parseLink(href);
    cache.set(href, p);
    for (const k of cache.keys()) {
      if (cache.size <= limit) break;
      cache.delete(k);
    }
    return p;
  };
  const go = (href: string, event: "enter" | "leave" | "click", link: string) => {
    let known = false;
    let ran = false;
    for (const p of parsed(href)) {
      const entry = keys.get(p.key);
      if (!entry) continue;
      known = true;
      const inv = entry[event]?.({ ...p, link });
      if (!inv) continue;
      ran = true;
      // A rejection is reported where the command failed; a link has nobody to tell.
      try {
        Promise.resolve(run(inv)).catch(() => {});
      } catch {}
    }
    if (event === "enter" && ran) entered.set(link, href);
    if (event === "leave") entered.delete(link);
    return known;
  };
  return {
    run: go,
    release(link) {
      const href = entered.get(link);
      if (href !== undefined) go(href, "leave", link);
    },
    dispose() {
      for (const [link, href] of [...entered]) go(href, "leave", link);
      cache.clear();
    },
  };
}
