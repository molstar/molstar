/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link grammar (seqspec:libraries#links): `!key=value&key=value`, URI-decoded, Mol*'s syntax. One href may carry
// Mol*'s keys and Seq*'s; each side runs the keys it knows. This reads link targets and finds no links: the host's
// renderer does.

export interface ParsedLink {
  /** The caller's id for one link on the page (empty from `parseLink`). */
  readonly link: string;
  readonly key: string;
  readonly value: string;
  /** Every key of the href, URI-decoded (range args, `section`, other viewers' keys). */
  readonly args: Readonly<Record<string, string>>;
}

/** One per key of a `!` href, in order; [] for anything else. A malformed escape skips its pair. */
export function parseLink(href: string, link = ''): readonly ParsedLink[] {
  if (!href.startsWith('!')) return [];
  const pairs: [string, string][] = [];
  const args: Record<string, string> = {};
  for (const part of href.slice(1).split('&')) {
    const eq = part.indexOf('=');
    try {
      const k = decodeURIComponent(eq < 0 ? part : part.slice(0, eq));
      const v = eq < 0 ? '' : decodeURIComponent(part.slice(eq + 1));
      pairs.push([k, v]);
      args[k] = v;
    } catch {}
  }
  return pairs.map(([key, value]) => ({ link, key, value, args }));
}

/** Commas stay readable (a value's `refs` are comma-separated); parentheses are escaped, so the href is a markdown link
 * destination as it is (an unbalanced `)` would end it). */
const esc = (s: string) =>
  encodeURIComponent(s)
    .replace(/%2C/g, ',')
    .replace(/[()]/g, (c) => (c === '(' ? '%28' : '%29'));

/** For builders and agents: `parseLink(formatLink(k, v, a))` gives back `k`, `v` and `a`. */
export function formatLink(key: string, value: string, args?: Readonly<Record<string, string | number>>): string {
  let out = `!${esc(key)}=${esc(value)}`;
  for (const [k, v] of Object.entries(args ?? {})) out += `&${esc(k)}=${esc(String(v))}`;
  return out;
}
