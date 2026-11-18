/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Links (_spec/README.md "Addressing", Links): an action and an address in Mol*'s link grammar, `!key=value&…`, each
// key and value percent-encoded as a whole (spaces, brackets, quotes, `#`, `&` and `=` too), so an address survives as
// one link destination. Keys of this viewer's prefix (`seq-` by default) are its actions; every other key is kept for
// the viewer that knows it. A link is data: parsed into actions and addresses, never into markup.

import { type Issue, issueOf, parseAddress } from '@molstar/core/repr/tree';
import { AddressIssues } from './address-issues.js';

/** The actions a link may carry: what a command on an address does in a view. */
export const LinkActions = ['focus', 'reveal', 'select', 'highlight', 'hover'] as const;
export type LinkActionName = (typeof LinkActions)[number];

export interface LinkAction {
  readonly action: LinkActionName;
  /** The address, as written (it parses). */
  readonly address: string;
}

export interface Link {
  readonly actions: readonly LinkAction[];
  /** The highlight channel the link names (`<prefix>channel`), if any. */
  readonly channel?: string;
  /** Keys of other prefixes, decoded, in order: another viewer's (Mol*'s own keys). */
  readonly other: readonly (readonly [key: string, value: string])[];
}

export interface LinkOptions {
  /** The prefix of this viewer's keys (`seq-`). */
  readonly prefix: string;
  /** The longest link read, in characters. */
  readonly maxLength: number;
  /** The most actions one link carries. */
  readonly maxActions: number;
}

export const DefaultLinkOptions: LinkOptions = { prefix: 'seq-', maxLength: 8192, maxActions: 16 };

const isAction = (a: string): a is LinkActionName => (LinkActions as readonly string[]).includes(a);

/** A link's text: `!` then each entry, key and value percent-encoded. Its addresses are printed as given. */
export function printLink(link: Partial<Link> & Pick<Link, 'actions'>, options: Partial<LinkOptions> = {}): string {
  const prefix = options.prefix ?? DefaultLinkOptions.prefix;
  const entries: (readonly [string, string])[] = [
    ...(link.other ?? []),
    ...link.actions.map((a): [string, string] => [`${prefix}${a.action}`, a.address]),
  ];
  if (link.channel !== undefined) entries.push([`${prefix}channel`, link.channel]);
  return `!${entries.map(([k, v]) => `${encodeURIComponent(k)}=${encodeURIComponent(v)}`).join('&')}`;
}

/**
 * A link's actions, or undefined when the text is not a link (no leading `!`). Entries split on `&` and the first `=`
 * before decoding, so an encoded `&` or `=` in an address stays in it. An action of the prefix this viewer lacks is
 * left out with a warning; an address that does not parse fails the link (nothing runs from a half-read link).
 */
export function parseLink(
  text: string,
  options: Partial<LinkOptions> = {},
): { readonly link?: Link; readonly issues: readonly Issue[] } | undefined {
  if (!text.startsWith('!')) return undefined;
  const { prefix, maxLength, maxActions } = { ...DefaultLinkOptions, ...options };
  const fail = (code: 'seq.address.link.syntax' | 'seq.address.link.limit', message: string, offset?: number) => ({
    issues: [issueOf(AddressIssues, code, message, offset === undefined ? {} : { offset })],
  });
  if (text.length > maxLength) return fail('seq.address.link.limit', `A link is longer than ${maxLength} characters`);
  const actions: LinkAction[] = [];
  const other: [string, string][] = [];
  const issues: Issue[] = [];
  let channel: string | undefined;
  let offset = 1;
  for (const entry of text.slice(1).split('&')) {
    const at = offset;
    offset += entry.length + 1;
    if (!entry.trim()) continue;
    const eq = entry.indexOf('=');
    let key: string;
    let value: string;
    try {
      key = decodeURIComponent(eq < 0 ? entry : entry.slice(0, eq)).trim();
      value = eq < 0 ? '' : decodeURIComponent(entry.slice(eq + 1)).trim();
    } catch {
      return fail('seq.address.link.syntax', 'An entry is not percent-encoded', at);
    }
    if (!key.startsWith(prefix)) {
      other.push([key, value]);
      continue;
    }
    const name = key.slice(prefix.length);
    if (name === 'channel') {
      channel = value;
      continue;
    }
    if (!isAction(name)) {
      issues.push(
        issueOf(AddressIssues, 'seq.address.link.action', 'A link names an action this viewer does not have', {
          offset: at,
          hint: AddressIssues.codes['seq.address.link.action'].hint,
        }),
      );
      continue;
    }
    const parsed = parseAddress(value);
    if (!parsed.pattern) return { issues: parsed.issues.map((i) => ({ ...i, offset: at })) };
    if (actions.length >= maxActions)
      return fail('seq.address.link.limit', `A link carries more than ${maxActions} actions`, at);
    actions.push({ action: name, address: value });
  }
  return { link: { actions, other, ...(channel !== undefined ? { channel } : {}) }, issues };
}
