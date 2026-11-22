/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Params the language declares and this version reads nowhere yet (_spec/README.md "Resolve", Unread params): one
// info issue each, naming the lane that will read it, so a spec using one is never drawn without it in silence.

import { type Issue, issueOf, type Node } from '@molstar/core/repr/tree';
import { ResolveIssues } from './issues.js';

/** A kind's unread params, each followed by the lane that reads it (`*`: the node itself). */
const Unread: Readonly<Record<string, string>> = {
  map: 'source FR',
  assembly: 'names GF',
};

/** The unread params of the spec into `issues`. */
export function unreadParams(spec: Node, issues: Issue[]): void {
  const say = (what: string, lane: string, address: string) =>
    issues.push(issueOf(ResolveIssues, 'seq.resolve.param.unread', `${what}: read by lane ${lane}`, { address }));
  const walk = (n: Node, at: string) => {
    const w = Object.hasOwn(Unread, n.kind) ? Unread[n.kind].split(' ') : [];
    for (let i = 0; i < w.length; i += 2)
      if (w[i] === '*') say(`A ${n.kind}`, w[i + 1], at);
      else if (n.params?.[w[i]] !== undefined) say(`\`${n.kind}.${w[i]}\``, w[i + 1], at);
    for (const c of n.children ?? []) if (c.id !== undefined) walk(c, `${at}/${c.id}`);
  };
  for (const c of spec.children ?? []) if (c.id !== undefined) walk(c, `/${c.id}`);
}
