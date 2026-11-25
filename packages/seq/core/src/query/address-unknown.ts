/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A node an address waits for that its source says it has none of (_spec/README.md "Addressing", Resolving: Waiting
// for nodes). Its own module: only whoever asks a source needs it (the viewer's lazy node seeking).

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import type { WaitedNode } from './address.js';
import { AddressIssues } from './address-issues.js';

/**
 * A node an address waits for whose source, asked where its key is, has none (`seq.address.node.unknown`): the address
 * asked, and the source as its place. Whoever asked the source says it.
 */
export function unknownNodeIssue(node: WaitedNode, address: string): Issue {
  const issue = issueOf(AddressIssues, 'seq.address.node.unknown', `Its source has no ${node.key}`, {
    address,
    hint: AddressIssues.codes['seq.address.node.unknown'].hint,
  });
  return { ...issue, places: [{ address: node.source }] };
}
