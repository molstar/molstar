/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The addressing codes (_spec/README.md "Addressing", Issues): what an address names once its refs are followed, and
// the links that carry addresses. An issue's `address` is the address asked; `places` name what it led to.

import { defineIssues } from '@molstar/core/repr/tree';

export const AddressIssues = defineIssues('seq.address', {
  'seq.address.missing': {
    severity: 'warning',
    description:
      'An address names nothing: no node of the spec at its ids, or a ref whose target is gone (reported with what it named)',
    hint: 'Point it at what is in the spec now; an address is never re-bound to something else',
  },
  'seq.address.waiting': {
    severity: 'info',
    description: 'An address names a node data makes that is not held yet: it resolves when the data arrives',
  },
  'seq.address.ref.cycle': {
    severity: 'error',
    description: 'A ref leads back to itself through other refs',
    hint: 'Point one of the refs at what it names, not at another ref',
  },
  'seq.address.ref.depth': {
    severity: 'error',
    description: 'A chain of refs is longer than the bound (`address.refDepth`)',
    hint: 'Point the ref at what it names, or raise the bound',
  },
  'seq.address.ref.set': {
    severity: 'error',
    description: 'A ref to a set (or a subtree) is followed by more segments, which would apply to each member',
    hint: 'Write the segments on each member of the set',
  },
  'seq.address.ref.empty': {
    severity: 'error',
    description: 'A ref names no address',
  },
  'seq.address.type': {
    severity: 'error',
    description: 'A set of addresses names results of different types (nodes, elements, a location)',
    hint: 'Keep one type of result in a set',
  },
  'seq.address.link.syntax': {
    severity: 'error',
    description: 'A link is not `!key=value&…` with percent-encoded keys and values',
  },
  'seq.address.link.action': {
    severity: 'warning',
    description: 'A link names an action of this viewer’s prefix that the viewer does not have; it is left out',
    hint: 'focus, reveal, select, highlight or hover',
  },
  'seq.address.link.limit': {
    severity: 'error',
    description: 'A link is longer than the bound (`links.maxLength`), or carries more actions than it allows',
  },
});
