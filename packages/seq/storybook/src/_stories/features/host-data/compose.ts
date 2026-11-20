/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { withAdapter } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { nodeAddress } from '@molstar/seqx-core/resolve';
import { MmcifRuntime } from '@molstar/seqx-core/runtimes/mmcif';
import { readChainTable } from '@molstar/seqx-examples/p53/structure-record';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { ChainTable } from '../../support/fake-structure/chain-table.js';
import type { FakeStructureViewer } from '../../support/fake-structure/fake-viewer.js';
import { fakeStructureAdapter } from '../../support/fake-structure-adapter.js';

// The Host data story's composition, DOM-free (the story and its headless twin call these): one spec naming PDB 2OCJ
// by its file (`format: mmcif`, a `uri`, no host named), composed with the fake structure viewer's adapter (the fake
// holds 2OCJ and serves it from memory) or with the default (the mmCIF runtime reading the file through the assets
// port). The same spec, the same nodes; "Drop chain" in the fake is the host's edit the sequences follow.

/** The entry's file (what the default adapter reads) and the fake's chain table of it, relative to the examples. */
export const HostFile = 'fixtures/structures/2OCJ.cif.gz';
export const HostTableFile = 'fixtures/p53/2OCJ-chains.json';
/** The source's address, and the nodes of it the spec reads. */
export const HostSource = '/data/2ocj';
export const HostComplex = nodeAddress(HostSource, 'complex');
/** The chains the view shows: the first and the last of the entry's four. */
export const HostChains = ['A', 'D'] as const;
/** A unit of the complex: its address. */
export const unitAddress = (chain: string) => `${HostComplex}/${chain}`;

/** The chain table of 2OCJ from its file's text (the story fetches it, a twin reads it from disk). */
export const hostTableOf = (text: string): ChainTable => readChainTable(text);

/**
 * The spec: the source by its file, and two chains as sections of their own, A and D: each a sequence of the spec that
 * names the source's unit (the source gives its length) with a B-factor track reading the source's `b_factor` channel
 * (`values.field`) over it. A spec's `values` node targets a system the spec holds, so a unit is named by a sequence of
 * the spec's; the source's own nodes (its sequences, molecules and complex) are the store's, listed beside the view.
 */
export function hostDataSpec(): Node {
  const b = createSeqSpec({ title: 'Host data: PDB 2OCJ' });
  const source = b.data.source('2ocj', { format: 'mmcif', uri: HostFile });
  const chains = HostChains.map((chain) => {
    const sequence = b.data.sequence(`chain-${chain}`, {
      alphabet: 'protein',
      name: chain,
      source,
    });
    return {
      chain,
      sequence,
      values: b.data.values(`bfactor-${chain}`, { target: sequence, source, field: 'b_factor' }),
    };
  });
  b.views.view('main', (v) => {
    for (const { chain, sequence, values } of chains)
      v.section(`chain-${chain}`, { of: sequence }, (s) => {
        s.track('ruler', { label: `Chain ${chain}` }).layer('ruler', { type: 'ruler' });
        s.track('bfactor', { label: 'B-factor' }).layer('bfactor', { type: 'values', data: values });
      });
  });
  return b.build();
}

/**
 * A viewer registry reading `format: mmcif` with the fake structure viewer's adapter (the mmCIF runtime its fallback:
 * the sources the fake does not hold are read from their files).
 */
export function withFakeHost<E extends { readonly kind: string }>(registry: readonly E[], fake: FakeStructureViewer) {
  return withAdapter(registry, fakeStructureAdapter(fake, MmcifRuntime));
}

/** A viewer registry with the default adapter of mmCIF: the runtime reading the file through the assets port. */
export function withDefaultMmcif<E extends { readonly kind: string }>(registry: readonly E[]) {
  return withAdapter(registry, MmcifRuntime);
}

/**
 * What the store holds of the host's source, as lines: each node by kind and key (a sequence with its length), the
 * complex with its units. What the story prints beside the view and a twin reads.
 */
export function heldLines(viewer: SeqViewer): string[] {
  const nodes = viewer.store.store.coverage(HostSource)?.nodes ?? [];
  if (!nodes.length) return ['nothing held yet'];
  return nodes.map(({ key, node }) => {
    const p = (node.params ?? {}) as Record<string, unknown>;
    if (node.kind === 'complex') {
      const units = (node.children ?? []).filter((c) => c.kind === 'unit').map((u) => u.id);
      return `complex ${key}: ${units.length} units ${units.join(' ')}`;
    }
    if (node.kind === 'sequence') return `sequence ${key}: ${String(p.letters ?? '').length} residues`;
    if (node.kind === 'molecule') return `molecule ${key}: ${String(p.name ?? '')}`;
    return `${node.kind} ${key}`;
  });
}
