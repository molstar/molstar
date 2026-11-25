/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The seq format's molecule, compound and complex tables (_spec/README.md, Data, The seq format): a source's nodes of
// the molecule kinds carry their typed tables through a file, JSON and MessagePack, and read back as nodes the
// language takes under `data`, checked by its rules as every dynamic node is.

import { MsgPackCodec } from '@molstar/core/repr/codecs/msgpack';
import { validateAndNormalize } from '@molstar/core/repr/schema';
import { hashNode, type Node } from '@molstar/core/repr/tree';
import {
  checkDataNode,
  createSeqSpec,
  identifiers,
  modified,
  SeqLanguage,
  SeqLanguageVersion,
  SeqValidateOptions,
} from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { decodeSeq, encodeSeq, readSeq, type SeqData, seqToNode } from '../seq-format.js';

/** What an mmCIF entry's reading gives, as dynamic nodes: an entity a sequence or a molecule, compounds, a complex. */
function entryNodes(): readonly Node[] {
  const b = createSeqSpec();
  b.data.sequence('1', {
    alphabet: 'rna',
    letters: 'GCGGAUUUAGCUCAG',
    modified: modified([
      { position: 9, component: '2MG' },
      { position: 12, component: 'H2U' },
    ]),
  });
  b.data.molecule(
    '2',
    {
      class: 'branched',
      name: 'NAG-NAG',
      residues: ['NAG', 'NAG'],
      atoms: { residue: [0, 0, 1], name: ['C1', 'O4', 'C1'], element: ['C', 'O', 'C'] },
      identifiers: identifiers({ wurcs: 'WURCS=2.0/1,2,1/…' }),
    },
    (m) =>
      m.links('links', {
        rows: { a_start: [1], a_end: [2], b_start: [0], b_end: [1], a_atom: ['C1'], b_atom: ['O4'] },
      }),
  );
  b.data.compound('2MG', { name: '2N-methylguanosine-5’-monophosphate', parents: ['G'], class: 'RNA LINKING' });
  b.data.complex('1abc', {}, (c) => {
    c.unit(
      'A',
      { of: '/data/src/(node :key 1)', chain: 'A', operator: '1_555', observed: { start: [0], end: [12] } },
      (u) =>
        u.numbering('auth', {
          method: 'labels',
          scheme: 'auth',
          labels: Array.from({ length: 15 }, (_, i) => `${i + 1}`),
        }),
    );
    c.unit('C', { of: '/data/src/(node :key 2)', chain: 'C' });
    c.links('links', {
      rows: { a_unit: ['A'], a_start: [3], a_end: [4], b_unit: ['C'], b_start: [0], b_end: [1], type: ['covalent'] },
    });
  });
  const data = b.build().children?.[0];
  return data?.children ?? [];
}

function fileOf(nodes: readonly Node[]): SeqData {
  const keyed = nodes.map((n) => ({ key: n.id ?? n.kind, node: n }));
  return {
    coverage: { systems: [{ name: '1', length: 15 }], shapes: [], levels: [1], nodes: keyed },
    pairs: [{ request: { system: '1', start: 0, end: 15, level: 1, shapes: [] }, answer: { ok: true, tables: {} } }],
  };
}

describe('the seq format’s molecule, compound and complex nodes', () => {
  const data = fileOf(entryNodes());

  it('round-trip through JSON and MessagePack, every table typed', () => {
    const text = encodeSeq(data);
    const json = decodeSeq(text);
    expect(json.issues).toEqual([]);
    expect(hashNode(seqToNode(json.data!))).toBe(hashNode(seqToNode(data)));
    const packed = readSeq(MsgPackCodec.decode(MsgPackCodec.encode(seqToNode(data))).node as Node);
    expect(packed.issues).toEqual([]);
    expect(hashNode(seqToNode(packed.data!))).toBe(hashNode(seqToNode(data)));
    const molecule = json.data!.coverage.nodes!.find((n) => n.key === '2')!.node;
    const atoms = molecule.params?.atoms as Record<string, unknown>;
    expect(atoms.residue).toBeInstanceOf(Int32Array);
  });

  it('read back as nodes the language takes under data, passing its rules (as the overlay checks them)', () => {
    for (const { node } of decodeSeq(encodeSeq(data)).data!.coverage.nodes!) {
      const spec: Node = {
        kind: 'spec',
        params: { version: SeqLanguageVersion },
        children: [{ kind: 'data', id: 'data', children: [node] }],
      };
      const read = validateAndNormalize(SeqLanguage, spec, SeqValidateOptions);
      expect(read.issues).toEqual([]);
      const said: string[] = [];
      checkDataNode(node, { target: () => undefined }, (code) => said.push(code));
      expect(said).toEqual([]);
    }
  });

  it('a host’s molecule with a hydrogen is said, as a spec’s is (EC-64)', () => {
    const b = createSeqSpec();
    b.data.molecule('w', { residues: ['HOH'], atoms: { residue: [0, 0], name: ['O', 'H1'], element: ['O', 'H'] } });
    const said: string[] = [];
    checkDataNode(b.build().children![0].children![0], { target: () => undefined }, (code) => said.push(code));
    expect(said).toEqual(['seq.lang.molecule.hydrogen']);
  });
});
