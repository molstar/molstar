/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The committed structure fixtures (./fixtures/structures; its README says where they come from) and the loaders the
// load-path tests read them through: PDBe's and AlphaFold DB's URIs as the services resolve them, P04637's FASTA from
// Seq*'s shared fixture (p53.json, read as data). No network.

import { readFileSync } from 'node:fs';
import { createSeqViewSpec } from '@molstar/seq-core/builder';
import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { createRegistries, KernelRegistry, type Loaders, MvsRegistry, register } from '@molstar/seq-core/format';
import { MolstarSvsEntry } from '../index.js';

const ref = (path: string) => new URL(`../../../../../../sequences/core/_test/fixtures/${path}`, import.meta.url);
export const bytesOf = (file: string) =>
  new Uint8Array(readFileSync(new URL(`./fixtures/structures/${file}`, import.meta.url)));
export const P53: { sequence: string } = JSON.parse(readFileSync(ref('p53.json'), 'utf8'));

const fasta = (header: string, seq: string) => `>${header}\n${seq.replace(/(.{60})/g, '$1\n')}\n`;
const FILES: Record<string, () => Uint8Array> = {
  'https://rest.uniprot.org/uniprotkb/P04637.fasta': () =>
    new TextEncoder().encode(fasta('sp|P04637|P53_HUMAN Cellular tumor antigen p53', P53.sequence)),
  ...Object.fromEntries(
    ['1tup', '2ocj', '2xwr', '1c26'].map((id) => [
      `https://www.ebi.ac.uk/pdbe/entry-files/download/${id}.bcif`,
      () => bytesOf(`${id}.bcif`),
    ]),
  ),
  'https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.bcif': () => bytesOf('AF-P04637-F1-model_v6.bcif'),
};

export function loaders(): Loaders {
  return {
    async load(request) {
      const read = FILES[request.uri ?? request.hostFile ?? ''];
      if (!read) throw new Error(`no fixture for ${request.uri ?? request.hostFile}`);
      return read();
    },
  };
}

export function registries() {
  const reg = createRegistries();
  register(reg, [KernelRegistry, MvsRegistry, MolstarSvsEntry]);
  return reg;
}

/**
 * P04637 with structures as complexes (`pdbe` ids, or the AlphaFold model `af`), each mapped onto it: a PDB entry's p53
 * entity by its SIFTS numbering, the model's by identity; the model's pLDDT by `molstar:plddt` and from its
 * `ma_qa_metric_local` category (case 8's route).
 */
export function p53Doc(structures: readonly { ref: string; pdb?: string; entity?: string }[]): SeqViewSpecDoc {
  const b = createSeqViewSpec({ title: 'p53' });
  const p53 = b.data.sequence({ ref: 'p53', source: { service: 'uniprot', resource: 'sequence', id: 'P04637' } });
  for (const s of structures) {
    const af = s.pdb === undefined;
    b.data.complex({
      ref: s.ref,
      provider: 'molstar:structure',
      source: af
        ? { service: 'afdb', resource: 'structure', id: 'AF-P04637-F1' }
        : { service: 'pdbe', resource: 'structure', id: s.pdb as string },
    });
    const to = `${s.ref}/entity:${s.entity ?? '1'}`;
    b.data.mapping(
      (af
        ? { ref: `${s.ref}-map`, from: p53, to, method: 'identity' }
        : { ref: `${s.ref}-map`, from: p53, to, method: 'by_numbering', numbering: 'uniprot:P04637' }) as never,
    );
    if (af) {
      b.data.computed({ ref: `${s.ref}-plddt`, provider: 'molstar:plddt', options: { of: s.ref } } as never);
      b.data.annotation({
        ref: `${s.ref}-qa`,
        geometry: 'values',
        from_source: s.ref,
        category_name: 'ma_qa_metric_local',
        numbering: 'label',
        fields: { target: 'label_asym_id', pos: 'label_seq_id', value: 'metric_value' },
      } as never);
    }
  }
  b.view({ ref: 'tracks', label: 'Tracks' }).track({ label: '' }).ruler();
  return b.getDocument();
}
