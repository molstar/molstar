/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { zip } from '@molstar/core/util/zip/zip';
import { CIF } from '@molstar/io/reader/cif';

// A tiny ModelCIF whose PAE is in an associated file (Model Archive's way, as ma-bak-cepc-0003): chain A, 3 residues;
// the PAE is 2·|i − j| + 1. The model names it in a zip with the MSA, or as a CIF of its own.

export const Residues = [
  [1, 'ALA'],
  [2, 'GLY'],
  [3, 'SER'],
] as const;
export const pae = (i: number, j: number) => 2 * Math.abs(i - j) + 1;
export const Zip =
  'https://modelarchive.org/api/projects/ma-test-0001?type=materials_procedures__accompanying_data_file_name';
export const Path = 'ma-test-0001_predicted_aligned_error_v1.cif';

const metrics = [
  'loop_',
  '_ma_qa_metric.id',
  '_ma_qa_metric.name',
  '_ma_qa_metric.type',
  '_ma_qa_metric.mode',
  '_ma_qa_metric.software_group_id',
  '1 PAE PAE local-pairwise .',
  '#',
];

/** The model, naming its PAE file: `zip` (an entry of an archive), `cif` (a file of its own, at `url`), or none. */
export function modelCif(associated?: { kind: 'zip' } | { kind: 'cif'; url: string }): string {
  const atoms = Residues.map(
    ([seq, comp], k) =>
      `ATOM ${k + 1} C CA . ${comp} A 1 ${seq} ? ${(3.8 * k).toFixed(3)} 0.000 0.000 1 90.0 ${seq} A 1`,
  );
  const files =
    associated?.kind === 'zip'
      ? [
          'loop_',
          ...['id', 'entry_id', 'file_url', 'file_type', 'file_format', 'file_content', 'details'].map(
            (f) => `_ma_entry_associated_files.${f}`,
          ),
          `1 ma-test-0001 '${Zip}' archive zip 'archive with multiple files' .`,
          '#',
          'loop_',
          ...['id', 'archive_file_id', 'file_path', 'file_format', 'file_content', 'description'].map(
            (f) => `_ma_associated_archive_file_details.${f}`,
          ),
          `1 1 ${Path} cif 'local pairwise QA scores' ?`,
          "2 1 A.a3m a3m 'multiple sequence alignments' ?",
          '#',
        ]
      : associated?.kind === 'cif'
        ? [
            'loop_',
            ...['id', 'entry_id', 'file_url', 'file_type', 'file_format', 'file_content'].map(
              (f) => `_ma_entry_associated_files.${f}`,
            ),
            `1 ma-test-0001 '${associated.url}' file cif 'local pairwise QA scores'`,
            '#',
          ]
        : [];
  return [
    'data_ma-test-0001',
    '_entry.id ma-test-0001',
    '#',
    'loop_',
    '_entity.id',
    '_entity.type',
    '1 polymer',
    '#',
    'loop_',
    ...[
      'group_PDB',
      'id',
      'type_symbol',
      'label_atom_id',
      'label_alt_id',
      'label_comp_id',
      'label_asym_id',
      'label_entity_id',
      'label_seq_id',
      'pdbx_PDB_ins_code',
      'Cartn_x',
      'Cartn_y',
      'Cartn_z',
      'occupancy',
      'B_iso_or_equiv',
      'auth_seq_id',
      'auth_asym_id',
      'pdbx_PDB_model_num',
    ].map((f) => `_atom_site.${f}`),
    ...atoms,
    '#',
    // The model lists the metric; its values are in the associated file.
    ...metrics,
    ...files,
    '',
  ].join('\n');
}

/** The associated file: the pairwise metric over the model's residues (`chain`: another model's chain, fitting none). */
export function paeCif(chain = 'A'): string {
  const pairs: string[] = [];
  let n = 0;
  Residues.forEach(([s1, r1], i) => {
    Residues.forEach(([s2, r2], j) => {
      pairs.push(`${++n} 1 ${chain} ${s1} ${r1} ${chain} ${s2} ${r2} 1 ${pae(i, j).toFixed(2)}`);
    });
  });
  return [
    'data_ma-test-0001_predicted_aligned_error_v1',
    ...metrics,
    'loop_',
    ...[
      'ordinal_id',
      'model_id',
      'label_asym_id_1',
      'label_seq_id_1',
      'label_comp_id_1',
      'label_asym_id_2',
      'label_seq_id_2',
      'label_comp_id_2',
      'metric_id',
      'metric_value',
    ].map((f) => `_ma_qa_metric_local_pairwise.${f}`),
    ...pairs,
    '#',
    '',
  ].join('\n');
}

const bytes = (s: string) => new TextEncoder().encode(s);
/** A zip of these files; `stored`: not compressed (jsdom's Blob has no stream for the compressor). */
export const archive = async (entries: Record<string, string>, stored = false) =>
  zip(SyncRuntimeContext, Object.fromEntries(Object.entries(entries).map(([k, v]) => [k, bytes(v)])), stored);

export async function frameOf(cif: string) {
  const parsed = await CIF.parse(cif).run();
  if (parsed.isError) throw new Error(parsed.message);
  return parsed.result.blocks[0];
}
