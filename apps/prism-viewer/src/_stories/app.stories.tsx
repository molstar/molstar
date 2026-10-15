/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { unzip } from '@molstar/core/util/zip/zip';
import { story } from '@molstar/prism-ui/testing/index';
import { RealViewerStory } from '@molstar/prism-viewer/testing/real-viewer';
import { FullPrismViewerSpec } from '../full-spec.js';

// The viewer as the app composes it (the viewer's defaults and Mol*'s extensions: `FullPrismViewerSpec`), for what
// needs an extension. These reach the network.

const meta: Meta = { title: 'Viewer/With extensions', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * A predicted structure with ligands: AlphaFill's model of CYP3A4 (AlphaFold's AF-P08684-F1 with ligands transplanted
 * from PDB homologs), the cartoon colored by pLDDT (the model-archive extension's annotation), the sequence's pLDDT track
 * in the same bands, linked to the 3D view. Network.
 */
export const Predicted: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => FullPrismViewerSpec()}
        setup={async (v) => {
          const loaded = await v.api.load.url!({ url: 'https://alphafill.eu/v1/aff/P08684', format: 'mmcif' });
          if (!loaded.ok) throw new Error(loaded.message);
          await v.api.annotations.apply({ annotation: 'plddt' });
        }}
      />
    )),
};

/**
 * The predicted aligned error of a model (the model-archive extension's PAE card): AlphaFold DB's model of calmodulin
 * (P0DP23, 149 residues: two lobes on a flexible linker, so two confident blocks) with its PAE written into it as ModelCIF's pairwise metric (`ma_qa_metric_local_pairwise`,
 * as ModelArchive files carry it; AlphaFold DB keeps it in a JSON beside the model), loaded, the PAE card open. Network.
 */
export const Pae: StoryObj = {
  name: 'PAE',
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => FullPrismViewerSpec()}
        setup={async (v) => {
          const [entry] = await (await fetch('https://alphafold.ebi.ac.uk/api/prediction/P0DP23')).json();
          const [cif, pae] = await Promise.all([
            fetch(entry.cifUrl).then((r) => r.text()),
            fetch(entry.paeDocUrl).then((r) => r.json()),
          ]);
          const matrix: number[][] = (Array.isArray(pae) ? pae[0] : pae).predicted_aligned_error;
          const loaded = await v.api.load.data!({ data: withPairwisePae(cif, matrix), format: 'mmcif' });
          if (!loaded.ok) throw new Error(loaded.message);
          await v.api.view.open({ card: 'model-archive.overlay.pae' });
        }}
      />
    )),
};

/**
 * The predicted aligned error of a complex: ModelArchive's ma-bak-cepc-0100, a predicted interaction between two yeast
 * ribosomal proteins (RPS4B, chain A, 261 residues; RPS24B, chain B, 135), its PAE taken from the entry's accompanying
 * archive (ModelArchive keeps it there, a partial ModelCIF of `ma_qa_metric_local_pairwise`) and put into the model,
 * loaded, the PAE card open: the plot's chain bands and boundaries, a box across them. Network.
 */
export const PaeComplex: StoryObj = {
  name: 'PAE (complex)',
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => FullPrismViewerSpec()}
        setup={async (v) => {
          const entry = 'https://modelarchive.org/api/projects/ma-bak-cepc-0100';
          const [cif, archive] = await Promise.all([
            fetch(`${entry}?type=basic__model_file_name`).then((r) => r.text()),
            fetch(`${entry}?type=materials_procedures__accompanying_data_file_name`).then((r) => r.arrayBuffer()),
          ]);
          const files = await unzip(SyncRuntimeContext, archive);
          const name = Object.keys(files).find((f) => f.endsWith('_predicted_aligned_error_v1.cif'));
          const file = name && files[name];
          if (!(file instanceof Uint8Array)) throw new Error('The archive has no PAE.');
          const pairwise = /loop_\n_ma_qa_metric_local_pairwise\.[\s\S]*?(?=\nloop_|\n#|$)/.exec(
            new TextDecoder().decode(file),
          );
          if (!pairwise) throw new Error('The PAE file has no ma_qa_metric_local_pairwise loop.');
          const loaded = await v.api.load.data!({ data: `${cif.trimEnd()}\n${pairwise[0]}\n#\n`, format: 'mmcif' });
          if (!loaded.ok) throw new Error(loaded.message);
          await v.api.view.open({ card: 'model-archive.overlay.pae' });
        }}
      />
    )),
};

/**
 * A ModelCIF with a PAE matrix as its pairwise metric: a row in its `_ma_qa_metric` loop (mode `local-pairwise`), and
 * an `_ma_qa_metric_local_pairwise` loop of every residue pair (one chain, `asym`; residues by `label_seq_id` from 1).
 */
function withPairwisePae(cif: string, matrix: readonly (readonly number[])[], asym = 'A'): string {
  const loop = /loop_\n((?:_ma_qa_metric\.[^\n]+\n)+)([\s\S]*?)\n#/.exec(cif);
  if (!loop) throw new Error('The model has no _ma_qa_metric loop.');
  const fields = loop[1]
    .trim()
    .split('\n')
    .map((f) => f.trim().slice('_ma_qa_metric.'.length));
  const ids = loop[2]
    .trim()
    .split('\n')
    .map((row) => Number(row.trim().split(/\s+/)[fields.indexOf('id')]));
  const id = Math.max(0, ...ids.filter(Number.isFinite)) + 1;
  const values: Record<string, string> = { id: String(id), mode: 'local-pairwise', name: 'PAE', type: 'PAE' };
  const metric = fields.map((f) => values[f] ?? '.').join(' ');
  const withMetric = cif.replace(loop[0], `loop_\n${loop[1]}${loop[2]}\n${metric}\n#`);
  const rows: string[] = [];
  let n = 0;
  matrix.forEach((row, i) => {
    row.forEach((value, j) => {
      rows.push(`${++n} 1 ${asym} ${i + 1} . ${asym} ${j + 1} . ${id} ${value.toFixed(2)}`);
    });
  });
  const pairwise = [
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
    ...rows,
    '#',
  ].join('\n');
  return `${withMetric.trimEnd()}\n${pairwise}\n`;
}
