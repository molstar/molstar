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
import type { PrismViewer } from '@molstar/prism-viewer';
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
  render: () => story(() => <RealViewerStory spec={() => FullPrismViewerSpec()} setup={loadPredicted} />),
};

/**
 * The embedded layout (`ui.layout.compactRail`, the app's `?hide-rail`): the rail hidden behind its button in the
 * canvas's corner, the canvas edge to edge; the button slides the rail in over the canvas, and a panel opens beside it.
 * For a viewer inside a page (the Seq* introduction's). The same structure as Predicted. Network.
 */
export const Embedded: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => {
          const spec = FullPrismViewerSpec();
          return { ...spec, ui: { ...spec.ui, layout: { ...spec.ui?.layout, compactRail: true } } };
        }}
        setup={loadPredicted}
      />
    )),
};

/** AlphaFill's CYP3A4 with ligands, colored by pLDDT. */
async function loadPredicted(v: PrismViewer) {
  const loaded = await v.api.load.url!({ url: 'https://alphafill.eu/v1/aff/P08684', format: 'mmcif' });
  if (!loaded.ok) throw new Error(loaded.message);
  await v.api.annotations.apply({ annotation: 'plddt' });
}

/**
 * The predicted aligned error of a model (the model-archive extension's PAE card): AlphaFold DB's model of calmodulin
 * (P0DP23, 149 residues: two lobes on a flexible linker, so two confident blocks), loaded by its address, its PAE
 * fetched from AlphaFold DB as the structure row's "Fetch PAE" does (`model-archive.fetch-pae`), the
 * PAE card open. Network.
 */
export const Pae: StoryObj = {
  name: 'PAE',
  render: () =>
    story(() => (
      <RealViewerStory
        spec={() => FullPrismViewerSpec()}
        setup={async (v) => {
          const [entry] = await (await fetch('https://alphafold.ebi.ac.uk/api/prediction/P0DP23')).json();
          const loaded = await v.api.load.url!({ url: entry.cifUrl, format: 'mmcif' });
          if (!loaded.ok) throw new Error(loaded.message);
          const fetched = await v.api['model-archive']!.fetchPae({});
          if (!fetched.ok) throw new Error(fetched.message);
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
