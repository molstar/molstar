/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/react-vite';
import { MolstarDemo } from './support/seqstar-molstar.js';
import { PrismPage } from './support/seqstar-prism.js';

/**
 * The Mol* adapter: P04637 in the Seq* viewer beside a Mol* canvas. The structures are complexes of the document,
 * mapped onto P04637 by their SIFTS numbering (the AlphaFold model by identity). The PDB coverage nests: collapsed, its
 * depth (a `summary` track: how many structures cover each position); open, a row per structure (the union of its
 * chains, `display: "union"`), each split by chain. Every coverage row names its structure in a JSON field, so the
 * viewer's hover card offers "Show in 3D" on it (a host action), and the track menu offers "Show … in 3D" and "Color
 * 3D by this track". What 3D shows is an MVS state composed from the viewer's state at that moment: the structure, the
 * colors of the chosen track as drawn (`color_from_uri` tables, sites also as ball and stick), the selection. SeqQL
 * reaches the structure through the link: `molql` runs a PyMOL selection in Mol* (`molstarQueryExtension`), and Mol*'s
 * selection is the `structure` input.
 */
const meta = { title: 'Seq*/Cases/Mol*', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

export const P53: Story = {
  name: 'p53: sequence and structure',
  render: (_args, { globals }) => (
    <PrismPage
      title="p53: sequence and structure"
      caption={
        <>
          <b>Heavily work in progress.</b> The Mol* link is a proof of concept: it will be reworked in Mol*'s own
          repository, and its behavior and look will change. UniProt P04637 with its features, variants, a hydropathy
          profile and the AlphaFold model's pLDDT. The PDB coverage, collapsed, shows how many structures cover each
          position; open it, and rest the pointer on a structure's bar: its card offers <b>Show in 3D</b> (open a
          structure too, to show one chain in front). Color the structure by any track that describes the sequence (a
          track's ⋯ menu, or the chooser over the canvas): domains, regions, hydropathy, the residues' letter scheme;
          sites (PTMs, binding sites, variants) also show as ball and stick. pLDDT describes the AlphaFold model only,
          so it is offered (and the default) there alone. Hover or select residues in either view: the other follows,
          and a selection focuses it. SeqQL queries reach the structure: <code>(molql "…")</code> runs a PyMOL selection
          in Mol* and gives its residues, <code>(input structure)</code> is what is selected in 3D. The structure in 3D
          has its PDB coverage entry opened and marked <b>· in 3D</b>: a row per chain (by its author's chain ID) where
          SIFTS puts its residues on P04637. A chain's residues light its own row alone, hovered in either view or found
          by a query; P04637's light every row. A result stays on the chains unless it is put <code>onto p53</code>,
          which means any chain; <b>Select in 3D</b> sends it back.
        </>
      }
    >
      <MolstarDemo theme={globals.theme === 'dark' || globals.theme === 'blue' ? globals.theme : 'light'} />
    </PrismPage>
  ),
};
