/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { dropStory } from './drop-root.js';
import { servicesStory } from './services-root.js';

// Services (Examples, Showcases): an application composing UniProt, PDBe (SIFTS) and AlphaFold DB from recorded answers
// with a switch to the services themselves, a lab's predictions onto a structure, addressing with a structure beside the
// view, and a view of your own file. Each runs on the fake structure viewer (the Mol* view waits for the new mol viewer).

const meta: Meta = { title: 'Examples/Services' };
export default meta;

/**
 * p53 from three services, each answered from its recorded file after a delay so the order of arrival is seen: UniProt's
 * entry and AlphaFold DB's confidence draw first, PDBe's mmCIF and SIFTS's mappings after; SIFTS is read by the
 * application into one `structure` map a chain, and the fake structure viewer is linked once they are all here. "Switch to
 * the live services" starts again with the same spec, the answers now from the services themselves through the asset
 * policy's hosts (`rest.uniprot.org`, `www.ebi.ac.uk`, `alphafold.ebi.ac.uk`, `https:` only): a service that does not
 * answer leaves what arrived drawn and says which failed. Hover R175: it lights on the four chains of 2OCJ. Hover a cell of a chain in the fake: it lights on P04637
 * and on the structure's B-factor track, which follows P04637 through SIFTS's map, so it reads under residues 94 to 312.
 */
export const ComposingServices: StoryObj = {
  render: () =>
    servicesStory({
      scenario: 'composed',
      hint: 'UniProt, AlphaFold DB, PDBe and SIFTS answer from their recorded files, each after 600 ms. Watch the view fill in, then hover a residue of P04637: it lights on the four chains of PDB 2OCJ in the fake. The switch asks the services themselves (the same spec).',
    }),
};

/**
 * A lab's scores on a sequence it indexes from 1 (made up: synthetic), mapped by an offset of 0 onto UniProt's P04637 and
 * by SIFTS's maps onto 2OCJ's chains: two maps in a chain. Hover a position of the lab's sequence: the residue lights on
 * every chain in the fake (the unmodeled ends light nothing). The cube icon in the "Predictions" track's header ("Color the structure by this track") colors the
 * fake's chains by the scores (a snapshot of the track; again restores).
 */
export const PredictionsOntoAStructure: StoryObj = {
  render: () =>
    servicesStory({
      scenario: 'predictions',
      delayMs: 300,
      hint: 'The scores are SYNTHETIC (generated from a seed, high at R175, R248 and R273 and over the DNA-binding domain). Hover a position of the lab’s sequence: it lights on the four chains in the fake. Click the cube icon (Color the structure by this track; the palette icon is the track’s own colors) in the "Predictions" track’s header: the fake’s chains are colored by the track.',
    }),
};

/**
 * A text whose links focus, highlight and select in the view and in the fake: the DNA-binding domain lights on the four
 * chains, the hotspots too; one residue is named three ways (UniProt's R175, residue 82 of the entity, 175 in chain A's
 * author numbers) and every link reaches the same place. "Copy link" turns the selection (or what the pointer is over)
 * back into a link.
 */
export const AddressingWithAStructure: StoryObj = {
  render: () =>
    servicesStory({
      scenario: 'addressing',
      delayMs: 300,
      text: true,
      hint: 'Click a link in the text: the view focuses or selects it and the fake lights it on every chain. Hover a link: its residues are highlighted in the view. Select residues in the view and press "Copy link": the link that names them.',
    }),
};

/**
 * Drop a file of your own on the view (or pick one): `.seq`, FASTA, Stockholm, aligned FASTA, TSV, CSV, JSON, MessagePack,
 * GFF3, mmCIF or GenBank (`.gb`, `.gbk`: each record a sequence with its features, a circular one a ring's sequence).
 * Nothing leaves the page: the viewer's assets port answers the dropped bytes at one location. The "Try" buttons load a
 * sample of each format as if dropped ("Try pUC19.gb": the examples' generated plasmid record).
 */
export const YourOwnFile: StoryObj = {
  render: () =>
    dropStory({
      hint: 'Drop a .seq, FASTA, Stockholm, TSV, CSV, JSON, GFF3, mmCIF or GenBank file on the box, or use a "Try" button. The file is read in the page and nothing else is.',
    }),
};
