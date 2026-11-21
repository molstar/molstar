/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CoordinatesViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// Coordinates (the coordinate model): one axis read through a map, on the gene and its protein. A section placed on one
// system follows another through the map between them, and reads backward or collapsed on its own terms. The
// datasets' figures with a gene are their examples.

const meta: Meta = { title: 'Features/Coordinates' };
export default meta;

const view = (id: string) => P53Views.find((v) => v.id === id)!;

/**
 * One gene section: the protein's letters a codon wide, its domains and the gene's letters, read along the protein on
 * its coding pieces. Each intron is a hairline with a notch (hover it: the hidden stretch); a click on a residue
 * selects it, its codon shown; Shift+R reverses the section, Shift+C shows the introns.
 */
export const ProteinAndGeneOnOneAxis: StoryObj = { render: () => exampleStory(view('one-axis')) };

/**
 * For contrast: the gene follows its protein in genomic order, its introns shown. A pan of the protein moves the gene
 * the other way (minus strand); in Examples › p53 › Protein on its gene both move together.
 */
export const GenomicFollower: StoryObj = { render: () => exampleStory(view('genomic-follower')) };

/**
 * The gene alone, reversed: positions count down, its first label says so, the letters are the other strand's (the
 * complement). Shift+C collapses the introns to hairlines; Shift+R turns it back.
 */
export const ReversedGene: StoryObj = { render: () => exampleStory(view('reversed-gene')) };

/**
 * A protein over its coding DNA (a short made-up gene, plus strand, genomic order, nothing collapsed): each residue sits
 * over its three nucleotides (the intron has no protein letters), the exons are coloured, the domain is cut at the
 * intron and joined across it by a connector, and a phosphosite marks residue 17. Hover a residue: it and its codon
 * light, the card names the residue and its codon positions.
 */
export const ProteinOverItsCodingDna: StoryObj = {
  name: 'Protein Over Its Coding DNA',
  render: () => exampleStory(CoordinatesViews.find((v) => v.id === 'coding-dna')!),
};
