/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GenomesViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { storyId } from '../../story-id.js';

const meta: Meta = { title: 'Examples/Genomes/Index' };
export default meta;

/** Each view's story export (`views.stories.ts`). */
const GenomesExports: Readonly<Record<string, string>> = {
  synteny: 'Synteny',
  'dot-plot': 'DotPlot',
  signal: 'Signal',
  samples: 'Samples',
  methylation: 'Methylation',
  transcripts: 'Transcripts',
  'by-name': 'ByName',
  names: 'NamesAliased',
  'contact-map': 'ContactMap',
  liftover: 'Liftover',
};

const Sources: readonly (readonly [title: string, summary: string])[] = [
  [
    'Generated (LicenseRef-seqx-generated)',
    'two toy genomes (A prefixed, B unprefixed), A′ (A revised), their alignment and map in the seq format, 150 genes as GFF3 and GTF, a bedGraph, a VCF of 40 samples (the four bgzipped with their tabix indexes), methylation, the dot plot’s and the Hi-C’s grid tiles: made by `generateGenomes` from one seed (core’s PCG and hash3), byte-identical to the committed files; the letters of both toys from an in-memory adapter, never stored.',
  ],
  [
    'Recorded (the formats’ stories)',
    'UCSC hg38 cytoBandIdeo, Ensembl 116’s GFF3 and GENCODE 50’s GTF around TP53, NCBI ClinVar of 2024-01-07 with its index: their licences and attributions in `fixtures/genomes/sources.json`.',
  ],
];

function item(id: string, title: string, summary: string): HTMLElement {
  const li = document.createElement('li');
  const a = document.createElement('a');
  a.href = `?path=/story/${storyId('Examples/Genomes/Views', GenomesExports[id])}`;
  a.target = '_top';
  a.textContent = title;
  li.append(a, ` ${summary}`);
  return li;
}

/** The views of the Genomes dataset, and its sources. */
export const Views: StoryObj = {
  render: () => {
    const root = document.createElement('div');
    root.style.cssText = 'font: 14px/1.5 system-ui, sans-serif; padding: 24px; max-width: 72ch;';
    const h = document.createElement('h2');
    h.textContent = 'Genomes';
    const ul = document.createElement('ul');
    ul.append(...GenomesViews.map((v) => item(v.id, v.title, v.summary)));
    const h2 = document.createElement('h3');
    h2.textContent = 'Sources';
    const sources = document.createElement('ul');
    for (const [title, summary] of Sources) {
      const li = document.createElement('li');
      li.append(`${title}: ${summary}`);
      sources.append(li);
    }
    const ring = document.createElement('p');
    ring.textContent = 'The genome ring and both genomes on one ring are the ring extension’s stories.';
    root.append(h, ul, h2, sources, ring);
    return root;
  },
};
