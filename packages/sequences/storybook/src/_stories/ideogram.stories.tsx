/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { useMemo } from 'react';
import f from './support/seqstar-format.module.css';
import {
  IDEOGRAM_PLACES,
  IDEOGRAM_REGISTRY,
  IDEOGRAM_SECTION,
  ideogramDoc,
  ideogramLoaders,
} from './support/seqstar-ideogram.js';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard } from './support/seqstar-viewer.js';
import { TextAction } from './ui/index.js';

/**
 * The ideogram navigator: each block's chromosome drawn as its karyotype strip from hg38's committed cytobands (UCSC's
 * cytoBandIdeo, 862 bands recorded once), the view's window over it. The data is the story's (a registry entry,
 * `ideogramNavigator(bands)`); the navigator knows no biology but the stains. Offline: the sequence is synthetic.
 */
const meta = { title: 'Seq*/Views and layouts/Navigators', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

function Panel() {
  const { context: ctx } = useSeqViewer();
  return (
    <div className={f.row}>
      {/* biome-ignore lint/suspicious/noFocusedTests: the viewport's verb, not a focused test. */}
      <TextAction onClick={() => ctx.viewports.get(IDEOGRAM_SECTION)?.fit()}>Whole genome</TextAction>
      {IDEOGRAM_PLACES.map(([label, block, a, e]) => (
        <TextAction key={label} onClick={() => ctx.viewports.get(IDEOGRAM_SECTION)?.zoomTo(block, a, e)}>
          {label}
        </TextAction>
      ))}
    </div>
  );
}

export const Karyotype: Story = {
  name: 'hg38 karyotype',
  render: () => {
    const doc = useMemo(ideogramDoc, []);
    const options = useMemo(
      () => ({ loaders: ideogramLoaders(), maxHeight: 420, navigator: 'ideogram', registry: IDEOGRAM_REGISTRY }),
      [],
    );
    return (
      <PrismPage
        title="The ideogram navigator: hg38's chromosomes as their bands"
        caption={
          <>
            The row above the tracks is the navigator: a chromosome a block, each drawn as its cytobands, Giemsa stains
            as grays from light (negative) to dark (positive 100), variable and stalk bands lighter, the centromere's
            two bands as red triangles meeting in the middle. The window is the view's: a press jumps it, dragging pans
            it, an end zooms, a double click zooms in on the pointer and back out. Zoom to a chromosome and its strip
            fills the row; the other chromosomes keep their share of it. The cytobands are a committed fixture of UCSC's
            cytoBandIdeo; the sequence, genes and variants are synthetic. It is a navigator registered as an entry (the
            same one the streaming genome browser case uses), so any document can name it.
          </>
        }
      >
        <ViewerCard doc={doc} options={options} label="24 chromosomes, one block each">
          <Panel />
        </ViewerCard>
      </PrismPage>
    );
  },
};
