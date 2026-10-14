/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import { useSeqState, useSeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { type ReactNode, useMemo } from 'react';
import {
  alignmentGallery,
  customDoc,
  HOST,
  manyDoc,
  NUCLEIC,
  PROTEIN,
  pickerDoc,
  rangeDoc,
  sharedDoc,
  valuesGallery,
} from './support/seqstar-docs-colors.js';
import f from './support/seqstar-format.module.css';
import { PrismCanvas, PrismPage, prism } from './support/seqstar-prism.js';
import { ViewerCard, type ViewerOptions } from './support/seqstar-viewer.js';
import { TextAction } from './ui/index.js';

/**
 * Color schemes: named schemes as data, picked in the document (`{"theme": "scheme"}`), in the builder
 * (`b.colors`), at runtime (a legend group's title, a track's Colors ▸) and by command (`seqv:colorScheme`). Every
 * scheme on a protein alignment (case 3's SH2 family), on an RNA alignment (case 4's tRNAs), on values and on
 * categories, light and dark; a live picker with a host's own scheme; and a document's own palettes.
 */
const meta = { title: 'Seq*/Themes/Color schemes', parameters: { layout: 'fullscreen' } } satisfies Meta;
// The pages follow the toolbar's theme; the viewers' canvases are pinned light and dark side by side (comparing a
// scheme's two variants is what these stories show), the live picker's follows the UI theme.
export default meta;
type Story = StoryObj;

/** What these stories are about, left out by `clean` chrome: the legend (`TrackActions` too where the track menu's Colors ▸ is used). */
const LEGEND = ['Legend'] as const;

/** The canvas pinned light or dark, whatever the UI theme (seqstar-prism.tsx). */
const lightCanvas = (viewer: ReactNode) => <PrismCanvas theme="light">{viewer}</PrismCanvas>;
const darkCanvas = (viewer: ReactNode) => <PrismCanvas theme="dark">{viewer}</PrismCanvas>;

/**
 * The same document on a light and a dark canvas: comparing a scheme's two variants is what these stories show, so the
 * canvases are pinned (Prism themes a canvas apart from the UI); the page around them follows the toolbar.
 */
function LightAndDark({ doc, title, options }: { doc: SeqViewSpecDoc; title: string; options?: ViewerOptions }) {
  return (
    <div className={prism.pair}>
      <ViewerCard doc={doc} label={`${title}, light canvas`} options={options} canvas={lightCanvas} />
      <ViewerCard doc={doc} label={`${title}, dark canvas`} options={options} canvas={darkCanvas} />
    </div>
  );
}

function Frame({ title, kicker, children }: { title: string; kicker: ReactNode; children: ReactNode }) {
  return (
    <PrismPage title={title} caption={kicker}>
      <div className={prism.stack}>{children}</div>
    </PrismPage>
  );
}

function GalleryCards() {
  const protein = useMemo(() => alignmentGallery('SH2 family (case 3)', 'PF00017.seed.sto', PROTEIN), []);
  const rna = useMemo(() => alignmentGallery('tRNA family (case 4)', 'RF00005.seed.sto', NUCLEIC), []);
  const values = useMemo(valuesGallery, []);
  // The legend expanded (`UI.Legend`, brought back by `parts`; collapsed by default): the gallery shows each scheme's.
  const rows = { maxHeight: 132, legend: 'expanded', parts: LEGEND } as const;
  return (
    <>
      <LightAndDark doc={protein} title="Protein letters" options={rows} />
      <LightAndDark doc={rna} title="Nucleotide letters" options={rows} />
      <LightAndDark doc={values} title="Values and categories" options={{ legend: 'expanded', parts: LEGEND }} />
    </>
  );
}

export const Gallery: Story = {
  name: 'Gallery',
  render: () => (
    <Frame
      title="Every v1 color scheme, on a light and a dark canvas"
      kicker={
        <>
          Each section, track or layer names a scheme: <code>{'{"theme": "scheme", "params": {"name": "zappo"}}'}</code>
          , in the builder <code>b.colors.zappo</code>. The protein alignment is case 3's SH2 family, the nucleotide one
          case 4's tRNAs (each section follows the first: pan one, all move); values are a score in 0..1 and a change in
          −1..1 on p53, categories seven kinds of regions. On a dark skin a scheme's dark variant draws by itself: Blues
          runs dark to bright (the more, the brighter), Red–blue has a dark center, Percent identity its brighter blues;
          letters and categories keep their colors (the letter ink follows each cell). Clustal X and Percent identity
          read each column's composition. Click a legend group's title to change it in place.
        </>
      }
    >
      <GalleryCards />
    </Frame>
  ),
};

/** Command links, as a description's or the Stories agent's: `seqv:colorScheme` with the layer and a scheme. */
function Commands() {
  const { context } = useSeqViewer();
  const picked = useSeqState(context.ui.colorSchemes);
  const run = (layer: string, scheme: string | null) =>
    void context.commands.run({ command: 'seqv:colorScheme', params: { layer, scheme } }, 'host:links');
  return (
    <div className={prism.actions}>
      <div className={f.row}>
        <TextAction tone="quiet" onClick={() => run('members', 'taylor')}>
          seqv:colorScheme members → taylor
        </TextAction>
        <TextAction tone="quiet" onClick={() => run('members', 'story:pastel')}>
          → story:pastel
        </TextAction>
        <TextAction tone="quiet" onClick={() => run('cons.color', 'magma')}>
          cons.color → magma
        </TextAction>
        <TextAction tone="quiet" onClick={() => run('members', 'viridis')}>
          members → viridis (refused: letters are not values)
        </TextAction>
        <TextAction tone="quiet" onClick={() => [run('members', null), run('cons.color', null)]}>
          reset
        </TextAction>
      </div>
      <span className={f.code}>ui.colorSchemes: {JSON.stringify(Object.fromEntries(picked))}</span>
    </div>
  );
}

function PickerCard() {
  const doc = useMemo(pickerDoc, []);
  return (
    <ViewerCard
      doc={doc}
      label="A legend title, or a track's Colors ▸"
      options={{ maxHeight: 220, registry: HOST, legend: 'expanded', parts: [...LEGEND, 'TrackActions'] }}
    >
      <Commands />
    </ViewerCard>
  );
}

export const Picker: Story = {
  name: 'Live picker',
  render: () => (
    <Frame
      title="Switch and compare: a legend title, a track's menu, and the command"
      kicker={
        <>
          Click the <b>Conservation</b> title in the legend (dotted underline; the menu checks the scheme shown), or
          open the <b>Members</b> track's ⋯ menu and its <b>Colors ▸</b> (its letters list nothing, so it has no legend
          row): each lists "document" (the document's own color) and the schemes that fit that layer's color (letters of
          its alphabet, values, categories), colorblind-safe ones marked. A pick is UI state (
          <code>ui.colorSchemes</code>), never written into the document unless the host sets{' '}
          <code>UI.ColorSchemeWrites</code>; it survives patches and a new document clears it. The host registered{' '}
          <code>story:pastel</code>, listed after the built-ins. The links below run{' '}
          <code>{'{"command": "seqv:colorScheme", "params": {"layer": "members", "scheme": "taylor"}}'}</code>, as a
          description's command link or the Stories agent would.
        </>
      }
    >
      <PickerCard />
    </Frame>
  ),
};

export const ManyTracks: Story = {
  name: 'Many tracks',
  render: () => (
    <Frame
      title="A compact legend over 25 tracks"
      kicker={
        <>
          Five genes, ten transcripts and ten domains, as a template gives them. A track with nothing to list has no
          legend row (its colors are in its ⋯ menu, <b>Colors ▸</b>); tracks whose legends are identical share one
          group, titled by their common name with the count. Click a group's title to change the colors of every track
          in it; past eight groups the legend ends in "+n more". Hover a track's name: its tooltip shows that track's
          own legend. (The legend is expanded here; by default it is collapsed.)
        </>
      }
    >
      <ManyCards />
    </Frame>
  ),
};
function ManyCards() {
  const doc = useMemo(manyDoc, []);
  return (
    <LightAndDark
      doc={doc}
      title="25 tracks"
      options={{ maxHeight: 260, legend: 'expanded', parts: [...LEGEND, 'TrackActions'] }}
    />
  );
}

export const CustomPalette: Story = {
  name: 'Your own palette',
  render: () => (
    <Frame
      title="A document's own palettes: custom letters, a derived scheme, a custom ramp"
      kicker={
        <>
          <code>b.colors.custom</code> writes an inline scheme into the document: letter groups with labels (the
          legend's), a default for the alphabet's other letters, and a dark default (on the dark canvas the gray is
          dark).
          <code> b.colors.from(b.colors.zappo, {"{ overrides: { C: '#ff00cc' } }"})</code> recolors one letter of a
          named scheme. The ramp is three stops, its dark variant derived (light-low, so reversed). A bad color or a
          letter the alphabet lacks is an error that names what to write instead.
        </>
      }
    >
      <CustomCards />
    </Frame>
  ),
};
function CustomCards() {
  const doc = useMemo(customDoc, []);
  return <LightAndDark doc={doc} title="Custom" options={{ maxHeight: 160, legend: 'expanded', parts: LEGEND }} />;
}

function SharedCards() {
  const docs = useMemo(() => (['own', 'group', 'explicit'] as const).map(sharedDoc), []);
  const titles = ['Each layer its own domain', 'One domain_group', 'A group, and one explicit domain'];
  const options = { legend: 'expanded', parts: LEGEND } as const;
  return (
    <>
      {docs.map((doc, i) => (
        <LightAndDark key={titles[i]} doc={doc} title={titles[i]} options={options} />
      ))}
    </>
  );
}

export const SharedScales: Story = {
  name: 'Shared scales',
  render: () => (
    <Frame
      title="domain_group: one number, one color, across layers"
      kicker={
        <>
          A score in 0..1 and a count in 0..8 on p53, as bars and as Viridis color. First each layer has its own domain,
          so the score's middle (0.5) and the count's middle (4) are the same color and the bars both fill their tracks,
          with two ramps in the legend. With <code>domain_group: "signal"</code> on all four layers the domain is the
          union of their data (0..8): the score now stays in the dark end of one shared ramp, the bars share one y
          scale, and the legend lists the ramp once. An explicit domain wins in its layer: the score's color asks for
          0..16 (<code>b.colors.viridis({'{ domain: [0, 16] }'})</code>) and the count's takes the group's. Check the
          legend ramps and the same numbers' colors in light and dark.
        </>
      }
    >
      <SharedCards />
    </Frame>
  ),
};

export const ColorRange: Story = {
  name: 'Color range',
  render: () => (
    <Frame
      title="The range of a color scale: wider or narrower than the data"
      kicker={
        <>
          One signal in 0..8 on p53, as three color strips over Viridis: the data's own range, a range wider than the
          data (0..20: the strip stays in the dark end), and one narrower (3..5: values beyond the ends take the end
          colors). Each ramp in the legend has its range in two number fields, and a thin bar under the ramp marks where
          the data sits. Type a new end (a range may be larger than the data's) and the strip recolors; the reset link
          returns to the document's range. The same from a command:{' '}
          <code>{'{"command": "seqv:colorRange", "params": {"layer": "…", "range": [0, 20]}}'}</code>.
        </>
      }
    >
      <LightAndDark doc={useMemo(rangeDoc, [])} title="Color range" options={{ legend: 'expanded', parts: LEGEND }} />
    </Frame>
  ),
};
