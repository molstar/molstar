/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SeqViewSpecDoc } from '@molstar/seq-core/format';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { type ReactNode, useMemo } from 'react';
import { sashimiDoc, scalesDoc, sizesDoc } from './support/seqstar-docs-sizes.js';
import f from './support/seqstar-format.module.css';
import { PrismPage } from './support/seqstar-prism.js';
import { ViewerCard, type ViewerOptions } from './support/seqstar-viewer.js';

/**
 * Size themes: a layer's size channels (a pair's width, apex height and its ends' radius, a feature's thickness, a
 * point's radius, a bar's thickness) take a theme as its color does: `uniform`, `rules`, or `field_scale`, a field
 * through a linear, sqrt or log scale into a range, or a string field through a map. Their legends: sample sizes drawn
 * as the channel draws them. The sashimi preset. The legend expanded; the theme is the toolbar's.
 */
const meta = { title: 'Seq*/Themes/Size themes', parameters: { layout: 'fullscreen' } } satisfies Meta;
export default meta;
type Story = StoryObj;

/** The legend (left out by `clean`) brought back and expanded (collapsed by default), so the size samples show. */
const EXPANDED: ViewerOptions = { legend: 'expanded', parts: ['Legend'] };
/** One viewer, in the toolbar's theme (light, dark, blue). */
function One({ doc, title }: { doc: SeqViewSpecDoc; title: string }) {
  return <ViewerCard doc={doc} label={title} options={EXPANDED} />;
}

function Frame({ title, kicker, children }: { title: string; kicker: ReactNode; children: ReactNode }) {
  return (
    <PrismPage title={title} caption={kicker}>
      <div className={f.stack}>{children}</div>
    </PrismPage>
  );
}

function Sizes() {
  const doc = useMemo(sizesDoc, []);
  return (
    <>
      <One doc={doc} title="Sized by a field" />
      <ViewerCard doc={doc} label="Sized by a field, the DOM UI" options={{ ...EXPANDED, ui: 'DOM' }} />
    </>
  );
}
function Scales() {
  const doc = useMemo(scalesDoc, []);
  return <One doc={doc} title="Scales" />;
}

export const BySize: Story = {
  name: 'Sized by a field',
  render: () => (
    <Frame
      title="Size themes: a field sets how large an element is"
      kicker={
        <>
          Contacts are stroked by how often they were seen (<code>width: field_scale</code> over <code>count</code>, 1
          to 8 px; the ends' radius too) and colored by the same count, so a heavy contact is thick and dark. Regions
          are as thick as their confidence (<code>thickness</code>, a share of the lane), points as large as their score
          (<code>radius</code>, px), and the bars as wide as their value (<code>thickness</code>, a share of the cell).
          A size changes how much of its slot an element fills, never the slot: hover and the layout do not move. The
          legend lists each scale's ends and middle, drawn as the channel draws them (a stroke, a bar, a dot), and in
          the color a ramp over the same field gives them (the contacts' widths in the count's blues). A matrix's height
          is the one channel not drawn: its bands are image tiles. It is reported and drawn at the viewer's size.
        </>
      }
    >
      <Sizes />
    </Frame>
  ),
};

export const ByScale: Story = {
  name: 'Linear, sqrt and log scales',
  render: () => (
    <Frame
      title="One score, three scales"
      kicker={
        <>
          The same twelve scores (1 to 100) as a point's radius, 2 to 8 px, through <code>scale</code>:{' '}
          <code>linear</code> (the large scores dominate), <code>sqrt</code> (area follows the value), and{' '}
          <code>log</code> (of one plus the distance from the domain's start: the small scores stay visible).
        </>
      }
    >
      <Scales />
    </Frame>
  ),
};

function SashimiCards() {
  const doc = useMemo(sashimiDoc, []);
  return <One doc={doc} title="Sashimi" />;
}

export const Sashimi: Story = {
  name: 'Sashimi',
  render: () => (
    <Frame
      title="Sashimi: splice junctions stroked by their read count"
      kicker={
        <>
          <code>b.preset.sashimi(section, junctions)</code> writes a <code>pairs</code> track per strand over a
          junctions table: rows whose <code>strand</code> is <code>+</code> rise from the line between the tracks,{' '}
          <code>-</code> hang below it (<code>side</code>, <code>filter</code>), each arc stroked by its{' '}
          <code>count</code> (<code>width: field_scale</code>, 1 to 8 px) and labeled with it under its apex where it
          fits (<code>labels: true</code>). Only generic params: the biology is the preset's. The exon-skipping junction
          (exon 1 to 3, 6 reads) is thin; the legend lists the widths for 2, 22 and 42 reads.
        </>
      }
    >
      <SashimiCards />
    </Frame>
  ),
};
