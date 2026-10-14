/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seq-viewer';
import { defineSeqViewer, type SeqViewerElement } from '@molstar/seq-viewer/dom';
import { InteractionsViewerRegistry } from '@molstar/seq-viewer/interactions/viewer';
import type { Meta, StoryObj } from '@storybook/react-vite';
import { useCallback, useEffect, useRef } from 'react';
import { fixtureLoaders, MolstarFixture } from '../../../core/_test/format/load-fixtures.js';
import { ComposedDom } from './support/seqstar-composed.js';
import { MarkdownDescription } from './support/seqstar-description.js';
import { describedCase10 } from './support/seqstar-docs-views.js';
import f from './support/seqstar-format.module.css';
import { PrismCanvas, PrismPage } from './support/seqstar-prism.js';
import { STORY_CHROME } from './support/seqstar-viewer.js';

/**
 * The viewer's DOM UI: the `<seq-viewer>` element on the canvas
 * surface, with its skin, parts, menus, tooltip and legend. The viewer has no markdown: a `markdown` description is shown
 * by a host's part through the `Description` slot, here case 10's, rendered with unified and `@molstar/seq-core/links`.
 */
interface Args {
  /** `clean`: the canvas with hover, selection and the tooltip; no tabs, header, menus or legend. */
  chrome: 'full' | 'clean';
}
const meta = {
  title: 'Seq*/Cases/Default UIs',
  parameters: { layout: 'fullscreen' },
  args: { chrome: STORY_CHROME },
  argTypes: {
    chrome: { control: 'inline-radio', options: ['full', 'clean'], description: "This UI's parts (a UI option)" },
  },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

/** Case 10 in a `<seq-viewer>`, its description by the host's markdown part. */
function Case10Links({ chrome }: Args) {
  const host = useRef<HTMLDivElement>(null);
  const el = useRef<SeqViewerElement | null>(null);
  // The element once per mount; its properties follow the args (`chrome` remounts the UI on the same context).
  useEffect(() => {
    defineSeqViewer();
    const e = document.createElement('seq-viewer');
    const base = DefaultSeqViewerSpec();
    e.spec = { ...base, registry: [...(base.registry ?? []), { svs: MolstarFixture }, InteractionsViewerRegistry] };
    e.ports = { loaders: fixtureLoaders() };
    e.ui = { slots: { Description: MarkdownDescription } };
    e.document = describedCase10();
    host.current?.append(e);
    el.current = e;
    return () => {
      e.remove();
      el.current = null;
    };
  }, []);
  useEffect(() => {
    if (el.current) el.current.chrome = chrome;
  }, [chrome]);
  return (
    <PrismPage
      title="Case 10: a description with links"
      caption={
        <>
          Case 10's description is markdown, which the viewer does not render: a host part through the Description slot
          does, with <code>@molstar/seq-core/links</code>. Hover a link to light its target, click to select or focus,
          step the snapshots: their descriptions replace this one and release its links.
        </>
      }
    >
      <section className={f.card}>
        <h2 className={f.heading}>
          {'<seq-viewer>'}, chrome {chrome}
        </h2>
        <PrismCanvas>
          <div ref={host} />
        </PrismCanvas>
      </section>
    </PrismPage>
  );
}
export const DescriptionLinks: Story = {
  name: 'DOM UI, description links',
  render: (a) => <Case10Links {...a} />,
};

/** Case 10 in a `clean` `<seq-viewer>`, its tabs, export, issues and legend mounted by the host around it. */
function ComposedElement() {
  const make = useCallback((e: SeqViewerElement) => {
    const base = DefaultSeqViewerSpec();
    e.spec = { ...base, registry: [...(base.registry ?? []), { svs: MolstarFixture }, InteractionsViewerRegistry] };
    e.ports = { loaders: fixtureLoaders() };
    e.document = describedCase10();
  }, []);
  return (
    <PrismPage
      title="Case 10: the chrome placed by the host"
      caption={
        <>
          The DOM parts are functions a host mounts into its own elements: here a <code>clean</code> element, and{' '}
          <code>ViewTabs</code>, <code>Export</code> and <code>IssueBadge</code> above it and the <code>Legend</code>{' '}
          below it, each given the element&apos;s context and its <code>mountedUI</code>.
        </>
      }
    >
      <section className={f.card}>
        <h2 className={f.heading}>{'<seq-viewer chrome=clean>'} and the host&apos;s parts</h2>
        <ComposedDom make={make} />
      </section>
    </PrismPage>
  );
}
export const DomComposed: Story = {
  name: 'DOM UI, the chrome placed by the host',
  render: () => <ComposedElement />,
};
