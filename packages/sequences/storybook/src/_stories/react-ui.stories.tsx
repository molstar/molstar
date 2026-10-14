/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultSeqViewerSpec } from '@molstar/seq-viewer';
import { InteractionsViewerRegistry } from '@molstar/seq-viewer/interactions/viewer';
import { SeqViewer } from '@molstar/seq-viewer/react';
import type { Meta, StoryObj } from './support/seqstar-story.js';
import { useMemo } from 'react';
import { case10 } from '../../../core/_test/format/cases.js';
import { fixtureLoaders, MolstarFixture } from '../../../core/_test/format/load-fixtures.js';
import '@molstar/seq-viewer/skin.css';
import { ComposedReact } from './support/seqstar-composed.js';
import f from './support/seqstar-format.module.css';
import { CANVAS_COLORS, type CanvasColor, PrismCanvas, PrismPage } from './support/seqstar-prism.js';
import { STORY_CHROME } from './support/seqstar-viewer.js';

/**
 * The viewer's React UI: `SeqViewer` from `/react` on the canvas
 * surface, with the skin, parts, menus, tooltip and legend, the same as the DOM UI's element. The `ui` is an inline
 * object each render: the provider compares it in two levels, so a toggle of `chrome` changes the parts and keeps the
 * surfaces.
 */
interface Args {
  /** `clean`: the canvas with hover, selection and the tooltip; no tabs, header, menus or legend. */
  chrome: 'full' | 'clean';
  /** The canvas's color: the UI theme's base, or a color of its own with the theme it implies (seqstar-prism.tsx). */
  canvas: CanvasColor;
}
const meta = {
  title: 'Seq*/Cases/Default UIs',
  parameters: { layout: 'fullscreen' },
  args: { chrome: STORY_CHROME, canvas: 'UI theme' },
  argTypes: {
    chrome: { control: 'inline-radio', options: ['full', 'clean'], description: "This UI's parts (a UI option)" },
    canvas: {
      control: 'inline-radio',
      options: Object.keys(CANVAS_COLORS),
      description: "The canvas's color: its chrome follows it, not the toolbar's UI theme",
    },
  },
} satisfies Meta<Args>;
export default meta;
type Story = StoryObj<Args>;

// The fixtures' structures and the interactions renderers, as the DOM UI's story; read at mount.
const base = DefaultSeqViewerSpec();
const spec = { ...base, registry: [...(base.registry ?? []), { svs: MolstarFixture }, InteractionsViewerRegistry] };
const ports = { loaders: fixtureLoaders() };

/** Case 10: the Stories drafts, with its snapshots. */
function Case10View({ chrome, canvas }: Args) {
  const doc = useMemo(case10, []);
  const color = CANVAS_COLORS[canvas] ?? undefined;
  return (
    <PrismPage
      title="Case 10: The Stories drafts"
      caption={
        <>
          The default React UI on case 10. <code>clean</code> chrome (the default) is the canvas with hover, selection
          and the tooltip; <code>full</code> adds tabs, section headers and their switches, a track&apos;s menu at the
          hovered header (or the keyboard cursor&apos;s, Shift+F10), the legend and the export menu. Toggle{' '}
          <code>chrome</code>: the parts change, the surfaces stay. The viewer draws on a canvas whose color and theme
          may differ from the UI&apos;s (<code>canvas</code>; the UI theme is the toolbar&apos;s): its chrome follows
          the canvas, as the chrome over the 3D canvas does in Mol* Prism.
        </>
      }
    >
      <section className={f.card}>
        <h2 className={f.heading}>
          SeqViewer, chrome {chrome}, canvas {color ? `${canvas} (${color})` : 'in the UI theme'}
        </h2>
        <PrismCanvas key={canvas} color={color}>
          <SeqViewer spec={spec} ports={ports} document={doc} ui={{ chrome }} />
        </PrismCanvas>
      </section>
    </PrismPage>
  );
}
export const Case10: Story = {
  name: 'React UI, with the chrome toggle',
  render: (a) => <Case10View {...a} />,
};

/** Case 10 with the chrome placed by the host: the parts are components, the viewer `clean`. */
function ComposedView({ canvas }: Args) {
  const doc = useMemo(case10, []);
  const color = CANVAS_COLORS[canvas] ?? undefined;
  return (
    <PrismPage
      title="Case 10: the chrome placed by the host"
      caption={
        <>
          Every piece of the viewer&apos;s chrome is a component a host places itself, against the viewer&apos;s
          context: here the viewer is <code>clean</code> (the canvas, hover and the tooltip) and the host puts{' '}
          <code>ViewTabs</code>, each section&apos;s <code>SectionSwitches</code>, <code>Export</code> and{' '}
          <code>IssueBadge</code> in a header row of its own, and the <code>Legend</code> below the viewer.{' '}
          <code>chrome: &quot;full&quot;</code> is the same parts in the viewer&apos;s own composition.
        </>
      }
    >
      <section className={f.card}>
        <h2 className={f.heading}>
          A host&apos;s composition, canvas {color ? `${canvas} (${color})` : 'in the UI theme'}
        </h2>
        <ComposedReact key={canvas} doc={doc} color={color} />
      </section>
    </PrismPage>
  );
}
export const Composed: Story = {
  name: 'React UI, the chrome placed by the host',
  render: (a) => <ComposedView {...a} />,
};
