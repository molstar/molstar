/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { PlaybackCard } from '../features/snapshots/model.js';
import { ExampleFiles, RealViewerStory } from '../testing/real-viewer.js';
import type { PrismViewerSpec } from '../spec.js';
import type { PrismViewer } from '../viewer.js';
import { DpdbFixtures } from './dynamic-pdb/data.js';
import { DpdbCard, dynamicPdbExtension } from './dynamic-pdb/extension.js';

// Every major layout of the shell, one story each, in a frame of a fixed size (the shell's mode follows its own width,
// not the window's: shell.md §1), to see at a glance where the rail, the panel, the cards, the canvas toolbar and the
// rail's button sit. Desktop (floating, pinned, docked), the compact rail (`compactRail`), embedded, expanded, and the
// narrow layout (a phone: the panel and the cards are bottom sheets). The real viewer over 1CBS.

const meta: Meta = { title: 'Viewer/Layouts', parameters: { layout: 'fullscreen' } };
export default meta;

type Size = 'desktop' | 'laptop' | 'tablet' | 'phone';
const Sizes: Record<Size, { readonly width: number; readonly height: number }> = {
  desktop: { width: 1280, height: 800 },
  laptop: { width: 980, height: 640 },
  tablet: { width: 700, height: 900 },
  phone: { width: 390, height: 760 },
};

type Layout = NonNullable<NonNullable<PrismViewerSpec['ui']>['layout']>;

interface FrameProps {
  readonly size: Size;
  readonly layout?: Layout;
  readonly setup?: (viewer: PrismViewer) => void | Promise<void>;
  readonly dpdb?: boolean;
}

/** The viewer in a frame of a device's size (at most), centred on the page; its layout options and a setup after the load. */
function Frame(props: FrameProps) {
  const { width, height } = Sizes[props.size];
  const spec = (base: PrismViewerSpec): PrismViewerSpec => {
    let s = base;
    if (props.layout) s = { ...s, ui: { ...s.ui, layout: { ...s.ui?.layout, ...props.layout } } } as PrismViewerSpec;
    if (props.dpdb) s = { ...s, extensions: [...(s.extensions ?? []), dynamicPdbExtension(DpdbFixtures)] };
    return s;
  };
  return (
    // The page's height, never more: the frame is the device's size at most and shrinks to fit (no scrolling to see it
    // whole). Its mode still follows its own width (a desktop frame stays desktop above `compactWidth`).
    <div style="height: 100vh; display: grid; grid-template: minmax(0, 1fr) / minmax(0, 1fr); place-items: center; padding: var(--msp-space-4); box-sizing: border-box; overflow: hidden; background: var(--msp-color-border)">
      <div
        style={`width: min(${width}px, 100%); height: min(${height}px, 100%); border-radius: 12px; overflow: hidden; box-shadow: var(--msp-shadow-overlay)`}
      >
        {RealViewerStory({
          ...(props.dpdb ? {} : { load: [ExampleFiles.cbs] }),
          spec,
          style: 'height: 100%',
          ...(props.setup ? { setup: props.setup } : {}),
        })}
      </div>
    </div>
  );
}

const frame = (props: FrameProps): StoryObj => ({ render: () => story(() => Frame(props)) });

/** Three snapshots, so the Story card has something to show and step through. */
async function withStory(v: PrismViewer) {
  const person = v.api.as({ kind: 'person' });
  const steps = [
    [
      'The fold',
      'Cellular retinoic acid-binding protein II (**1CBS**): a ten-stranded β-barrel with a helix-turn-helix lid.',
    ],
    ['The ligand', 'Retinoic acid sits inside the barrel, its carboxylate toward the polar residues at the bottom.'],
    ['The lid', 'The two helices close the barrel; they move to let the ligand in.'],
  ];
  for (const [name, description] of steps) await person.run('snapshot.add', { name, description });
  await person.view.open({ card: PlaybackCard });
}

const showPanel = (tab: string) => (v: PrismViewer) => void v.views.shell?.show(tab);
const unpinned = (v: PrismViewer) => v.views.shell?.setPinned('right', false);

async function loadDpdb(v: PrismViewer) {
  const person = v.api.as({ kind: 'person' });
  const r = await person.run('dynamic-pdb.load', { entry: 'dpdb_gijqn2no' });
  if (!r.ok) throw new Error(r.message);
  await person.view.open({ card: DpdbCard });
}

// Desktop: the rail on the right, the panel floating (pinned by default: the canvas ends at it) or docked.

export const DesktopRail: StoryObj = { ...frame({ size: 'desktop' }), name: 'Desktop: rail, no panel' };
export const DesktopPinned: StoryObj = {
  ...frame({ size: 'desktop', setup: showPanel('structure') }),
  name: 'Desktop: panel pinned',
};
export const DesktopFloating: StoryObj = {
  ...frame({
    size: 'desktop',
    setup: (v) => {
      unpinned(v);
      showPanel('structure')(v);
    },
  }),
  name: 'Desktop: panel floating',
};
export const DesktopDocked: StoryObj = {
  ...frame({ size: 'desktop', layout: { docks: { right: { floating: false } } }, setup: showPanel('structure') }),
  name: 'Desktop: panel docked',
};
export const DesktopStoryPinned: StoryObj = {
  ...frame({ size: 'desktop', setup: withStory }),
  name: 'Desktop: Story card, pinned',
};
export const DesktopStoryFloating: StoryObj = {
  ...frame({
    size: 'desktop',
    setup: async (v) => {
      unpinned(v);
      await withStory(v);
    },
  }),
  name: 'Desktop: Story card, floating',
};
export const DesktopDpdb: StoryObj = {
  ...frame({ size: 'desktop', dpdb: true, setup: loadDpdb }),
  name: 'Desktop: Dynamic PDB card',
};
export const DesktopExpanded: StoryObj = {
  ...frame({ size: 'desktop', setup: (v) => void v.views.shell?.setExpanded(true) }),
  name: 'Desktop: expanded',
};

// The compact rail: the rail's column hidden behind the button in the canvas's corner (`compactRail`).

export const CompactRailHidden: StoryObj = {
  ...frame({ size: 'laptop', layout: { compactRail: true } }),
  name: 'Compact rail: hidden',
};
export const CompactRailShown: StoryObj = {
  ...frame({ size: 'laptop', layout: { compactRail: true }, setup: (v) => v.views.shell?.toggleRail() }),
  name: 'Compact rail: shown',
};
export const CompactRailStory: StoryObj = {
  ...frame({ size: 'laptop', layout: { compactRail: true }, setup: withStory }),
  name: 'Compact rail: Story card',
};
export const CompactRailStoryRailShown: StoryObj = {
  ...frame({
    size: 'laptop',
    layout: { compactRail: true },
    setup: async (v) => {
      await withStory(v);
      v.views.shell?.toggleRail();
    },
  }),
  name: 'Compact rail: Story card, rail shown',
};

// Embedded: a viewer in a page (the compact rail, the panel floating, never pinned).

export const Embedded: StoryObj = { ...frame({ size: 'laptop', layout: { mode: 'embedded' } }), name: 'Embedded' };
export const EmbeddedStory: StoryObj = {
  ...frame({ size: 'laptop', layout: { mode: 'embedded' }, setup: withStory }),
  name: 'Embedded: Story card',
};

// Narrow (compact mode, below `compactWidth`): the panel and the cards are bottom sheets, the tabs along the bottom.

export const Tablet: StoryObj = { ...frame({ size: 'tablet' }), name: 'Tablet' };
export const Phone: StoryObj = { ...frame({ size: 'phone' }), name: 'Phone' };
export const PhonePanel: StoryObj = {
  ...frame({ size: 'phone', setup: showPanel('structure') }),
  name: 'Phone: panel',
};
export const PhoneStory: StoryObj = { ...frame({ size: 'phone', setup: withStory }), name: 'Phone: Story card' };
export const PhoneDpdb: StoryObj = {
  ...frame({ size: 'phone', dpdb: true, setup: loadDpdb }),
  name: 'Phone: Dynamic PDB card',
};
