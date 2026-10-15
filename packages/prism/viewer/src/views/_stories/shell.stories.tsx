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
import { ExampleFiles, RealViewerStory } from '../../testing/real-viewer.js';
import type { PrismViewer } from '../../viewer.js';

// The viewer mode: the real viewer (`createPrismViewer`, the viewer's default spec) over the repository's example
// files. A few stories reach the network (an EMDB map, a state server): said on them.

const meta: Meta = { title: 'Shell/Viewer mode', parameters: { layout: 'fullscreen' } };
export default meta;

/** Clicks a button of the shell by its accessible name (the shell is in a shadow root). */
function click(v: PrismViewer, label: string) {
  const root = [...document.querySelectorAll('*')].find((e) => e.shadowRoot)?.shadowRoot;
  void v;
  root?.querySelector<HTMLButtonElement>(`button[aria-label="${label}"]`)?.click();
}

/** Nothing loaded: Home is open (`ui.homeWhenEmpty`), where loading starts. */
export const Empty: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory spec={(base) => ({ ...base, options: { ...base.options, ui: { homeWhenEmpty: true } } })} />
    )),
};

/** A structure loaded (1CBS): the sequence band, the canvas, the Structure tab. */
export const Loaded: StoryObj<{ compact: boolean }> = {
  args: { compact: false },
  render: (args) =>
    story(() => (
      <div style="height: 100vh; display: grid; place-items: center; background: var(--msp-color-border)">
        <RealViewerStory
          load={[ExampleFiles.cbs]}
          style={args.compact ? 'width: 390px; height: 760px' : 'width: 100vw; height: 100vh'}
          setup={(v) => v.views.shell?.show('structure')}
        />
      </div>
    )),
};

/** The command menu open (⌘K / Ctrl+K anywhere, or the rail's ⌘): everything by name. */
export const CommandMenuOpen: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => click(v, 'Command menu')} />),
};

/** Several structures: 1CBS, crambin and gramicidin, each its row in the tree. */
export const ManyStructures: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs, ExampleFiles.crambin, ExampleFiles.grm]}
        setup={(v) => v.views.shell?.show('structure')}
      />
    )),
};

/** Playback: a trajectory's frames (the trajectory bar) and two snapshots (the snapshots card). */
export const Playback: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.trajectory?.({
            model: ExampleFiles.trajectoryModel,
            modelFormat: 'pdb',
            modelBinary: false,
            coordinates: ExampleFiles.trajectoryCoordinates,
            coordinatesFormat: 'dcd',
            coordinatesBinary: true,
          });
          await v.api.snapshot?.add({ name: 'The whole protein' });
          await v.api.snapshot?.add({ name: 'Again, for the playback' });
        }}
      />
    )),
};

/**
 * An EMDB entry and its fitted model (EMD-2984 and 5A1A): the map is a volume next to the structure, picked: its levels
 * (the histogram with the iso value to drag, the color, the opacity) in the inspector. Network.
 */
export const EmdbMap: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'pdb', id: '5A1A' });
          await v.api.load.entry?.({ source: 'emdb', id: 'EMD-2984' });
          v.views.shell?.show('structure');
          // The map picked: its levels (the histogram, the iso value, the color) in the inspector.
          v.views.scene.pick(v.views.scene.roots.value.find((n) => n.kind === 'volume'));
        }}
      />
    )),
};

/**
 * The map as a direct volume (EMD-2984), its row picked: the transfer function in the inspector (the opacity curve over
 * the map's histogram, the colors under it, the color stops). Network.
 */
export const EmdbDirectVolume: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        setup={async (v) => {
          await v.api.load.entry?.({ source: 'emdb', id: 'EMD-2984' });
          await v.api.volume.representation({ type: 'direct-volume' });
          v.views.shell?.show('structure');
          const map = v.views.scene.roots.value.find((n) => n.kind === 'volume');
          v.views.scene.pick(map?.children.value.find((n) => n.label.value.toLowerCase().includes('direct')) ?? map);
        }}
      />
    )),
};

/** The canvas background: the canvas' chrome (the toolbar, the band, the status) follows its color. */
export const CanvasBackground: StoryObj<{ background: string }> = {
  args: { background: '#14181f' },
  argTypes: { background: { control: 'color' } },
  render: (args) =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={(v) => void v.api.canvas.background({ color: args.background })}
      />
    )),
};

/** A command opened the Export card at a section (Session & state). */
export const ExportAtSection: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={(v) => void v.api.view.open({ card: 'core.overlay.export', section: 'export-session.section.session' })}
      />
    )),
};

/** The panel pinned: the canvas gives it room instead of running under it. */
export const Pinned: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={(v) => {
          v.views.shell?.show('structure');
          v.views.shell?.setPinned('right', true);
        }}
      />
    )),
};

/** Remote states on in the spec: a Remote states card next to Export, the server's states on Home. Network. */
export const RemoteStates: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={(base) => ({
          ...base,
          options: { ...base.options, ui: { showRemoteStates: true, homeWhenEmpty: true } },
        })}
      />
    )),
};

/** A narrow host: the panel becomes a bottom sheet with its tabs along the bottom. */
export const Compact: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh; display: grid; place-items: center; background: var(--msp-color-border)">
        <div style="width: 390px; height: 760px; border-radius: 16px; overflow: hidden; box-shadow: var(--msp-shadow-overlay)">
          <RealViewerStory load={[ExampleFiles.cbs]} style="height: 100%" />
        </div>
      </div>
    )),
};

/**
 * Expand or fullscreen, by context: in a frame that allows fullscreen the button is Fullscreen; in one that does not,
 * Expand (the canvas fills the viewer).
 */
export const ExpandByContext: StoryObj = {
  parameters: { layout: 'fullscreen' },
  render: () =>
    story(() => {
      const src = 'iframe.html?id=prism-shell-viewer-mode--loaded&viewMode=story';
      const caption = 'color: var(--msp-color-text-muted); font-size: var(--msp-font-size-sm);';
      const frame =
        'width: 100%; height: 100%; border: 0; border-radius: var(--msp-radius-lg); box-shadow: inset 0 0 0 1px var(--msp-glass-edge);';
      return (
        <div style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-3); height: 100vh; padding: var(--msp-space-4); box-sizing: border-box; background: var(--msp-color-bg);">
          <div style={`${caption} display: flex; gap: var(--msp-space-4);`}>
            <span>Hover Expand at the top right of each canvas.</span>
            <a href={src} target="_blank" rel="noopener" style="color: var(--msp-color-accent)">
              Open as a page (Expand, Fullscreen beside it)
            </a>
          </div>
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: var(--msp-space-4); min-height: 0;">
            <figure style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-2); margin: 0; min-height: 0;">
              <figcaption style={caption}>In an iframe that allows fullscreen: the button is Fullscreen</figcaption>
              <iframe title="Allows fullscreen" src={src} allow="fullscreen" style={frame} />
            </figure>
            <figure style="display: grid; grid-template-rows: auto 1fr; gap: var(--msp-space-2); margin: 0; min-height: 0;">
              <figcaption style={caption}>In an iframe that does not: Expand, Fullscreen not possible</figcaption>
              {/* Same-origin frames may go fullscreen by default ('self'): forbidden explicitly, as a page on another origin is. */}
              <iframe title="Does not allow fullscreen" src={src} allow="fullscreen 'none'" style={frame} />
            </figure>
          </div>
        </div>
      );
    }),
};

/** Help: the keys and the mouse, and the viewer's words. */
export const Help: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => v.views.shell?.show('help')} />),
};
