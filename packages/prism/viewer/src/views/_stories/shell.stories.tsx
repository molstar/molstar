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
import { RealPluginShellStory, StructureAlphaFill } from '../../testing/real-plugin-story.js';
import { ViewerShellStory } from '../../testing/shell-story.js';

const meta: Meta = { title: 'Shell/Viewer mode', parameters: { layout: 'fullscreen' } };
export default meta;

/** The command menu open (⌘K / Ctrl+K anywhere, or the rail's ⌘): everything by name. */
export const CommandMenuOpen: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} commandMenu={true} />}</div>),
};

/** Nothing loaded: Home is the open tab. */
export const Empty: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={false} />}</div>),
};

/** A structure loaded: the Structure tab, a status message, a running task and the selection chip. */
export const Loaded: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} />}</div>),
};

/** Several structures: four CYP3A4 structures superposed on 1TQN, collapsed but the first. */
/**
 * A session with snapshots and a multi-model trajectory: playback appears at the top left only because the data has
 * something to play (plain 1TQN has neither).
 */
export const SessionWithPlayback: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} playback={{ snapshots: true, trajectory: true }} />
      </div>
    )),
};

export const ManyStructures: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} preset="many" />}</div>),
};

/** An EMDB entry: the map is a top-level volume next to the fitted model (EMD-2984 and 5A1A). */
export const EmdbMap: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} preset="em" />}</div>),
};

/**
 * The chrome over the canvas follows the canvas background automatically (light glass on a light background, dark on
 * a dark one), whatever the UI theme. Pick a background color in the controls.
 */
export const CanvasBackground: StoryObj<{ background: string }> = {
  args: { background: '#14181f' },
  argTypes: { background: { control: 'color' } },
  render: (args) =>
    story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} background={args.background} />}</div>),
};

/**
 * A real Mol* plugin on the canvas, loading 1TQN from RCSB (needs the network), under the same mock UI (structure rows,
 * selection bar and the rest are still mocks); it starts with no panel open. The sequence is the real Seq* viewer over
 * the plugin's structure, in a band of its own at the top of the canvas area (the 3D canvas starts below it), the
 * canvas' background behind it: hover and select residues in either view (they light up in both), the text version is
 * centered, scrolls to a focus and fades at its cut edges. The canvas stack's Sequence button (always there) opens a
 * menu: show or hide, Text / Tracks, the tracks' rail, the colors and the legend. The canvas background is Canvas
 * settings › Background (the stack's monitor button): the chrome, the band and the letters follow it; the `background`
 * arg only starts it. Hovering the 3D canvas shows its label next to the cursor, hovering the sequence the sequence's.
 */
export const RealPlugin: StoryObj<{ background: string }> = {
  args: { background: 'white' },
  argTypes: { background: { control: 'select', options: ['white', 'light', 'gray', 'dark', 'blue'] } },
  render: (args) =>
    story(() => (
      <div style="height: 100vh">
        <RealPluginShellStory background={args.background} />
      </div>
    )),
};

/**
 * The real plugin on a predicted structure with ligands (needs the network): AlphaFill's CYP3A4 model, the cartoon colored
 * by pLDDT (Mol*'s confidence theme, the AlphaFold bands) with the heme and erythromycin; in Tracks a pLDDT track under the
 * letters in the same bands, softened, and the legend (the Sequence menu) lists them. Everything else as the 1TQN story.
 */
export const RealPluginPredicted: StoryObj<{ background: string }> = {
  args: { background: 'white' },
  argTypes: { background: { control: 'select', options: ['white', 'light', 'gray', 'dark', 'blue'] } },
  render: (args) =>
    story(() => (
      <div style="height: 100vh">
        <RealPluginShellStory background={args.background} structure={StructureAlphaFill} />
      </div>
    )),
};

/** A command opened the Export panel at a section ("Export animation…" in the command menu). */
export const ExportAtSection: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} openExportAt="animation" />
      </div>
    )),
};

/**
 * The panel pinned (the pin in its header): it sits beside the canvas, which makes room, instead of floating over it;
 * for workflows that keep a panel open the whole time.
 */
export const Pinned: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} pinned={true} />}</div>),
};

/** Remote states turned on in the spec: a Remote states card next to Export, and the server's states on Home. */
export const RemoteStates: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} remoteStates={true} />}</div>),
};

/** A busy state server (50 states): the lists scroll on their own, on Home and in the Remote states card. */
export const ManyRemoteStates: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={false} remoteStates={50} />}</div>),
};

/** A narrow host: the panel becomes a bottom sheet with its tabs along the bottom. */
export const Compact: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh; display: grid; place-items: center; background: var(--msp-color-border)">
        <div style="width: 390px; height: 760px; border-radius: 16px; overflow: hidden; box-shadow: var(--msp-shadow-overlay)">
          <ViewerShellStory loaded={true} />
        </div>
      </div>
    )),
};

/**
 * Expand by context, in real frames of the Loaded story (hover the canvas's Expand button, at the top right):
 * - an iframe that allows fullscreen (`allow="fullscreen"`): the button is Fullscreen, since expanding would only fill
 *   the frame; Escape leaves it;
 * - an iframe that does not: Expand fills the frame, and Fullscreen beside it says why it is not possible;
 * - a page (open it below): Expand fills the page, Fullscreen beside it.
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

/**
 * Help, a tab at the rail's foot that stays open while you try things on the canvas (the rail's ?, the "Keyboard and mouse" command on ?, the command menu's ?): the canvas's mouse and
 * touch bindings (the mode's), the keys that work everywhere, and every command with a shortcut, from the registry as
 * it is; two commands on one key would be flagged first.
 */
export const Help: StoryObj = {
  render: () => story(() => <div style="height: 100vh">{<ViewerShellStory loaded={true} openTab="help" />}</div>),
};
