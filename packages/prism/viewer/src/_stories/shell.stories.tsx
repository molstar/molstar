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
import { ExampleFiles, RealViewerStory } from '../testing/real-viewer.js';
import type { PrismViewer } from '../viewer.js';

// The whole app: the real viewer (`createPrismViewer`, the viewer's default spec) over the repository's example files,
// for what is about the shell (its layouts, the command menu, the canvas's chrome, cards and pages). What a panel or a
// kind of data shows has its own stories beside (Viewer › Home, Structure panel, Volumes, Playback, Data).

const meta: Meta = { title: 'Viewer/Shell', parameters: { layout: 'fullscreen' } };
export default meta;

/** Clicks a button of the shell by its accessible name (the shell is in a shadow root). */
function click(v: PrismViewer, label: string) {
  const root = [...document.querySelectorAll('*')].find((e) => e.shadowRoot)?.shadowRoot;
  void v;
  root?.querySelector<HTMLButtonElement>(`button[aria-label="${label}"]`)?.click();
}

/** A structure loaded (1CBS): the sequence band, the canvas, the Structure panel. */
export const Loaded: StoryObj = {
  render: () =>
    story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => v.views.shell?.show('structure')} />),
};

/** The command menu open (⌘K / Ctrl+K anywhere, or the rail's ⌘): everything by name. */
/** The hover label beside the pointer (`options.ui.hoverLabel: 'cursor'`); the default is the canvas' bottom right. */
export const HoverLabelAtCursor: StoryObj = {
  name: 'Hover label at the cursor',
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        spec={(base) => ({ ...base, options: { ...base.options, ui: { ...base.options?.ui, hoverLabel: 'cursor' } } })}
      />
    )),
};

export const CommandMenuOpen: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => click(v, 'Command menu')} />),
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
      const src = 'iframe.html?id=prism-viewer-shell--loaded&viewMode=story';
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
