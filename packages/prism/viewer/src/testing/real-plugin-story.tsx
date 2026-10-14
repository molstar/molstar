/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canvasInk, Choices, observableState, onCleanup } from '@molstar/prism-ui';
import { ShellStory, type ShellStoryProps } from '@molstar/prism-ui/testing/index';
import { createViewerPlugin } from '../sources/index.js';
import { mockViewerMode } from './viewer-mode.js';

// The viewer mode's mock over a real Mol* plugin (prism-sequence.md §3.1): the plugin is mounted on the shell's canvas,
// loads 1TQN from RCSB (the stories may use the network; tests never do), and its canvas props drive the shell's ink and
// the sequence model's skin. A small control above the app changes the canvas background.

/** The canvas backgrounds the story offers: light and dark among them, so the ink and the skin follow both ways. */
export const Backgrounds = [
  { id: 'light', label: 'Light', title: 'Prism light base', color: '#eeece7' },
  { id: 'white', label: 'White', title: 'White', color: '#ffffff' },
  { id: 'gray', label: 'Gray', title: 'Mid gray', color: '#8a9099' },
  { id: 'dark', label: 'Dark', title: 'Prism dark base', color: '#111318' },
  { id: 'blue', label: 'Blue', title: 'Prism blue base', color: '#2d3e50' },
] as const;

export const Structure1TQN = { url: 'https://models.rcsb.org/1TQN.bcif', format: 'bcif', title: '1TQN' } as const;

export function RealPluginShellStory(props: Omit<ShellStoryProps, 'loaded' | 'ink'> & { background?: string }) {
  const initial = Backgrounds.find((b) => b.id === props.background) ?? Backgrounds[0];
  const host = (<div style="position: absolute; inset: 0" />) as HTMLElement;
  const plugin = createViewerPlugin(host, { background: initial.color });
  const ink = canvasInk(plugin.canvasProps);
  const active = observableState<string>(initial.id);
  void plugin.load(Structure1TQN);
  onCleanup(() => {
    ink.dispose();
    plugin.dispose();
  });
  return (
    <div style="display: grid; grid-template-rows: auto 1fr; height: 100%; background: var(--msp-color-bg); color: var(--msp-color-text); font-size: var(--msp-font-size-sm)">
      <div style="display: flex; align-items: center; gap: var(--msp-space-4); padding: var(--msp-space-2) var(--msp-space-4)">
        <span style="color: var(--msp-color-text-muted)">Canvas background</span>
        <div style="width: 420px">
          <Choices
            choices={Backgrounds}
            active={active.asReadonly}
            onPick={(id) => {
              active.set(id);
              plugin.setBackground(Backgrounds.find((b) => b.id === id)?.color ?? initial.color);
            }}
          />
        </div>
        <span style="color: var(--msp-color-danger)">{plugin.error}</span>
      </div>
      <div style="min-height: 0">
        <ShellStory {...props} loaded={true} ink={ink} mode={mockViewerMode({ plugin: { plugin, host, ink } })} />
      </div>
    </div>
  );
}
