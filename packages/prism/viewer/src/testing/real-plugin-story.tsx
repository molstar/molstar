/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { canvasInk, onCleanup } from '@molstar/prism-ui';
import { ShellStory, type ShellStoryProps } from '@molstar/prism-ui/testing/index';
import { createViewerPlugin } from '../sources/index.js';
import { mockViewerMode } from './viewer-mode.js';

// The viewer mode's mock over a real Mol* plugin (prism-sequence.md §3.1): the plugin is mounted on the shell's canvas,
// loads 1TQN from RCSB (the stories may use the network; tests never do), and its canvas props drive the shell's ink and
// the sequence model's skin. The canvas background is changed in the app's own Canvas settings (Background), which
// drives the plugin; it starts with nothing open.

/** The canvas backgrounds the story starts with (the `background` arg): light and dark among them. */
export const Backgrounds = {
  // White is the default: next to the glass panels (the UI theme's warm base) a cool gray reads as a clash, white does
  // not. The cool gray stays here for the comparison (the `background` arg).
  light: '#f3f5f8',
  white: '#ffffff',
  gray: '#8a9099',
  dark: '#111318',
  blue: '#2d3e50',
} as const;

export const Structure1TQN = { url: 'https://models.rcsb.org/1TQN.bcif', format: 'bcif', title: '1TQN' } as const;

/**
 * AlphaFill's model of CYP3A4 (AlphaFold's AF-P08684-F1 with ligands transplanted from PDB homologs): mmCIF (served as
 * text), pLDDT per residue. Several transplanted ligands overlap in space; the story shows the heme and erythromycin (a
 * large substrate in the same pocket), the other ligands are in the sequence but not drawn.
 */
export const StructureAlphaFill = {
  url: 'https://alphafill.eu/v1/aff/P08684',
  format: 'mmcif',
  title: 'CYP3A4 (AlphaFill)',
  ligands: ['HEM', 'ERY'],
} as const;

export function RealPluginShellStory(
  props: Omit<ShellStoryProps, 'loaded' | 'ink'> & {
    background?: string;
    structure?: typeof Structure1TQN | typeof StructureAlphaFill;
  },
) {
  const color = Backgrounds[props.background as keyof typeof Backgrounds] ?? props.background ?? Backgrounds.white;
  const host = (<div style="position: absolute; inset: 0" />) as HTMLElement;
  const plugin = createViewerPlugin(host, { background: color });
  const ink = canvasInk(plugin.canvasProps);
  void plugin.load(props.structure ?? Structure1TQN);
  onCleanup(() => {
    ink.dispose();
    plugin.dispose();
  });
  return (
    <ShellStory
      {...props}
      loaded={true}
      noPanel={true}
      ink={ink}
      mode={mockViewerMode({ plugin: { plugin, host, ink } })}
    />
  );
}
