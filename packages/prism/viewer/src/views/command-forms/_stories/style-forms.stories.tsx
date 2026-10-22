/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import type { SceneNode } from '../../../model/scene.js';
import type { PrismViewer } from '../../../viewer.js';
import { ExampleFiles, findNode, openAt, RealViewerStory } from '../../../testing/real-viewer.js';

// The command menu's own forms for a component's look (commands.md §5.1) in the real viewer over 1CBS: opened from the
// Structure tab's component inspector (its Add representation… and Set coloring… rows), and from the command menu.

const meta: Meta = { title: 'Viewer/Command menu/Own forms', parameters: { layout: 'fullscreen' } };
export default meta;

const isPolymer = (n: SceneNode) => n.kind === 'component' && n.label.value === 'Polymer';

/** The Polymer component picked, then one of its inspector rows run (as a click on it). */
function fromInspector(row: 'add-repr' | 'coloring') {
  return (v: PrismViewer) => {
    openAt(v, isPolymer);
    findNode(v.views.scene.roots.value, isPolymer)
      ?.menu?.find((i) => i.id === row)
      ?.run?.();
  };
}

/** Add representation… from the inspector: the types Polymer supports (Cartoon marked shown), one it does not show picked, its options under them. */
export const ShowAs: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={fromInspector('add-repr')} />),
};

/** Set coloring… from the inspector: the themes that apply, painted with their colors on 1CBS. */
export const ColorByTheme: StoryObj = {
  render: () => story(() => <RealViewerStory load={[ExampleFiles.cbs]} setup={fromInspector('coloring')} />),
};

/** The uniform theme picked: its color's swatches (the whole picker behind the chevron). */
export const ColorByThemeUniform: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        load={[ExampleFiles.cbs]}
        setup={(v) => {
          fromInspector('coloring')(v);
          const form = v.managers.commands.active.value;
          form?.set('theme', 'uniform', { stay: true });
          form?.setParam('color', '#e15759');
        }}
      />
    )),
};

/** Show as… from the command menu, nothing named: the types the structure supports, the component to name in the sentence. */
export const ShowAsFromCommandMenu: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory load={[ExampleFiles.cbs]} setup={(v) => void v.managers.commands.open('representation.set')} />
    )),
};
