/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { story } from '@molstar/prism-ui/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import type { PrismViewerSpec } from '../../../spec.js';
import { ExampleFiles, RealViewerStory } from '../../../testing/real-viewer.js';
import { MoveSectionId, StructureMoveFeature } from '../index.js';

// Moving a structure by hand (structure-move.md): the Selection tab's Move section open, the transform handles on 1CBS.
// Drag an arrow, a square, a ring or the center; each release is one step (Undo puts it back); Reset, where it was
// loaded.

const meta: Meta = { title: 'Viewer/Move structures', parameters: { layout: 'fullscreen' } };
export default meta;

const withMove = (base: PrismViewerSpec): PrismViewerSpec => ({
  ...base,
  extensions: [...(base.extensions ?? []), StructureMoveFeature],
});

/** 1CBS with the Move section open: the handles on it */
export const MoveSection: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={withMove}
        load={[ExampleFiles.cbs]}
        setup={async (v) => {
          await v.api.view.open({ tab: 'selection', section: MoveSectionId });
        }}
      />
    )),
};
