/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { FrameCounts, onCleanup, restartable, Show } from '../../../dom/index.js';
import { state } from '../../../state.js';
import { Button } from '../../button/button.js';
import { Stage } from '../../../stage/stage.js';
import { StageModel } from '../../../stage/stage-model.js';
import { StagePlaceholder } from '../../../stage/placeholders.js';
import { GlassSplash } from '../splash.js';

const meta: Meta = { title: 'Components/Splash' };
export default meta;

/** The stage empty: after its wait the glass assembles; press it, it shatters and gathers again. Load and Clear. */
export const Glass: StoryObj = {
  render: () =>
    story(() => {
      const stage = new StageModel({ background: '#f4f5f7' });
      onCleanup(() => stage.dispose());
      return (
        <div style="position: relative; display: flex; height: 100vh">
          <Stage
            model={stage}
            canvas={<StagePlaceholder />}
            lead={
              <div style="display: flex; gap: var(--msui-space-2)">
                <Button size="sm" onClick={() => stage.setEmpty(false)}>
                  Load
                </Button>
                <Button size="sm" onClick={() => stage.setEmpty(true)}>
                  Clear
                </Button>
              </div>
            }
          />
        </div>
      );
    }),
};

/** At rest: once assembled it asks for no frames (the count stops); hover it and it moves, then rests again. */
export const Idle: StoryObj = {
  render: () =>
    story(() => {
      const shown = state(true);
      const counted = state(0);
      // Read once a second (a timer, not a frame: reading must not keep frames coming).
      const read = restartable(() => {
        counted.set(FrameCounts.requested);
        read.start(1000);
      });
      read.start(0);
      return (
        <div style="display: grid; gap: var(--msui-space-4); justify-items: center; padding: calc(24 * var(--msui-unit)) var(--msui-space-5)">
          <Show when={shown}>{() => <GlassSplash />}</Show>
          <span style="color: var(--msui-text-muted)" class="msui-num">
            frames asked for on the page: {counted}
          </span>
          <Button size="sm" onClick={() => shown.set(!shown.value)}>
            Show or hide
          </Button>
        </div>
      );
    }),
};
