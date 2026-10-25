/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ActionFooter, domView, effect, ParamsForm, SectionList, Show } from '@molstar/prism-ui';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { createAnimationModel } from './model.js';
import type { WiggleState } from './source.js';

import { derive } from '@molstar/prism-ui/reactive/state';
// The Animation card (parity.md §4; the mock's `registerMockAnimation`): a list of sections in the `animation` slot,
// like Export's. Play chooses an animation of the plugin's (with its own params) and plays or stops it; Wiggle is
// Mol*'s procedural animation. Both act through the person's API; other extensions add sections to the same slot.

/** What moves now, under the Wiggle section's Apply. */
function wiggleInfo(w: WiggleState): string {
  if (w.mode === 'off') return w.tumble > 0 ? 'Tumbling' : '';
  return `Wiggling by ${w.mode}${w.tumble > 0 ? ', tumbling' : ''}`;
}

export function animationUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api.animation;
  if (!area) return {};
  // One model for the card's views: the sections share what plays.
  const model = createAnimationModel(ctx.api, area);
  return {
    overlays: [
      {
        id: 'animation.overlay.animation',
        title: 'Animation',
        render: domView(({ ctx: view }) => <SectionList shell={view.shell} slot="animation" focus={view.focus} />),
      },
    ],
    sections: [
      {
        id: 'animation.section.play',
        slot: 'animation',
        title: 'Play',
        icon: 'play',
        order: 0,
        render: domView(() => {
          // Built again if the animations or their params changed since it last showed (a structure to unwind).
          model.refresh();
          const { form } = model;
          if (!form) return <p style="margin: 0">No animations are registered.</p>;
          return (
            <div style="display: grid; gap: var(--msp-space-3)">
              <ParamsForm form={form} />
              {/* Play, or Stop while playing: the footer's action follows. */}
              <Show
                when={model.isPlaying}
                fallback={() => (
                  <ActionFooter
                    icon="play"
                    label="Play"
                    info={derive(model.blocked, (b) => b ?? '')}
                    disabled={derive(model.blocked, (b) => b !== undefined)}
                    onReset={() => form.reset()}
                    onApply={() => model.play()}
                  />
                )}
              >
                {() => <ActionFooter icon="square" label="Stop" onApply={() => model.stop()} />}
              </Show>
            </div>
          );
        }),
      },
      {
        id: 'animation.section.wiggle',
        slot: 'animation',
        title: 'Wiggle',
        icon: 'orbit',
        keywords: ['wiggle', 'dynamics', 'uncertainty', 'b-factor', 'tumble', 'motion'],
        order: 1,
        // Mol*'s wiggle: motion from noise (dynamics) or from model uncertainty, and whole copies tumbling (the scene's
        // animation options: `look`'s defaults); Clear stops both.
        render: domView(() => {
          // The scene's wiggle, followed (set elsewhere: a command, a session, Undo); a paused section reads it again
          // when it shows (`effect` resumes with its scope).
          effect(model.wiggle, () => model.refreshWiggle());
          return (
            <div style="display: grid; gap: var(--msp-space-3)">
              <ParamsForm form={model.wiggleForm} />
              <ActionFooter
                icon="orbit"
                label="Apply"
                info={derive(model.wiggle, wiggleInfo)}
                secondary={{ icon: 'x', label: 'Clear', onClick: () => model.clearWiggle() }}
                onApply={() => model.applyWiggle()}
              />
            </div>
          );
        }),
      },
    ],
    toolbar: [
      {
        id: 'animation.toolbar.animation',
        icon: 'clapperboard',
        title: 'Animation',
        group: 'scene',
        overlay: 'animation.overlay.animation',
      },
    ],
  };
}
