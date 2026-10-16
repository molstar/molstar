/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandDefinition } from '@molstar/prism-ui/core';
import { ApiResult } from '../../api/area.js';
import type { ViewerCommandContext } from '../../api/managers.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { AnimationApi } from './operations.js';

declare module '../../api/area.js' {
  // The `animation` area (operations.md §2): the plugin's animations and the wiggle as operations.
  interface BuiltInApis {
    animation: AnimationApi;
  }
}

// The `animation` feature (viewer-spec.md §3.1): the `animation` area and the Animation card (play, procedural
// animation). The views are a chunk of their own. Stepping through the models of a trajectory is the core's
// `trajectory` area. One command plays an animation, its form listing those registered when it is shown (an animation
// registered later is in it); there is no command per animation, which would be made once, as the feature starts.

export const AnimationFeature: ViewerExtension = {
  id: 'animation',
  title: 'Animation',
  builtIn: true,
  api: { animation: (ctx) => new AnimationApi(ctx) },
  commands: (api) => {
    const a = api.animation;
    if (!a) return [];
    const commands: CommandDefinition<ViewerCommandContext>[] = [
      operationCommand(api, a.play, {
        category: 'Animation',
        icon: 'play',
        title: 'Play an animation…',
        alias: ['animate', 'animation', 'spin', 'rock', 'turntable', 'rotate', 'unwind', 'explode', 'play animation'],
        action: 'Play',
        // The form takes the animation, among those registered when it is shown (one registered later is in it); its
        // params are the card's (`animation.section.play`).
        values: { params: {} },
        form: (params) => {
          const animations = a.state.value.animations;
          if (!animations.length) return params;
          const first = animations.find((x) => x.available === true) ?? animations[0];
          return {
            ...params,
            name: PD.Select(
              first.name,
              animations.map(
                (x) => [x.name, x.available === true ? x.title : `${x.title} (not now)`] as [string, string],
              ),
              { label: 'Animation' },
            ),
          };
        },
      }),
      operationCommand(api, a.stop, { category: 'Animation', icon: 'square', alias: ['pause', 'stop animation'] }),
      operationCommand(api, a.wiggle, {
        category: 'Animation',
        icon: 'orbit',
        title: 'Wiggle…',
        alias: ['wiggle', 'dynamics', 'uncertainty', 'b-factor', 'procedural'],
        action: 'Apply',
      }),
      operationCommand(api, a.wiggle, {
        id: 'animation.wiggle.off',
        title: 'Stop the wiggle',
        category: 'Animation',
        icon: 'x',
        alias: ['clear wiggle', 'no wiggle'],
        values: { mode: 'off' },
      }),
      {
        id: 'animation.open',
        title: 'Show the Animation card',
        category: 'Animation',
        icon: 'clapperboard',
        alias: ['animation', 'animate', 'playback'],
        // Opens the card: a view, not a change.
        async run() {
          const r = await api.view.open({ card: 'animation.overlay.animation' });
          if (ApiResult.isFailure(r)) throw new Error(r.message);
          return { message: 'Opened the Animation card.', ops: 'view.open({ card: "animation.overlay.animation" })' };
        },
      },
    ];
    return commands;
  },
  ui: async (ctx) => (await import('./ui.js')).animationUI(ctx),
};
