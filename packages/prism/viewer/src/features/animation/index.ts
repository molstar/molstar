/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition, IconName } from '@molstar/prism-ui/core';
import { ApiResult } from '../../api/area.js';
import type { ViewerCommandContext } from '../../api/managers.js';
import { operationCommand } from '../../commands/operations.js';
import type { ViewerExtension } from '../../spec.js';
import { AnimationApi } from './operations.js';

declare module '../../api.js' {
  // The `animation` area (operations.md §2): the plugin's animations and the wiggle as operations.
  interface ExtensionApis {
    animation: AnimationApi;
  }
}

// The `animation` feature (viewer-spec.md §3.1): the `animation` area and the Animation card (play, procedural
// animation). The views are a chunk of their own. Stepping through the models of a trajectory is the core's
// `trajectory` area.

/** Animations with a command of their own, by registered name: a preset of `animation.play` with its title. */
const Presets: readonly { name: string; title: string; icon: IconName; alias: readonly string[] }[] = [
  {
    name: 'built-in.animate-camera-spin',
    title: 'Spin the camera',
    icon: 'rotate-3d',
    alias: ['spin', 'turntable', 'rotate'],
  },
  { name: 'built-in.animate-camera-rock', title: 'Rock the camera', icon: 'rotate-3d', alias: ['rock', 'sway'] },
  {
    name: 'built-in.animate-structure-spin',
    title: 'Spin the structure',
    icon: 'rotate-3d',
    alias: ['spin structure'],
  },
  {
    name: 'built-in.animate-assembly-unwind',
    title: 'Unwind the assembly',
    icon: 'play',
    alias: ['unwind', 'assembly'],
  },
  { name: 'built-in.animate-units-explode', title: 'Explode the units', icon: 'play', alias: ['explode'] },
];

export const AnimationFeature: ViewerExtension = {
  id: 'animation',
  title: 'Animation',
  builtIn: true,
  api: { animation: (ctx) => new AnimationApi(ctx) },
  commands: (api) => {
    const a = api.animation;
    if (!a) return [];
    const registered = new Set(a.state.value.animations.map((x) => x.name));
    const commands: CommandDefinition<ViewerCommandContext>[] = [
      operationCommand(api, a.play, {
        category: 'Animation',
        icon: 'play',
        title: 'Play an animation…',
        alias: ['animate', 'animation'],
        action: 'Play',
        // The form takes the animation; its params are the card's (`animation.section.play`).
        values: { params: {} },
      }),
      ...Presets.filter((p) => registered.has(p.name)).map((p) =>
        operationCommand(api, a.play, {
          id: `animation.play.${p.name.replace('built-in.animate-', '')}`,
          title: p.title,
          category: 'Animation',
          icon: p.icon,
          alias: p.alias,
          values: { name: p.name, params: {} },
        }),
      ),
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
