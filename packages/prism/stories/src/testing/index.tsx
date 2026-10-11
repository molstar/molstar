/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, observableState, SnapshotCard } from '@molstar/prism-ui';
import { CommandManager } from '@molstar/prism-ui/commands/manager';
import type { CommandDefinition } from '@molstar/prism-ui/commands/types';
import type { SnapshotStepperModel } from '@molstar/prism-ui/model/playback';
import { CanvasPlaceholder, type MockMode } from '@molstar/prism-ui/testing/index';
import { of } from 'rxjs';
import type { StoryModel, StoryScene } from '../model.js';
import { StoryTab } from '../views/story-tab.js';

const Tour: StoryScene[] = [
  {
    id: 'overview',
    title: 'Cytochrome P450 3A4',
    caption: 'The enzyme that clears **half of all drugs** from the body.',
  },
  { id: 'heme', title: 'The heme', caption: 'Iron at the centre of a porphyrin ring, held by `Cys442`.' },
  {
    id: 'pocket',
    title: 'A large, flexible pocket',
    caption: 'Room for many different drugs, and sometimes two at once.',
  },
  { id: 'channel', title: 'Substrate channel', caption: 'How a drug gets in from the membrane side.' },
];

/** A story in memory (the CYP3A4 tour): scenes to select, add and remove, and a player over them. */
export function mockStory(post: (text: string) => void): StoryModel {
  const scenes = observableState<readonly StoryScene[]>(Tour);
  const current = observableState<string | undefined>(Tour[0].id);
  const isPlaying = observableState(false);
  let timer: ReturnType<typeof setInterval> | undefined;
  let added = 0;
  const index = () => scenes.value.findIndex((s) => s.id === current.value);
  const go = (i: number) => {
    const list = scenes.value;
    if (!list.length) return;
    current.set(list[(i + list.length) % list.length].id);
  };
  // The player's states, kept in step with the scenes and the current one.
  const count = observableState(Tour.length);
  const at = observableState(0);
  const title = observableState<string | undefined>(Tour[0].title);
  const titles = observableState<readonly string[]>(Tour.map((t) => t.title));
  const description = observableState<string | undefined>(Tour[0].caption);
  const sync = () => {
    const i = index();
    const scene = scenes.value[i];
    count.set(scenes.value.length);
    at.set(i);
    title.set(scene?.title);
    titles.set(scenes.value.map((x) => x.title));
    description.set(scene?.caption);
  };
  scenes.asReadonly.subscribe(sync);
  current.asReadonly.subscribe(sync);
  const player: SnapshotStepperModel = {
    count: count.asReadonly,
    index: at.asReadonly,
    title: title.asReadonly,
    titles: titles.asReadonly,
    description: description.asReadonly,
    isPlaying: isPlaying.asReadonly,
    previous: () => go(index() - 1),
    next: () => go(index() + 1),
    goTo: go,
    togglePlay: () => {
      isPlaying.set(!isPlaying.value);
      clearInterval(timer);
      if (isPlaying.value) timer = setInterval(() => go(index() + 1), 2500);
    },
    openLink: (href) => post(`Link: ${href}`),
  };
  return {
    title: observableState('How the liver clears drugs').asReadonly,
    scenes: scenes.asReadonly,
    current: current.asReadonly,
    select: (id) => current.set(id),
    addScene: () => {
      const scene = { id: `scene-${++added}`, title: `New scene ${added}` };
      const list = scenes.value.slice();
      list.splice(index() + 1, 0, scene);
      scenes.set(list);
      current.set(scene.id);
      post(`Added “${scene.title}” from the current view`);
    },
    removeScene: (id) => scenes.set(scenes.value.filter((s) => s.id !== id)),
    player,
  };
}

/**
 * Story mode in the mock shell (a first sketch, its UX still to be designed): the Story tab (the scenes, the current
 * caption), the scenes playing on the canvas (the playback card), and story commands instead of Mol*'s: no state tree,
 * no structure tools; what a story does edits the story.
 */
export function mockStoryMode(): MockMode {
  return {
    id: 'story',
    title: 'Story mode',
    icon: 'play',
    activate: ({ shell, status, background }) => {
      const post = (text: string) => status.post(text);
      const story = mockStory(post);
      shell.registerTabs([
        {
          id: 'story',
          compactBar: true,
          title: 'Story',
          icon: 'play',
          order: 0,
          place: 'head',
          render: domView(() => <StoryTab story={story} />),
        },
      ]);
      const done = (message: string) => ({ message });
      const commands = new CommandManager<object>({
        context: observableState({}).asReadonly,
        formats: { get: () => undefined },
      });
      const definitions: CommandDefinition<object>[] = [
        {
          id: 'story.add-scene',
          title: 'Add scene',
          alias: ['new scene', 'capture'],
          category: 'Story',
          icon: 'plus',
          run: () => {
            story.addScene();
            return done('');
          },
        },
        {
          id: 'story.play',
          title: 'Play story',
          alias: ['present', 'run'],
          category: 'Story',
          icon: 'play',
          run: () => {
            story.player.togglePlay();
            return done('Playing the story');
          },
        },
        {
          id: 'story.next',
          title: 'Next scene',
          category: 'Story',
          icon: 'step-forward',
          shortcut: { key: 'ArrowRight', alt: true },
          run: () => {
            story.player.next();
            return done('');
          },
        },
        {
          id: 'story.previous',
          title: 'Previous scene',
          category: 'Story',
          icon: 'step-back',
          shortcut: { key: 'ArrowLeft', alt: true },
          run: () => {
            story.player.previous();
            return done('');
          },
        },
        {
          id: 'story.export',
          title: 'Export story',
          alias: ['mvsx', 'download', 'publish'],
          category: 'Story',
          icon: 'download',
          run: () => done('Saved story.mvsx'),
        },
      ];
      commands.register(definitions);
      return {
        start: 'story',
        commands,
        canvas: () => <CanvasPlaceholder empty={false} background={background} />,
        playback: () => <SnapshotCard model={story.player} />,
        empty: of(false),
      };
    },
  };
}
