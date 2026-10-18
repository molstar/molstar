/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, For, Icon, Markdown, ScrollArea, Show } from '@molstar/prism-ui';
import type { StoryModel } from '../model.js';
import { classes } from './_generated/story-tab.module.css.js';

import { derive } from '@molstar/prism-ui/reactive/state';
/**
 * Story mode's main tab (a first sketch): the story's scenes in order, the current one picked, a click shows it; the
 * current scene's caption under the list; Add scene takes what the canvas shows.
 */
export function StoryTab(props: { story: StoryModel }) {
  const { story } = props;
  const current = derive([story.scenes, story.current], (scenes, id) => scenes.find((s) => s.id === id));
  return (
    <div class={classes.tab}>
      <header class={classes.header}>
        <span class={classes.title}>{story.title}</span>
        <Button size="sm" variant="ghost" icon={<Icon name="plus" size={14} />} onClick={() => story.addScene()}>
          Add scene
        </Button>
      </header>
      <ScrollArea edge class={classes.scenes} role="listbox" aria-label="Scenes">
        <For each={story.scenes} key={(s) => s.id}>
          {(scene, index) => (
            <div
              role="option"
              class={classes.scene}
              aria-selected={derive(story.current, (id) => id === scene.id)}
              onClick={() => story.select(scene.id)}
            >
              <span class={classes.number}>{derive(index, (i) => i + 1)}</span>
              <span class={classes.sceneTitle} data-tip-overflow="">
                {scene.title}
              </span>
            </div>
          )}
        </For>
      </ScrollArea>
      <Show when={current}>
        {(scene) => (
          <section class={classes.caption} aria-label="Caption">
            <ScrollArea class={classes.captionScroll} contentClass={classes.captionContent}>
              {scene.caption ? <Markdown source={scene.caption} /> : <p class={classes.empty}>No caption yet.</p>}
            </ScrollArea>
          </section>
        )}
      </Show>
    </div>
  );
}
