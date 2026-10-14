/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the Seq* stories need of the host's preview, picked up by Mol*'s Storybook (storybook/packages.ts, `preview`):
// Prism's tokens and the page defaults (the toolbar's theme sets `data-msp-theme` on <html>, as for Mol*'s own
// stories), and the decorator that mounts a Seq* story. The stories are React and the host is Storybook for web
// components (`render` returns a React element, which that renderer cannot show), so each is rendered into a node of
// its own by `react-dom/client` and unmounted with the story.

import { createElement, type ReactElement } from 'react';
import { createRoot, type Root } from 'react-dom/client';
import { addons } from 'storybook/preview-api';
import { STORY_CHANGED } from 'storybook/internal/core-events';
import './_stories/ui/tokens.css';

/** The sidebar group of every Seq* story and docs page (their titles all start with it). */
export const SEQ_STORIES_GROUP = 'Seq*/';

interface SeqStoryContext {
  id: string;
  title: string;
  args: unknown;
  originalStoryFn: (args: any, context: any) => unknown;
}

/**
 * The story's `render` as a component, so React calls it (the stories' renders use hooks, as `@storybook/react` allows):
 * it gets the args and context of the decorator's call.
 */
function SeqStory({ context }: { context: SeqStoryContext }) {
  return context.originalStoryFn(context.args, context) as ReactElement;
}

/** The story shown now: its node and React root. */
let live: { id: string; host: HTMLElement; root: Root } | undefined;
let listening = false;

function unmountLive() {
  live?.root.unmount();
  live = undefined;
}

function mountReact(context: SeqStoryContext): HTMLElement {
  if (!listening) {
    listening = true;
    // Another story or a docs page is chosen (or the story is gone): what it showed is unmounted. (Not the decorator
    // hooks' effects: Storybook runs their cleanup only for a story whose render it saw complete.)
    addons.getChannel().on(STORY_CHANGED, unmountLive);
  }
  if (live?.id !== context.id) {
    unmountLive();
    const host = document.createElement('div');
    live = { id: context.id, host, root: createRoot(host) };
  }
  // An args change calls the decorator again: the same root renders the new element, keeping the components' state.
  live.root.render(createElement(SeqStory, { context }));
  return live.host;
}

/**
 * A decorator for Storybook's web components renderer: shows a Seq* story's React element, leaves other stories be. It
 * renders the story's own `render` and none of the decorators inside it (the Seq* stories have none).
 */
export function seqReactDecorator<R>(storyFn: () => R, context: SeqStoryContext): R | HTMLElement {
  return context.title.startsWith(SEQ_STORIES_GROUP) ? mountReact(context) : storyFn();
}

/** The decorators the host's preview applies (`StoryPackage.preview`). */
export const decorators = [seqReactDecorator];
