/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { createElement } from 'react';
import { createRoot } from 'react-dom/client';
import { onCleanup } from '../../dom/scope.js';
import { observableState } from '../../reactive/observable-state.js';
import { mockBlankMode, ShellStory } from '../../testing/shell.js';
import { story } from '../../testing/story-root.js';
import { preactPanelView } from '../_stories/preact/preact-parts.js';
import { ReactApp, ReactPanel } from '../_stories/react-parts.js';
import { reactView } from '../react-view.js';

const meta: Meta = { title: 'Internals/React and Preact', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * Contributions written in React and in Preact, in the shell next to the built-in tabs (the rail's React and Preact
 * buttons): both mount through `reactView` and read and write the same prism-ui state with `useObservable`, and post to
 * the status line. Add one in either; the other, and the status line, follow.
 */
export const ReactAndPreactTabs: StoryObj = {
  render: () =>
    story(() => {
      const count = observableState(0);
      return (
        <div style="height: 100vh">
          <ShellStory
            loaded={true}
            mode={mockBlankMode()}
            addOn={({ status, shell }) => {
              const post = (text: string) => status.post(text);
              queueMicrotask(() => shell.show('react'));
              return {
                tabs: [
                  {
                    id: 'react',
                    title: 'React',
                    icon: 'orbit',
                    order: 2.8,
                    render: reactView(() => createElement(ReactPanel, { count, post }), { createRoot }),
                  },
                  {
                    id: 'preact',
                    title: 'Preact',
                    icon: 'hexagon',
                    order: 2.9,
                    render: preactPanelView({ count, post }),
                  },
                ],
              };
            }}
          />
        </div>
      );
    }),
};

/**
 * prism-ui inside a React app: the header is React (with the same state as a prism-ui view would read), the shell is a
 * prism-ui element that React hosts with `MspElementHost`, its model set as a property.
 */
export const InAReactApp: StoryObj = {
  render: () =>
    story(() => {
      const count = observableState(0);
      const host = document.createElement('div');
      const root = createRoot(host);
      root.render(createElement(ReactApp, { count, shell: () => <ShellStory loaded={true} mode={mockBlankMode()} /> }));
      onCleanup(() => root.unmount());
      return host;
    }),
};
