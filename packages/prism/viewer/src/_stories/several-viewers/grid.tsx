/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from '@molstar/prism-ui/dom/jsx-types';
import { onCleanup } from '@molstar/prism-ui/dom/scope';
import type { UITheme } from '@molstar/prism-ui/shell/theme';
import { DefaultPrismViewerSpec } from '../../default-spec.js';
import { createPrismViewer } from '../../page.js';
import { memoryStore, webStore } from '../../persist/stores.js';
import type { PrismViewer } from '../../viewer.js';

// Several real viewers on one page (README decision 29): a grid of hosts, each its own `createPrismViewer` with a persist
// key of its own. The stories' shared part.

export interface ViewerCell {
  /** The viewer's persist key: its stored state (theme, layout, kept records) is kept under it; one per live viewer. */
  readonly key: string;
  readonly title: string;
  readonly theme?: UITheme;
  /** After the viewer is made: load what it shows. */
  readonly setup: (viewer: PrismViewer) => void | Promise<void>;
}

export interface ViewerGridProps {
  readonly cells: readonly ViewerCell[];
  /** The grid's smallest column, in px: the columns that fit. */
  readonly minWidth?: number;
  /** The preferences go to localStorage under each key (a viewer keeps its theme and layout); else they are not kept. */
  readonly keep?: boolean;
  /** Above the grid (a story's controls). */
  readonly header?: Child;
  /** Once every viewer is made and set up; what it returns runs when the story goes. */
  readonly ready?: (viewers: readonly PrismViewer[]) => void | (() => void);
}

// A story that is shown again makes its viewers while the last one's are going: a persist key is free only once the
// viewer that had it is disposed, so the next waits for the disposals.
let disposing: Promise<unknown> = Promise.resolve();

const cellStyle = 'display: flex; flex-direction: column; min-width: 0; min-height: 0; border: 1px solid #8884';
const titleStyle =
  'font: 12px/1.6 system-ui, sans-serif; padding: 2px 8px; display: flex; justify-content: space-between; gap: 8px';

/** A grid of viewers, each in a host of its own. */
export function ViewerGrid(props: ViewerGridProps) {
  const root = document.createElement('div');
  root.style.cssText =
    'height: 100vh; display: flex; flex-direction: column; gap: 4px; box-sizing: border-box; padding: 4px';
  const grid = document.createElement('div');
  grid.style.cssText = `flex: 1; min-height: 0; display: grid; gap: 4px; grid-template-columns: repeat(auto-fit, minmax(${props.minWidth ?? 420}px, 1fr)); grid-auto-rows: minmax(320px, 1fr); overflow: auto`;
  if (props.header !== undefined) {
    const header = document.createElement('div');
    header.append(props.header as Node);
    root.append(header);
  }
  root.append(grid);

  const hosts = props.cells.map((cell) => {
    const frame = document.createElement('div');
    frame.style.cssText = cellStyle;
    frame.dataset.viewer = cell.key;
    const title = document.createElement('div');
    title.style.cssText = titleStyle;
    const name = document.createElement('strong');
    name.textContent = cell.title;
    const key = document.createElement('code');
    key.textContent = `persist.key: ${cell.key}`;
    title.append(name, key);
    const host = document.createElement('div');
    host.style.cssText = 'flex: 1; min-height: 0; position: relative';
    frame.append(title, host);
    grid.append(frame);
    return host;
  });

  const viewers: PrismViewer[] = [];
  let gone = false;
  let stop: (() => void) | void;
  onCleanup(() => {
    gone = true;
    stop?.();
    const own = viewers.splice(0);
    disposing = Promise.all([disposing, ...own.map((v) => v.dispose())]);
  });

  // Mounted first (a viewer measures its host), then made, one after the other.
  queueMicrotask(async () => {
    await disposing;
    for (const [i, cell] of props.cells.entries()) {
      if (gone) return;
      const base = DefaultPrismViewerSpec();
      const viewer = await createPrismViewer({
        target: hosts[i],
        spec: {
          ...base,
          ...(cell.theme ? { ui: { ...base.ui, theme: cell.theme } } : {}),
          persist: {
            key: cell.key,
            preferences: props.keep
              ? webStore('local', { prefix: cell.key, version: 'stories', budget: 64 * 1024 })
              : memoryStore({ version: 'stories', budget: 64 * 1024 }),
          },
        },
      });
      if (gone) return void viewer.dispose();
      viewers.push(viewer);
    }
    await Promise.all(
      props.cells.map(async (cell, i) => {
        try {
          await cell.setup(viewers[i]);
        } catch (e) {
          console.error(`[story] ${cell.key}:`, e);
        }
      }),
    );
    if (!gone) stop = props.ready?.(viewers);
  });
  return root;
}
