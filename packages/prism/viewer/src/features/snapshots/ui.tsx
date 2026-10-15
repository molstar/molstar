/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, onCleanup, Show, SnapshotCard } from '@molstar/prism-ui';
import { map } from 'rxjs';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { SnapshotsPanel } from '../../views/state/snapshots-panel.js';
import { createSnapshotStepper, createSnapshotsModel } from './model.js';

// The Snapshots tab (parity.md §4: the classic State snapshots) and the playback card. The shell has no place for a
// card at the canvas's top left that an extension can fill, so the playback card is an overlay opened from a toolbar
// button; the views read the snapshot area's state and act through the person's API.

/** Saves a blob as a file through the page. */
function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  a.style.display = 'none';
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function snapshotsUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api.snapshot;
  if (!area) return {};
  return {
    tabs: [
      {
        id: 'snapshots.tab.snapshots',
        title: 'Snapshots',
        icon: 'images',
        order: 1.5,
        render: domView(() => {
          const model = createSnapshotsModel(ctx.api, area, downloadBlob);
          onCleanup(() => model.dispose());
          return <SnapshotsPanel model={model} />;
        }),
      },
    ],
    overlays: [
      {
        id: 'snapshots.overlay.playback',
        title: 'Snapshots',
        render: domView(() => {
          const model = createSnapshotStepper(ctx.api, area);
          return (
            <Show
              when={model.count.pipe(map((c) => c >= 2))}
              fallback={() => <p style="margin: 0">Add two or more snapshots to play them.</p>}
            >
              {() => <SnapshotCard model={model} />}
            </Show>
          );
        }),
      },
    ],
    toolbar: [
      {
        id: 'snapshots.toolbar.playback',
        icon: 'play',
        title: 'Play snapshots',
        group: 'scene',
        overlay: 'snapshots.overlay.playback',
        visible: area.state.pipe(map((s) => s.entries.length >= 2)),
      },
    ],
  };
}
