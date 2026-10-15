/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Icon, onActive, onCleanup, Show, SnapshotCard } from '@molstar/prism-ui';
import { map } from 'rxjs';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { SnapshotsPanel } from '../../views/state/snapshots-panel.js';
import { createSnapshotStepper, createSnapshotsModel } from './model.js';

// The Snapshots tab (parity.md §4: the classic State snapshots) and the playback card. The shell has no place for a
// card at the canvas's top left that an extension can fill, so the playback card is an overlay opened from a toolbar
// button; the views read the snapshot area's state and act through the person's API.

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
          const stepper = createSnapshotStepper(ctx.api, area);
          onActive((active) => model.setActive(active));
          onCleanup(() => model.dispose());
          // Play or stop the slideshow from the tab, when there is one (two or more snapshots); the panel fills the rest.
          return (
            <div style="display: flex; flex-direction: column; flex: 1 1 0; min-height: 0; height: 100%">
              <Show when={stepper.count.pipe(map((c) => c >= 2))}>
                {() => (
                  <div style="flex: none; display: flex; padding: 0 var(--msp-panel-inset) var(--msp-space-1)">
                    <Button
                      size="sm"
                      variant="ghost"
                      icon={
                        <Show when={stepper.isPlaying} fallback={() => <Icon name="play" />}>
                          {() => <Icon name="pause" />}
                        </Show>
                      }
                      onClick={() => stepper.togglePlay()}
                    >
                      <span>{stepper.isPlaying.pipe(map((p) => (p ? 'Stop' : 'Play')))}</span>
                    </Button>
                  </div>
                )}
              </Show>
              <SnapshotsPanel model={model} />
            </div>
          );
        }),
      },
    ],
    overlays: [
      {
        id: 'snapshots.overlay.playback',
        title: 'Playback',
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
