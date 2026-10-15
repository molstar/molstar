/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Icon, onActive, onCleanup, Show, SnapshotCard } from '@molstar/prism-ui';
import { distinctUntilChanged, map } from 'rxjs';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { SnapshotsPanel } from '../../views/state/snapshots-panel.js';
import { AudioBar } from '../../views/state/audio-bar.js';
import { createAudioModel } from './audio.js';
import { createMarkdownHost } from './markdown.js';
import { createSnapshotStepper, createSnapshotsModel, openPlaybackOnLoad, PlaybackCard, presentable } from './model.js';

// The Snapshots tab (parity.md §4: the classic State snapshots) and the playback card. The shell has no place for a
// card at the canvas's top left that an extension can fill, so the playback card is an overlay opened from a toolbar
// button, and by itself when a story is loaded (`openPlaybackOnLoad`); the views read the snapshot area's state and act through
// the person's API.

export function snapshotsUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api.snapshot;
  if (!area) return {};
  // As long as the viewer: the run events complete with it (the views have no dispose of their own).
  openPlaybackOnLoad(ctx.api, area);
  // A story's audio: its bar in the playback card, and a canvas button while there is a player (the card may be shut).
  const audio = createAudioModel(area.markdown);
  return {
    tabs: [
      {
        id: 'snapshots.tab.snapshots',
        title: 'Snapshots',
        icon: 'snapshots',
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
        id: PlaybackCard,
        title: 'Playback',
        render: domView(() => {
          const model = createSnapshotStepper(ctx.api, area, createMarkdownHost(ctx.api, area.markdown));
          return (
            <div style="display: grid; gap: var(--msp-space-2)">
              <Show when={audio.present}>{() => <AudioBar model={audio} />}</Show>
              <Show
                when={model.count.pipe(map((c) => c >= 1))}
                fallback={() => <p style="margin: 0">There are no snapshots: add some in the Snapshots tab.</p>}
              >
                {() => <SnapshotCard model={model} />}
              </Show>
            </div>
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
        overlay: PlaybackCard,
        visible: area.state.pipe(map(presentable), distinctUntilChanged()),
      },
      {
        id: 'snapshots.toolbar.audio',
        icon: 'play',
        title: 'Play audio',
        group: 'view',
        place: 'canvas',
        visible: audio.present,
        face: audio.playing.pipe(
          map((p) =>
            p ? { icon: 'pause' as const, title: 'Pause audio' } : { icon: 'play' as const, title: 'Play audio' },
          ),
        ),
        run: () => audio.toggle(),
        menu: [
          { id: 'stop', label: 'Stop', icon: 'square', run: () => audio.stop() },
          { id: 'close', label: 'Close the audio', icon: 'x', run: () => audio.close() },
        ],
      },
    ],
  };
}
