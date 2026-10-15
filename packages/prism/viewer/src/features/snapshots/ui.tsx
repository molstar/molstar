/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  Button,
  domView,
  Icon,
  onActive,
  onCleanup,
  Show,
  SnapshotCard,
  SnapshotRow,
  TransitionRow,
} from '@molstar/prism-ui';
import { distinctUntilChanged, map } from 'rxjs';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { downloadBlob } from '../../ui/download.js';
import { SnapshotsPanel } from '../../views/state/snapshots-panel.js';
import { AudioBar } from '../../views/state/audio-bar.js';
import { createAudioModel } from './audio.js';
import { createMarkdownHost } from './markdown.js';
import {
  createSnapshotStepper,
  createSnapshotTransition,
  createSnapshotsModel,
  openPlaybackOnLoad,
  PlaybackCard,
  presentable,
} from './model.js';

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
          // Play or stop the slideshow (two or more snapshots): an icon at the header's end, beside the options.
          return (
            <SnapshotsPanel
              model={model}
              actions={() => (
                <Show when={stepper.count.pipe(map((c) => c >= 2))}>
                  {() => (
                    <Button
                      size="sm"
                      variant="ghost"
                      icon={
                        <Show when={stepper.isPlaying} fallback={() => <Icon name="play" />}>
                          {() => <Icon name="pause" />}
                        </Show>
                      }
                      title={stepper.isPlaying.pipe(map((p) => (p ? 'Stop the slideshow' : 'Play the snapshots')))}
                      onClick={() => stepper.togglePlay()}
                    />
                  )}
                </Show>
              )}
            />
          );
        }),
      },
    ],
    // The snapshots' rows in the canvas's playback card: step through them, and the current one's transition.
    sections: [
      {
        id: 'snapshots.section.playback',
        slot: 'playback',
        title: 'Snapshots',
        render: domView(() => (
          <>
            <SnapshotRow model={createSnapshotStepper(ctx.api, area)} />
            <TransitionRow model={createSnapshotTransition(ctx.api, area)} />
          </>
        )),
      },
    ],
    overlays: [
      {
        id: PlaybackCard,
        title: 'Story',
        // Read beside the scene: a click on the canvas or a panel leaves it open; the description fills the card.
        persistent: true,
        fill: true,
        render: domView(() => {
          const model = createSnapshotStepper(ctx.api, area, createMarkdownHost(ctx.api, area.markdown));
          return (
            <div style="display: flex; flex: 1 1 0; flex-direction: column; gap: var(--msp-space-2); min-height: 0">
              <Show when={audio.present}>{() => <AudioBar model={audio} />}</Show>
              <Show
                when={model.count.pipe(map((c) => c >= 1))}
                fallback={() => <p style="margin: 0">There are no snapshots: add some in the Snapshots tab.</p>}
              >
                {() => <SnapshotCard model={model} bare />}
              </Show>
            </div>
          );
        }),
      },
    ],
    toolbar: [
      {
        id: 'snapshots.toolbar.playback',
        icon: 'book-open',
        title: 'Story',
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
