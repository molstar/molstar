/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { CommandDefinition } from '@molstar/ui/commands/types';
import { mountKeysPanel } from '@molstar/seqx-viewer/mount-keys-panel';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The keys panel (_spec/README.md "Behaviors", keysPanel): `help.keys` (`?`, and the toolbar's "Keys and gestures (?)")
// opens a panel over the canvas listing what the view answers to (`describeKeys`, generated); `help.keys` again or Esc
// closes it. Whether it is open is the reader's session, never the spec: the command changes nothing in the models.

export const KeysPanelCommandId = 'help.keys';

export interface KeysPanelBehavior {
  /** Whether the panel is open. */
  readonly open: ReadonlyState<boolean>;
  readonly commands: readonly CommandDefinition<unknown>[];
  /** Closes the panel; true when it was open (Esc's). */
  close(): boolean;
  /** The mounted canvas's element and its viewer: the panel is drawn over it while open. */
  attach(container: HTMLElement, viewer: SeqViewer): () => void;
}

export function keysPanelBehavior(): KeysPanelBehavior {
  const open = new State(false);
  return {
    open: open.asReadonly,
    commands: [
      {
        id: KeysPanelCommandId,
        title: 'Keys and gestures',
        category: 'View',
        icon: 'keyboard',
        scope: 'canvas',
        run: () => {
          open.set(!open.value);
          return { message: open.value ? 'Keys and gestures shown' : 'Keys and gestures hidden' };
        },
      },
    ],
    close() {
      if (!open.value) return false;
      open.set(false);
      return true;
    },
    attach(container, viewer) {
      return mountKeysPanel(container, viewer, open.asReadonly, () => open.set(false));
    },
  };
}

/** The keys panel (`keysPanel`, in the default spec). */
export const KeysPanelProvider: BehaviorProvider<undefined, KeysPanelBehavior> = {
  name: 'keysPanel',
  create: () => keysPanelBehavior(),
};
