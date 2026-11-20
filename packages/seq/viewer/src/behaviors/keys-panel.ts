/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { CommandDefinition } from '@molstar/ui/commands/types';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The keys panel (_spec/README.md "Behaviors", keysPanel): `help.keys` (`?`, and the toolbar's "Keys and gestures (?)")
// opens a panel over the canvas listing what the view answers to (`describeKeys`, generated); `help.keys` again or Esc
// closes it. Whether it is open is the reader's session, never the spec: the command changes nothing in the models.
// The panel is non-modal: it takes the focus when it opens (its Close button) and gives it back to what opened it when
// it closes; the focus may leave it (no trap). Its keys are its own: none reaches the canvas under it. Its DOM is a
// lazy chunk (`keys-panel-dom.ts`), asked for at the first open.

export const KeysPanelCommandId = 'help.keys';

export interface KeysPanelBehavior {
  /** Whether the panel is open. */
  readonly open: ReadonlyState<boolean>;
  readonly commands: readonly CommandDefinition<unknown>[];
  /** Closes the panel; true when it was open (Esc's). */
  close(): boolean;
  /**
   * The mounted canvas's element and its viewer: the panel is drawn over it while open (its DOM a lazy chunk, asked
   * for at the first open).
   */
  attach(container: HTMLElement, viewer: SeqViewer): () => void;
  /** Resolves once the panel's DOM asked for landed (`settle` waits for it). */
  loaded(): Promise<void>;
}

function keysPanelBehavior(): KeysPanelBehavior {
  const open = new State(false);
  /** The panel's DOM asked for (the last mount's). */
  let landing: Promise<void> | undefined;
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
      let off: (() => void) | undefined;
      let asked = false;
      let gone = false;
      // the DOM at the first open: it draws the state it finds when it lands
      const ask = () => {
        if (asked) return;
        asked = true;
        landing = import('./keys-panel-dom.js').then((m) => {
          if (!gone) off = m.mountKeysPanel(container, viewer, open.asReadonly, () => open.set(false));
        });
      };
      if (open.value) ask();
      const stop = open.listen((shown) => shown && ask());
      return () => {
        gone = true;
        stop();
        off?.();
      };
    },
    loaded: () => landing ?? Promise.resolve(),
  };
}

/** The keys panel (`keysPanel`, in the default spec). */
export const KeysPanelProvider: BehaviorProvider<undefined, KeysPanelBehavior> = {
  name: 'keysPanel',
  create: () => keysPanelBehavior(),
};
