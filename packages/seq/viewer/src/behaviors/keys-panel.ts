/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ReadonlyState, State } from '@molstar/core/reactive';
import type { CommandDefinition } from '@molstar/ui/commands/types';
import { describeKeys } from '@molstar/seqx-viewer/keys';
import type { BehaviorProvider } from '@molstar/seqx-viewer/viewer-spec';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

// The keys panel (_spec/README.md "Behaviors", keysPanel): `help.keys` (`?`, and the toolbar's "Keys and gestures (?)")
// opens a panel over the canvas listing what the view answers to (`describeKeys`, generated); `help.keys` again or Esc
// closes it. Whether it is open is the reader's session, never the spec: the command changes nothing in the models.
// The panel is non-modal: it takes the focus when it opens (its Close button) and gives it back to what opened it when
// it closes; the focus may leave it (no trap). Its keys are its own: none reaches the canvas under it.

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

function keysPanelBehavior(): KeysPanelBehavior {
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

function mountKeysPanel(
  container: HTMLElement,
  viewer: SeqViewer,
  open: ReadonlyState<boolean>,
  close: () => void,
): () => void {
  const doc = container.ownerDocument;
  const panel = doc.createElement('div');
  panel.className = 'msui-surface seqx-keys-panel';
  panel.setAttribute('role', 'dialog');
  panel.setAttribute('aria-label', 'Keys and gestures');
  Object.assign(panel.style, {
    position: 'absolute',
    display: 'none',
    top: 'var(--msui-space-2)',
    right: 'var(--msui-space-2)',
    width: 'min(360px, calc(100% - 2 * var(--msui-space-2)))',
    // past a short canvas's bottom (it floats over what follows), never past the window's
    maxHeight: 'min(420px, 80vh)',
    minHeight: '120px',
    overflow: 'auto',
    zIndex: 'var(--msui-z-overlay)',
    boxShadow: '0 0 0 1px var(--msui-line), var(--msui-shadow-float)',
  });
  // the panel's events are its own: the canvas under it takes none of them, keys included (an arrow in the panel moves
  // no cursor, Enter on Close runs no command); Esc is the panel's own: it closes it, wherever the focus came from
  const own = (e: Event) => e.stopPropagation();
  const types = ['pointerdown', 'pointermove', 'pointerup', 'click', 'dblclick', 'wheel', 'keyup'];
  for (const t of types) panel.addEventListener(t, own);
  const onKey = (e: KeyboardEvent) => {
    e.stopPropagation();
    if (e.key === 'Escape') {
      e.preventDefault();
      close();
    } else if ((e.key === 'Enter' || e.key === ' ') && e.target instanceof HTMLButtonElement) {
      // the button's own press, at once (the platform's click after it is prevented: it would close twice)
      e.preventDefault();
      e.target.click();
    }
  };
  panel.addEventListener('keydown', onKey);
  container.append(panel);
  /** What had the focus when the panel opened: it gets it back when the panel closes (non-modal: no trap). */
  let opener: Element | null = null;

  const draw = () => {
    const header = doc.createElement('div');
    header.className = 'msui-header';
    header.style.justifyContent = 'space-between';
    const title = doc.createElement('span');
    title.textContent = 'Keys and gestures';
    const x = doc.createElement('button');
    x.type = 'button';
    x.className = 'msui-control';
    x.textContent = 'Close';
    x.title = 'Close (Esc)';
    x.addEventListener('click', () => close());
    header.append(title, x);
    const out: HTMLElement[] = [header];
    for (const g of describeKeys(viewer)) {
      const h = doc.createElement('div');
      h.className = 'msui-header';
      h.setAttribute('role', 'heading');
      h.setAttribute('aria-level', '3');
      h.textContent = g.title;
      const list = doc.createElement('ul');
      list.setAttribute('aria-label', g.title);
      Object.assign(list.style, { listStyle: 'none', margin: '0', padding: '0' });
      for (const e of g.entries) {
        const row = doc.createElement('li');
        row.className = 'msui-row msui-row-plain';
        // a long list of gestures wraps in its own column, never over the name
        Object.assign(row.style, {
          height: 'auto',
          minHeight: 'var(--msui-row)',
          columnGap: 'var(--msui-space-3)',
          gridTemplateColumns: 'minmax(0, 1fr) minmax(0, 1fr)',
          paddingBlock: '2px',
        });
        const name = doc.createElement('span');
        name.textContent = e.name;
        const how = doc.createElement('span');
        how.className = 'msui-label';
        how.style.textAlign = 'end';
        how.textContent = [...e.keys, ...e.gestures].join(', ');
        row.append(name, how);
        list.append(row);
      }
      out.push(h, list);
    }
    panel.replaceChildren(...out);
  };

  const sync = (shown: boolean) => {
    if (shown) {
      opener = doc.activeElement;
      draw();
      panel.style.display = '';
      panel.querySelector<HTMLButtonElement>('button')?.focus({ preventScroll: true });
      return;
    }
    const hadFocus = panel.contains(doc.activeElement);
    panel.style.display = 'none';
    // the focus goes back to what opened it (the canvas, a toolbar's button), else to the canvas
    if (hadFocus) {
      const back = opener instanceof HTMLElement && opener.isConnected ? opener : container;
      back.focus({ preventScroll: true });
    }
    opener = null;
  };
  if (open.value) sync(true);
  const stop = open.listen(sync);
  return () => {
    stop();
    for (const t of types) panel.removeEventListener(t, own);
    panel.removeEventListener('keydown', onKey);
    panel.remove();
  };
}
