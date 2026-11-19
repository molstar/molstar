/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/core/reactive';
import { describeKeys } from './keys.js';
import type { SeqViewer } from './viewer.js';

// The keys panel's part (_spec/README.md "Behaviors", keysPanel): `@molstar/ui`'s primitives over the canvas's corner: a
// surface with a header and a close button, a header a group and a row an entry (`msui-row`: its name, then its keys and
// gestures in the trailing column, muted). Drawn from `describeKeys` each time it opens; its events are its own (a click
// in it is no click on the canvas). Everything it adds is taken away by the function it returns.

export function mountKeysPanel(
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
  // the panel's events are its own: the canvas under it takes none of them
  const own = (e: Event) => e.stopPropagation();
  const types = ['pointerdown', 'pointermove', 'pointerup', 'click', 'dblclick', 'wheel'];
  for (const t of types) panel.addEventListener(t, own);
  container.append(panel);

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
    x.addEventListener('click', () => {
      close();
      container.focus({ preventScroll: true });
    });
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
    if (shown) draw();
    panel.style.display = shown ? '' : 'none';
  };
  sync(open.value);
  const stop = open.listen(sync);
  return () => {
    stop();
    for (const t of types) panel.removeEventListener(t, own);
    panel.remove();
  };
}
