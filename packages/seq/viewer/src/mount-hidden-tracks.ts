/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import {
  Toolbar,
  type ToolbarFace,
  type ToolbarItem,
  type ToolbarMenuRow,
  type ToolbarMenuSide,
} from '@molstar/ui/components/toolbar/toolbar';
import { installTooltips } from '@molstar/ui/components/tooltip/tooltip';
import { render } from '@molstar/ui/dom/render';
import { onCleanup } from '@molstar/ui/dom/scope';
import type { SeqViewer } from './viewer.js';

// The `HiddenTracks` part (_spec/README.md "Hidden tracks"): a part the host places outside the view (a toolbar, a
// panel), never drawn on the canvas, which draws nothing for a hidden track so a hide moves nothing. One button with a
// menu: each hidden track, then "Show all". It is always there, disabled with the reason while no track is hidden, so
// the host's layout never shifts on a hide. Thin over `viewer.hiddenTracks` and `track.show`. Everything it adds is
// taken away by the function it returns.

export interface HiddenTracksPartOptions {
  /** Where the menu opens from the button (default `below`). */
  readonly menuSide?: ToolbarMenuSide;
}

/** Mounts the hidden tracks' list of a viewer into `element` (after whatever is in it); the function it returns takes it away. */
export function mountHiddenTracks(
  viewer: SeqViewer,
  element: HTMLElement,
  options: HiddenTracksPartOptions = {},
): () => void {
  const doc = element.ownerDocument;
  const host = doc.createElement('div');
  host.className = 'seqx-hidden-tracks';
  host.style.position = 'relative';
  element.append(host);
  const hidden = viewer.hiddenTracks;
  const face = derive(
    hidden,
    (list): ToolbarFace => ({
      icon: 'eye-off',
      title: list.length ? `Hidden tracks (${list.length}): show` : 'Hidden tracks',
    }),
    (a, b) => a.title === b.title,
  );
  const rows = derive(hidden, (list): readonly ToolbarMenuRow[] => {
    const sections = new Set(list.map((t) => t.section));
    const out: ToolbarMenuRow[] = list.map((t) => ({
      id: t.address,
      label: t.label,
      icon: 'eye',
      // sections told apart only when there is more than one
      ...(sections.size > 1 ? { group: t.section, heading: t.section } : {}),
      run: () => void viewer.execute('track.show', { address: t.address }),
    }));
    if (list.length > 1)
      out.push({
        id: 'show-all',
        label: 'Show all',
        group: '\u0001all',
        run: () => void viewer.execute('track.show', { all: true }),
      });
    return out;
  });
  const item: ToolbarItem = {
    id: 'hidden-tracks',
    icon: 'eye-off',
    title: 'Hidden tracks',
    face,
    disabled: derive(hidden, (list) => (list.length ? false : 'No hidden tracks')),
    menu: rows,
  };
  const dispose = render(() => {
    onCleanup(installTooltips(host));
    return Toolbar({
      groups: [{ id: 'hidden-tracks', label: 'Hidden tracks', items: [item] }],
      'aria-label': 'Hidden tracks',
      ...(options.menuSide ? { menuSide: options.menuSide } : {}),
    });
  }, host);
  // the count as a badge on the button (a face's title is read by a tooltip only): shown with a track hidden, so the
  // button looks different with and without; it takes no room (positioned over the button's corner)
  const badge = doc.createElement('span');
  badge.className = 'seqx-hidden-tracks-badge';
  badge.setAttribute('aria-hidden', 'true');
  badge.style.cssText =
    'position: absolute; top: 2px; right: 2px; min-width: 14px; height: 14px; padding: 0 3px; box-sizing: border-box; border-radius: 7px; font: 600 10px/14px var(--msui-font-sans, sans-serif); text-align: center; pointer-events: none; background: var(--msui-accent, #2563eb); color: var(--msui-on-accent, #fff);';
  host.append(badge);
  const showCount = (list: readonly unknown[]) => {
    badge.textContent = list.length ? String(list.length) : '';
    badge.style.display = list.length ? '' : 'none';
    host.dataset.count = String(list.length);
  };
  showCount(hidden.value);
  const stop = hidden.listen(showCount);
  return () => {
    stop();
    dispose();
    host.remove();
  };
}
