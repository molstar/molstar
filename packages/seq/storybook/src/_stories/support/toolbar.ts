/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mountHiddenTracks } from '@molstar/seqx-viewer/mount-hidden-tracks';
import { mountToolbar } from '@molstar/seqx-viewer/mount-toolbar';
import type { ToolbarMount } from '@molstar/seqx-viewer/viewer-spec';

// The default viewer's toolbar as the stories compose it: the viewer's own (back, forward, keys) and, beside it, the
// application's list of the hidden tracks (`mountHiddenTracks`): the view draws nothing for a hidden track, so the
// host that wants them reachable places the list (track-header-parts.md, "Hidden tracks are the application's"). The
// minimal viewer mounts none.

/** `ui.toolbar` for a story: `mountToolbar` with the hidden tracks' button in the same bar. */
export const toolbarWithHiddenTracks: ToolbarMount = (viewer, target) => {
  const off = mountToolbar(viewer, target);
  const bar = target.querySelector<HTMLElement>('.seqx-toolbar');
  if (!bar) return off;
  Object.assign(bar.style, { display: 'flex', alignItems: 'flex-start', gap: 'var(--msui-space-2, 8px)' });
  const slot = target.ownerDocument.createElement('div');
  bar.append(slot);
  const offList = mountHiddenTracks(viewer, slot);
  return () => {
    offList();
    slot.remove();
    off();
  };
};
