/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { VideoEncoderUI } from '@molstar/video-export-extension/ui';

export const VideoExport = PluginBehavior.create<{}>({
  // Preserve the serialized behavior id for existing snapshots.
  name: 'extension-mp4-export',
  category: 'misc',
  display: {
    name: 'Video Animation Export',
  },
  ctor: class extends PluginBehavior.Handler<{}> {
    register(): void {
      this.ctx.customStructureControls.set('video-export', VideoEncoderUI as any);
    }

    update() {
      return false;
    }

    unregister() {
      this.ctx.customStructureControls.delete('video-export');
    }
  },
  params: () => ({}),
});
