/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { domView, onCleanup, ParamsForm, SectionList } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { IconName } from '@molstar/prism-ui/core';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { type CanvasSection, CanvasSections, sectionValues } from './operations.js';

// The Canvas card (parity.md §4; the mock's `registerMockCanvasSettings`): one section per group of Canvas3D's params,
// each a form over the canvas's props, applying through `canvas.settings` as the person; it follows the canvas when it
// changes elsewhere (Undo, a session, a command). Other extensions add sections to the same slot.

const Sections: readonly { id: CanvasSection; title: string; icon: IconName }[] = [
  { id: 'background', title: 'Background', icon: 'palette' },
  { id: 'lighting', title: 'Lighting & effects', icon: 'lightbulb' },
  { id: 'camera', title: 'Camera', icon: 'camera' },
  { id: 'clipping', title: 'Clipping', icon: 'scan' },
  { id: 'advanced', title: 'Advanced', icon: 'sliders-horizontal' },
];

export function canvasSettingsUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const form = (section: CanvasSection) => {
    const props = ctx.sources.canvas.value?.props;
    let timer: ReturnType<typeof setTimeout> | undefined;
    onCleanup(() => clearTimeout(timer));
    return new ParamsFormModel(
      CanvasSections[section].params,
      (props ? sectionValues(section, props) : undefined) as never,
      {
        onChange: (values) => {
          clearTimeout(timer);
          timer = setTimeout(() => void ctx.api.run('canvas.settings', { section, values }), 120);
        },
      },
    );
  };
  return {
    overlays: [
      {
        id: 'canvas-settings.overlay.canvas',
        title: 'Canvas',
        render: domView(({ ctx: view }) => (
          <SectionList shell={view.shell} slot="canvas-settings" focus={view.focus} />
        )),
      },
    ],
    sections: Sections.map((s, order) => ({
      id: `canvas-settings.section.${s.id}`,
      slot: 'canvas-settings' as const,
      title: s.title,
      icon: s.icon,
      order,
      render: domView(() => {
        const f = form(s.id);
        // The canvas changed elsewhere: the form follows, without reporting it as an edit.
        const sub = ctx.sources.canvas.subscribe((c) => c && f.load(sectionValues(s.id, c.props) as never));
        onCleanup(() => sub.unsubscribe());
        return <ParamsForm form={f} />;
      }),
    })),
    toolbar: [
      {
        id: 'canvas-settings.toolbar.canvas',
        icon: 'monitor-cog',
        title: 'Canvas settings',
        group: 'scene',
        overlay: 'canvas-settings.overlay.canvas',
      },
    ],
  };
}
