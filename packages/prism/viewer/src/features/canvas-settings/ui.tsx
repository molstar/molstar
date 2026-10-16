/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { bind, domView, onCleanup, ParamsForm, SectionList } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { IconName } from '@molstar/prism-ui/core';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { debounce } from '../../ui/debounce.js';
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
  // A section's form: edits wait for the debounce, then go through `canvas.settings`; what the canvas has is loaded
  // into the form unless the person's edit is still on its way (waiting, or sent and not yet applied), so a drag is
  // not pulled back to the value of a moment ago.
  const form = (section: CanvasSection) => {
    let sending = 0;
    const show = () => {
      const props = ctx.sources.canvas.value?.props;
      if (props) f.load(sectionValues(section, props) as never);
    };
    const send = debounce((values: Record<string, unknown>) => {
      sending++;
      void ctx.api.run('canvas.settings', { section, values }).finally(() => {
        sending--;
        if (!sending && !send.pending) show();
      });
    });
    onCleanup(() => send.dispose());
    const props = ctx.sources.canvas.value?.props;
    const f = new ParamsFormModel(
      CanvasSections[section].params,
      (props ? sectionValues(section, props) : undefined) as never,
      { onChange: (values) => send(values) },
    );
    return { form: f, editing: () => send.pending || sending > 0, show };
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
        // Bound, not subscribed: a folded section holds no subscription (prism AGENTS.md, hidden UI).
        bind(ctx.sources.canvas, (c) => c && !f.editing() && f.form.load(sectionValues(s.id, c.props) as never));
        return <ParamsForm form={f.form} />;
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
