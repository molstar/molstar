/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type State } from '@molstar/core/reactive';
import { type Child, onCleanup } from '../../dom/index.js';
import type { ToolbarGroup, ToolbarMenuRow } from '../../components/toolbar/toolbar.js';
import { installTooltips } from '../../components/tooltip/tooltip.js';
import { state } from '../../state.js';
import { BandPlaceholder, structureLines, residueLines, StagePlaceholder } from '../placeholders.js';
import type { StageModel } from '../stage-model.js';

// What the stage's stories share: a frame with tooltips, the placeholders wired to a model, a demo toolbar.

/** A story's box for a stage: the page's size less a margin, or a fixed one. */
export function StageBox(props: { children: Child; width?: string; height?: string }) {
  return (
    <div
      style={{
        position: 'relative',
        display: 'flex',
        width: props.width ?? '100%',
        height: props.height ?? 'calc(100vh - 2 * var(--msui-space-4))',
        margin: 'var(--msui-space-4) auto',
        maxWidth: props.width ? undefined : 'calc(100% - 2 * var(--msui-space-4))',
        borderRadius: 'var(--msui-radius-lg)',
        overflow: 'hidden',
        boxShadow: 'var(--msui-edge)',
      }}
      ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
    >
      {props.children}
    </div>
  );
}

/** The canvas stand-in, feeding the model's hover lines. */
export function DemoCanvas(props: { stage: StageModel }) {
  return <StagePlaceholder onHover={(i) => props.stage.setHover(i === undefined ? [] : structureLines(i))} />;
}

/** The band stand-in, feeding the model's hover lines. */
export function DemoBand(props: { stage: StageModel; rows?: number }) {
  return (
    <BandPlaceholder
      rows={props.rows ?? 2}
      length={160}
      selection={[11, 18]}
      onHover={(hit) => props.stage.setHover(hit ? residueLines(hit.row, hit.index) : [])}
    />
  );
}

/** A quiet pill for a slot's content in the stories (the status components are their own). */
export function Pill(props: { children: Child }) {
  return (
    <span
      style={{
        display: 'inline-flex',
        alignItems: 'center',
        height: 'var(--msui-control-height)',
        padding: '0 var(--msui-space-2)',
        borderRadius: 'var(--msui-radius-md)',
        background: 'var(--msui-surface-raised)',
        boxShadow: 'inset 0 0 0 1px var(--msui-line), var(--msui-shadow-float)',
        color: 'var(--msui-text-muted)',
        whiteSpace: 'nowrap',
      }}
    >
      {props.children}
    </span>
  );
}

function choice(
  value: State<string>,
  group: string,
  options: readonly (readonly [string, string])[],
): ToolbarMenuRow[] {
  return options.map(([id, label]) => ({
    id,
    label,
    group,
    choice: true,
    pressed: derive(value, (v) => v === id),
    run: () => value.set(id),
  }));
}

/** A viewer-like toolbar: the camera (an action, a toggle, a menu with its own action), lighting with a submenu, view. */
export function demoGroups(ran: (what: string) => void = () => {}): ToolbarGroup[] {
  const spin = state(false);
  const projection = state('perspective');
  const lighting = state('studio');
  const outline = state(false);
  const style = state('cartoon');
  return [
    {
      id: 'camera',
      label: 'Camera',
      items: [
        {
          id: 'reset',
          icon: 'refresh-cw',
          title: 'Reset the camera',
          run: () => ran('Reset the camera'),
          menu: choice(projection, 'projection', [
            ['perspective', 'Perspective'],
            ['orthographic', 'Orthographic'],
          ]),
        },
        { id: 'spin', icon: 'rotate-3d', title: 'Spin', pressed: spin, run: () => spin.set(!spin.value) },
        { id: 'screenshot', icon: 'camera', title: 'Screenshot', run: () => ran('Screenshot') },
      ],
    },
    {
      id: 'look',
      label: 'Look',
      items: [
        {
          id: 'lighting',
          icon: 'lightbulb',
          title: 'Lighting',
          menu: [
            ...choice(lighting, 'lighting', [
              ['flat', 'Flat'],
              ['studio', 'Studio'],
              ['soft', 'Soft shadows'],
            ]),
            {
              id: 'outline',
              label: 'Outline',
              group: 'extra',
              pressed: outline,
              run: () => outline.set(!outline.value),
            },
            {
              id: 'style',
              label: 'Style',
              icon: 'palette',
              group: 'extra',
              items: choice(style, 'style', [
                ['cartoon', 'Cartoon'],
                ['surface', 'Surface'],
                ['spacefill', 'Spacefill'],
              ]),
            },
          ],
        },
        { id: 'sequence', icon: 'panel-top', title: 'Sequence', run: () => ran('Sequence') },
      ],
    },
    {
      id: 'view',
      label: 'View',
      items: [
        { id: 'expand', icon: 'maximize', title: 'Expand', run: () => ran('Expand') },
        {
          id: 'fullscreen',
          icon: 'fullscreen',
          title: 'Fullscreen',
          disabled: 'The page does not allow it here',
          run: () => {},
        },
      ],
    },
  ];
}
