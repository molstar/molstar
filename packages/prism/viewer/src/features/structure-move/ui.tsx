/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Icon } from '@molstar/prism-ui';
import { For, Show } from '@molstar/prism-ui/dom/control-flow';
import { effect, onActive, onCleanup } from '@molstar/prism-ui/dom/scope';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { ViewerExtensionUI, ViewerUIContext } from '../../spec.js';
import { MoveSectionId } from './index.js';
import type { MoveApi } from './operations.js';

// The Move section of the Selection tab (structure-move.md): which structure the handles are on (the selection's, or
// the first; a row each when there are several), Reset; the handles show while the section does.

const caption = 'margin: 0; font-size: var(--msp-font-size-xs); color: var(--msp-color-text-muted)';

export function structureMoveUI(ctx: ViewerUIContext): ViewerExtensionUI {
  return {
    sections: [
      {
        id: MoveSectionId,
        slot: 'selection',
        title: 'Move',
        icon: 'move-3d',
        order: 0.7,
        keywords: ['move', 'rotate', 'turn', 'translate', 'transform', 'place', 'handles', 'gizmo'],
        render: domView(() => {
          const move = ctx.api.move;
          return move ? <MoveSection ctx={ctx} move={move} /> : <></>;
        }),
      },
    ],
  };
}

function MoveSection(props: { ctx: ViewerUIContext; move: MoveApi }) {
  const { ctx, move } = props;
  const gizmos = move.gizmos;
  // the structure the handles are on: kept while it is loaded; else the selection's, or the first
  const pick = () => {
    const list = move.structures.value;
    const now = gizmos.target.value;
    if (now && list.some((s) => s.ref === now)) return;
    gizmos.setTarget(move.selectedStructure() ?? list[0]?.ref);
  };
  let active = false;
  onActive((now) => {
    active = now;
    if (now) {
      pick();
      move.showGizmos(ctx.api);
    } else gizmos.hide();
  });
  // a structure loaded while the section is open (the first one), or the target removed: picked again
  effect(move.structures, () => {
    if (active) pick();
  });
  onCleanup(() => gizmos.hide());
  const many = derive(move.structures, (s) => s.length > 1);
  const target = gizmos.target;
  return (
    <div style="display: grid; gap: var(--msp-space-2); min-width: 0">
      <Show
        when={derive(move.structures, (s) => s.length > 0)}
        fallback={() => <p style={caption}>Nothing to move.</p>}
      >
        {() => (
          <>
            <p style={caption}>
              Drag the handles in 3D: an arrow moves it along its axis, a square in its plane, a ring turns it, the
              center moves it in the screen's plane. Esc cancels a drag.
            </p>
            <Show when={many}>
              {() => (
                <div style="display: grid; gap: 2px; min-width: 0" role="radiogroup" aria-label="The structure to move">
                  <For each={move.structures} by={(s) => s.ref}>
                    {(s) => (
                      <Button
                        size="xs"
                        variant="ghost"
                        block
                        pressed={derive(target, (t) => t === s.ref)}
                        title={`The handles on ${s.label}`}
                        onClick={() => gizmos.setTarget(s.ref)}
                      >
                        <span style="flex: 1; text-align: start">{s.label}</span>
                      </Button>
                    )}
                  </For>
                </div>
              )}
            </Show>
            <Button
              size="xs"
              variant="ghost"
              block
              icon={<Icon name="rotate-ccw" />}
              title="Back to where it was loaded"
              disabled={derive(target, (t) => !t)}
              onClick={() => {
                const t = target.value;
                if (t) void ctx.api.run('structure.transform-reset', { structure: t });
              }}
            >
              <span style="flex: 1; text-align: start">Reset its place</span>
            </Button>
          </>
        )}
      </Show>
    </div>
  );
}
