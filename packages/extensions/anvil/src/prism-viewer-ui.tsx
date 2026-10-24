/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { bind, Button, domView, onCleanup, ParamsForm, Show } from '@molstar/prism-ui';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import { mapState } from '@molstar/prism-ui/reactive/map-state';
import { type SceneNode, structureOfGroup } from '@molstar/prism-viewer/model/scene';
import { SceneRow } from '@molstar/prism-viewer/views/structure/scene-row';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import type { MembraneRowState } from './operations.js';

// The membrane's row (parity.md §6; v5's generic representation control for ANVIL's membrane orientation): under each
// structure ANVIL applies to, its eye shows and hides the membrane (`membrane.show`, `.hide`), opened it has how the
// membrane looks (`membrane.options`) and Fit…, which opens `membrane.fit`'s form in the command menu (ANVIL's
// settings). It reads the row's state from the `membrane` area and acts through the person's API.

export const MembraneRow = 'anvil.row.membrane';

const Look = {
  planes: PD.Boolean(true, { label: 'Planes' }),
  rims: PD.Boolean(true, { label: 'Rims' }),
  dashed: PD.Boolean(false, { label: 'Dashed rims' }),
  opacity: PD.Numeric(0.5, { min: 0.05, max: 1, step: 0.05 }, { label: 'Opacity' }),
};

const describe = (s: MembraneRowState) =>
  !s.shown ? 'Off' : s.thickness !== undefined ? `${s.thickness} Å thick` : 'Shown';

export function membraneUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api.membrane;
  if (!area) return {};
  return {
    rows: [
      {
        id: MembraneRow,
        // In the Other group of each structure ANVIL applies to (a row whose `when` holds takes a slot).
        when: (node) => {
          const structure = node.kind === 'other' ? structureOfGroup(node) : undefined;
          return !!structure && area.applies(structure);
        },
        render: domView(({ ctx: view }) => {
          const structure = structureOfGroup(view.subject as SceneNode)!;
          const state = area.rowState(structure);
          return (
            <Show when={mapState(state, (s) => s.applicable)}>
              {() => (
                <SceneRow
                  label="Membrane"
                  description={mapState(state, describe)}
                  shown={mapState(state, (s) => s.shown)}
                  onToggle={() => void (state.value.shown ? area.hide({ structure }) : area.show({ structure }))}
                  body={() => <MembraneControls ctx={ctx} structure={structure} state={state} />}
                />
              )}
            </Show>
          );
        }),
      },
    ],
  };
}

function MembraneControls(props: {
  ctx: ViewerUIContext;
  structure: string;
  state: ReturnType<NonNullable<ViewerUIContext['api']['membrane']>['rowState']>;
}) {
  const { ctx, structure, state } = props;
  const area = ctx.api.membrane!;
  let timer: ReturnType<typeof setTimeout> | undefined;
  // A pending call must not outlive the controls: the row folded, the membrane hidden or the structure removed.
  onCleanup(() => clearTimeout(timer));
  const form = new ParamsFormModel(Look, PD.getDefaultValues(Look), {
    // The slider changes at pointer rate: one call a moment after the last.
    onChange: (v, name) => {
      clearTimeout(timer);
      const run = () =>
        void area.options({
          structure,
          ...(name === 'opacity' ? { opacity: v.opacity } : {}),
          ...(name === 'planes' ? { planes: v.planes } : {}),
          ...(name === 'rims' ? { rims: v.rims } : {}),
          ...(name === 'dashed' ? { dashed: v.dashed } : {}),
        });
      if (name === 'opacity') timer = setTimeout(run, 150);
      else run();
    },
  });
  // The form follows the membrane (Undo, a command, the agent), without reporting it as an edit.
  bind(state, (s) => s.look && form.load(s.look));
  return (
    <>
      <Show
        when={mapState(state, (s) => s.shown)}
        fallback={() => (
          <p style="margin: 0; color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)">
            Shown, it is placed by ANVIL from the residues' hydrophobicity (computed here, a moment for a large
            structure).
          </p>
        )}
      >
        {() => <ParamsForm form={form} />}
      </Show>
      <div style="display: flex; gap: var(--msp-space-2)">
        <Button size="sm" variant="ghost" onClick={() => void ctx.managers.commands.open('membrane.fit')}>
          Fit…
        </Button>
      </div>
    </>
  );
}
