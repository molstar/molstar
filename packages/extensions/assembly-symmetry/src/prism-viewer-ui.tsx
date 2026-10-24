/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Button, domView, Select, Show } from '@molstar/prism-ui';
import { mapState } from '@molstar/prism-ui/reactive/map-state';
import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import { type SceneNode, structureOfGroup } from '@molstar/prism-viewer/model/scene';
import { SceneRow } from '@molstar/prism-viewer/views/structure/scene-row';
import type { ViewerExtensionUI, ViewerUIContext } from '@molstar/prism-viewer/spec';
import type { SymmetryInfo, SymmetryRowState } from './operations.js';

// The assembly symmetry's row (parity.md §6; v5's Assembly Symmetry panel): under each structure that is a PDB entry's
// assembly, its eye shows and hides the symmetry axes (`assembly-symmetry.show`, `.hide`; the first that is not C1),
// opened it has the symmetries RCSB PDB lists (choosing one shows its axes) and Color by cluster
// (`assembly-symmetry.color`). The list shows once an operation has fetched it: the row never asks RCSB PDB by itself.

export const SymmetryRow = 'assembly-symmetry.row.symmetry';

const Quiet = 'margin: 0; color: var(--msp-color-text-muted); font-size: var(--msp-font-size-xs)';

const nameOf = (x: SymmetryInfo) =>
  `${x.symbol} · ${x.kind.replace(/ Symmetry$/, '').toLowerCase()} · ${x.stoichiometry.join(', ')}`;

/** The symmetry to show or color by: the chosen one, unless it is C1 or not known yet (0: the first that is not C1). */
const chosen = (s: SymmetryRowState) => {
  const x = s.symmetries?.find((y) => y.index === s.selected);
  return x && x.symbol !== 'C1' ? x.index : 0;
};

function describe(s: SymmetryRowState): string {
  const chosen = s.symmetries?.find((x) => x.index === s.selected);
  if (s.shown && chosen) return `${chosen.symbol} axes`;
  if (s.symmetries) return `${s.symmetries.length} symmetr${s.symmetries.length === 1 ? 'y' : 'ies'}`;
  return 'From RCSB PDB';
}

export function assemblySymmetryUI(ctx: ViewerUIContext): ViewerExtensionUI {
  const area = ctx.api['assembly-symmetry'];
  if (!area) return {};
  return {
    rows: [
      {
        id: SymmetryRow,
        // In the Other group of the structures it is for: a row whose `when` holds takes a slot in the
        // (virtualized) tree.
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
                  label="Assembly symmetry"
                  description={mapState(state, describe)}
                  shown={mapState(state, (s) => s.shown)}
                  onToggle={() =>
                    void (state.value.shown
                      ? area.hide({ structure })
                      : area.show({ structure, symmetry: chosen(state.value) }))
                  }
                  body={() => <SymmetryControls ctx={ctx} structure={structure} state={state} />}
                />
              )}
            </Show>
          );
        }),
      },
    ],
  };
}

function SymmetryControls(props: { ctx: ViewerUIContext; structure: string; state: ReadonlyState<SymmetryRowState> }) {
  const { structure, state } = props;
  const area = props.ctx.api['assembly-symmetry']!;
  // Only the symmetries with axes or clusters can be chosen (C1 has neither).
  const choices = mapState(state, (s) => s.symmetries?.filter((x) => x.symbol !== 'C1'));
  return (
    <>
      <Show
        when={choices}
        fallback={() => <p style={Quiet}>Show the axes or color by cluster to get the symmetries from RCSB PDB.</p>}
      >
        {(list) =>
          list.length ? (
            <Select
              size="sm"
              aria-label="Symmetry"
              options={list.map((x) => [String(x.index), nameOf(x)] as const)}
              value={mapState(state, (s) => String(chosen(s) || list[0].index))}
              onChange={(v) => void area.show({ structure, symmetry: Number(v) })}
            />
          ) : (
            <p style={Quiet}>RCSB PDB lists no symmetry but C1 for this assembly.</p>
          )
        }
      </Show>
      <div style="display: flex; gap: var(--msp-space-2)">
        <Button size="sm" variant="ghost" onClick={() => void area.color({ structure, symmetry: chosen(state.value) })}>
          Color by cluster
        </Button>
      </div>
    </>
  );
}
