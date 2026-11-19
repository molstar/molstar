/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Patch } from '@molstar/core/repr/patch';
import type { Node } from '@molstar/core/repr/tree';
import { Data, entryLigands, ligand } from '@molstar/seqx-examples/complex/data';
import { hiddenData, ligandId } from '@molstar/seqx-examples/complex/ligands';
import { plddtColor } from '@molstar/seqx-examples/complex/confidence';
import { cellReadout, rectReadout } from '@molstar/seqx-examples/complex/readout';
import type { EditingBehavior } from '@molstar/seqx-viewer/behaviors/editing';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import type { EditingSession } from '../../features/editing/compose.js';

// The Complex's editing views, each an edit of one patch, one step of the host's stack (Undo and Redo walk them): the
// interface bookmarked as a ref; the confidence palette's bins; the ligands of entries past a resolution hidden (a
// filter in each layer's data); the copies merged or split. And the readout of the grid over the complex, the host's
// panel: the hovered cell and the selected rect named, the rect's PAE summarized. DOM-free: the stories and their
// headless twins call the same functions.

const editing = (s: EditingSession) => s.viewer.behaviors.get('editing') as EditingBehavior;

/** One patch as one step: a preview ended at once. */
async function step(s: EditingSession, label: string, patch: Patch): Promise<string> {
  if (!editing(s).preview(patch)) return 'Refused';
  await editing(s).endPreview(label);
  return label;
}

/** The spec's node at a path of ids (`['views', 'main', 'complex']`), as written. */
export function nodeAt(root: Node, path: readonly string[]): Node | undefined {
  let n: Node | undefined = root;
  for (const id of path) n = n?.children?.find((c) => c.id === id);
  return n;
}

/** Bookmarks 1FIN's interface (both sides' interface residues) as the ref `interface`; said when it is there. */
export function bookmarkInterface(s: EditingSession): Promise<string> {
  if (s.viewer.addresses.refs.value.some((r) => r.id === 'interface')) return Promise.resolve('Already bookmarked');
  return Promise.resolve(
    s.viewer.execute('refs.add', {
      ref: { id: 'interface', to: ['/data/cdk2-interface', '/data/cyclin-interface'], label: 'The interface (1FIN)' },
    }),
  );
}

/** The confidence view's pLDDT layers' bins: AlphaFold's edges, or others (one theme patch over both chains). */
export function setPlddtBins(s: EditingSession, edges: readonly number[]): Promise<string> {
  const ops = (['cdk2', 'cyclin'] as const).map((chain) => ({
    op: 'set' as const,
    at: ['views', 'main', chain, 'plddt', 'plddt'],
    key: 'color',
    value: plddtColor(edges),
  }));
  return step(s, 'Set the confidence bins', { ops });
}

/** The ligands whose entry's resolution is worse than `resolution` (Å) hidden: their layers' data filtered out. */
export function hideLigandsPast(s: EditingSession, resolution: number): Promise<string> {
  const all = [ligand('ATP'), ...entryLigands()];
  const ops = all.map((l) => {
    const pairs = `/data/${ligandId(l.comp)}-contacts`;
    return {
      op: 'set' as const,
      at: ['views', 'main', 'ligands', 'contacts', ligandId(l.comp)],
      key: 'data',
      value: l.resolution > resolution ? hiddenData(pairs) : pairs,
    };
  });
  return step(s, `Hide ligands past ${resolution} Å`, { ops });
}

/** The ligands shown by the spec as shown: those whose layer's data is not filtered out. */
export function shownLigands(s: EditingSession): string[] {
  const track = nodeAt(s.viewer.spec.shown.value, ['views', 'main', 'ligands', 'contacts']);
  return (track?.children ?? [])
    .filter((l) => typeof l.params?.data === 'string' && !(l.params.data as string).includes('(not'))
    .map((l) => l.id as string);
}

/** A section's copies merged (`×N`) or split, one patch. */
export function setCopies(s: EditingSession, section: string, copies: 'merge' | 'each'): Promise<string> {
  return step(s, copies === 'merge' ? 'Merge the copies' : 'Split the copies', {
    ops: [{ op: 'set', at: ['views', 'main', section], key: 'copies', value: copies }],
  });
}

/** What the readout panel shows: the hovered cell, the selected rect and its PAE, each a line ('' when none). */
export function readout(viewer: SeqViewer): { readonly hover: string; readonly rect: string } {
  const layout = viewer.frame?.layout;
  const model = viewer.spec.model.value;
  const view = viewer.view.resolved.value;
  if (!layout || !model || !view) return { hover: '', rect: '' };
  const h = viewer.hover.current;
  const v = { model, view };
  const hover =
    h.system !== undefined && h.ySystem !== undefined && h.position !== undefined && h.yPosition !== undefined
      ? cellReadout(v, {
          system: h.system,
          position: h.position,
          ySystem: h.ySystem,
          yPosition: h.yPosition,
        })
      : '';
  return { hover, rect: rectReadout(v, viewer.selection.value.value.location) };
}

/** The data the views hold, for a story's text (the chains, their copies). */
export const ComplexData = Data;
