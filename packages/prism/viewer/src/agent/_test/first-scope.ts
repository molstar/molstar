/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureElement } from '@molstar/model/model/structure';
import type { PrismViewer } from '../../viewer.js';
import type { EvalCheck, EvalTask } from '../eval.js';

// The first eval (plan 4.5): tasks from the agent's first scope (README §3.3) on 1CBS, offline (the file is loaded
// first), each checked on the scene. The unit test runs them with a scripted model; `evals.live.ts` with real ones.

const Red = 0xff0000;

const structures = (v: PrismViewer) => v.plugin.managers.structure.hierarchy.current.structures;
const components = (v: PrismViewer) => structures(v).flatMap((s) => s.components);
const representationTypes = (v: PrismViewer) =>
  components(v).flatMap((c) =>
    c.representations.map((r) => (r.cell.transform.params as { type: { name: string } }).type.name),
  );
const label = (c: ReturnType<typeof components>[number]) => c.cell.obj?.label ?? '';

const check = (name: string, pass: boolean, detail: () => string): EvalCheck =>
  pass ? { name, pass } : { name, pass, detail: detail() };

/** Whether a red is set on the ligand: its overpaint layers, or its representations' color theme. */
function ligandIsRed(v: PrismViewer): boolean {
  const ligand = components(v).find((c) => /ligand/i.test(label(c)));
  if (!ligand) return false;
  const cells = [...v.plugin.state.data.cells.values()].filter(
    (cell) =>
      cell.transform.ref === ligand.cell.transform.ref || isUnder(v, cell.transform.ref, ligand.cell.transform.ref),
  );
  // The ligand's representations, and the overpaint layers on them (`color` on the target "ligand").
  const red = new RegExp(`\\b${Red}\\b`);
  return cells.some((cell) => red.test(JSON.stringify(cell.transform.params)));
}

function isUnder(v: PrismViewer, ref: string, ancestor: string): boolean {
  const tree = v.plugin.state.data.tree;
  let at = tree.transforms.get(ref);
  while (at && at.parent !== at.ref) {
    if (at.parent === ancestor) return true;
    at = tree.transforms.get(at.parent);
  }
  return false;
}

function selectedResidues(v: PrismViewer): number {
  let n = 0;
  for (const s of structures(v)) {
    const data = s.cell.obj?.data;
    if (!data) continue;
    const loci = v.plugin.managers.structure.selection.getLoci(data);
    if (StructureElement.Loci.is(loci)) n += StructureElement.Stats.ofLoci(loci).residueCount;
  }
  return n;
}

export const FirstScope: readonly EvalTask[] = [
  {
    id: 'surface',
    title: 'A style: the protein as a surface',
    files: ['1cbs_updated.cif'],
    prompts: ['Show the protein as a molecular surface.'],
    undoes: true,
    check: (v) => [
      check(
        'a surface is shown',
        representationTypes(v).some((t) => /surface/.test(t)),
        () => representationTypes(v).join(', '),
      ),
    ],
  },
  {
    id: 'ligand-red',
    title: 'A color on a target: the ligand red',
    files: ['1cbs_updated.cif'],
    prompts: ['Color the ligand red.'],
    undoes: true,
    check: (v) => [check('the ligand is red', ligandIsRed(v), () => 'no red on the ligand')],
  },
  {
    id: 'select-range',
    title: 'A selection by residue range',
    files: ['1cbs_updated.cif'],
    prompts: ['Select residues 10 to 20 of chain A.'],
    check: (v) => {
      const n = selectedResidues(v);
      return [check('11 residues are selected', n === 11, () => `${n} selected`)];
    },
  },
  {
    id: 'remove-water',
    title: 'A destructive change, allowed: the water removed',
    files: ['1cbs_updated.cif'],
    prompts: ['Remove the water from the scene.'],
    undoes: true,
    check: (v) => [
      check('no water component', !components(v).some((c) => /water/i.test(label(c))), () =>
        components(v).map(label).join(', '),
      ),
      check(
        'the protein stays',
        components(v).some((c) => /polymer|protein/i.test(label(c))),
        () => components(v).map(label).join(', '),
      ),
    ],
  },
];
