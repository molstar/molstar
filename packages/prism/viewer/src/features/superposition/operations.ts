/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '@molstar/core/math/linear-algebra';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { QueryContext, type Structure, StructureElement, StructureSelection } from '@molstar/model/model/structure';
import { alignAndSuperpose, superpose } from '@molstar/model/model/structure/structure/util/superposition';
import {
  DefaultLigandMccsOptions,
  isSingleLigandLoci,
  superposeLigandsByMccs,
} from '@molstar/model/model/structure/structure/util/superposition-ligand';
import { alignAndSuperposeWithSIFTSMapping } from '@molstar/model/model/structure/structure/util/superposition-sifts-mapping';
import { tmAlign } from '@molstar/model/model/structure/structure/util/tm-align';
import { SIFTSMapping } from '@molstar/model/props/sequence/sifts-mapping';
import { trace } from '@molstar/plugin/state/queries/structure/structure';
import { polymer } from '@molstar/plugin/state/queries/structure/type';
import type { StructureRef as HierarchyStructure } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { TransformStructureConformation } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { StateTransform } from '@molstar/core/state';
import type { PluginContext } from '@molstar/plugin/context';
import type { ReadonlyState } from '@molstar/prism-ui/core';
import { type ApiAreaContext, ApiArea, type ApiContext, type ApiFailure, ApiResult } from '../../api/area.js';
import { OptionalParam } from '../../api/params.js';
import { isFailure, labelOf, loadedStructures, pickNamed } from '../../api/scope.js';
import type { ResolvedPart } from '../../api/targets.js';
import { LanguageParam, type TargetValue } from '../../api/target-param.js';
import { structureNames } from '../../sources/naming.js';
import { polymerPart } from '../../sources/selection.js';
import { sameItems } from '../../sources/state.js';
import { watchedState } from '../../sources/watched.js';
import type { QueryLanguage, Scoped } from '../../targets/types.js';

// The `superpose` area (api.md §2, state.md: "the algorithms move out of the view into the API"): each operation
// superposes structures onto a reference (the first) and records the transforms as one history step, so Undo puts them
// back. The algorithms are the classic Superposition panel's, DOM-free; the transform is Mol*'s own
// `TransformStructureConformation` tagged `SuperpositionTransform`, so the panel and this API leave the same state.
//
// Targets: one per structure, in order; the first is the reference. A target is any target (targets.md §1) or a
// `{ target, structure }` that names its structure; a bare one applies to the structure at its position among the
// loaded ones. No targets: the selection, one structure each (the structures that have one). Positions come from the
// structures as loaded (their untransformed coordinates), so superposing again with the same inputs gives the same
// transforms, which is "already so" (no history step).

const SuperpositionTag = 'SuperpositionTransform';

/** Fewest atoms that fix an orientation; with fewer the rotation is not determined. */
const MinAtoms = 3;

/** Transforms equal within this (Å for the translation) are the same: re-running reproduces them to ~1e-9. */
const Epsilon = 1e-4;

export type SuperposeMethod = 'chains' | 'tm-align' | 'ligands' | 'atoms' | 'db';

/** One target of a superposition: a target, or a target with the structure it is on. */
export type SuperposeTarget = TargetValue | Scoped;

/** The options of the chain-based methods (the classic panel's `StructureSuperpositionParams`). */
export const SuperpositionOptions = {
  alignSequences: PD.Boolean(true, {
    description:
      'Perform a sequence alignment of the chains and use the aligned residue pairs to guide the 3D superposition.',
  }),
  traceOnly: PD.Boolean(true, {
    description: 'Base the superposition only on CA (and equivalent) atoms.',
  }),
};

const TargetsParam = (what: string, extra = '') =>
  PD.Value<readonly SuperposeTarget[]>([], {
    description: `${what}, one per structure; the first is the reference. A target, or { target, structure } to name its structure (else the one at its position among the loaded). Empty: the selection, one per structure that has one.${extra}`,
  });

/** What one superposed structure got. */
export interface SuperposedStructure {
  readonly structure: { readonly ref: string; readonly label: string };
  /** Å, over the atoms the method paired. */
  readonly rmsd: number;
  /** How many atoms (or residues, for TM-align) were paired. */
  readonly paired: number;
  /** TM-align: the score normalized by the reference's length, and by this structure's. */
  readonly tmScore?: number;
  readonly tmScoreMobile?: number;
  /** Ligands: how they were matched. */
  readonly matched?: 'atom-name' | 'mccs';
  /** Its transform changed (false: it already was there). */
  readonly moved: boolean;
}

/** What a superposition gives back. */
export interface SuperposeValue {
  readonly method: SuperposeMethod;
  readonly reference: { readonly ref: string; readonly label: string };
  readonly results: readonly SuperposedStructure[];
}

// --- the targets ------------------------------------------------------------------------------------------------

/** One resolved target: the structure it is on and the loci. */
interface Member {
  readonly structure: HierarchyStructure;
  readonly data: Structure;
  readonly label: string;
  readonly loci: StructureElement.Loci;
  /** What the target selected; absent for a whole structure (by UniProt). */
  readonly part?: ResolvedPart;
}

/** What a method needs of each target, for the messages ("a chain", "a ligand"). */
interface Wanted {
  readonly what: string;
  readonly hint: string;
}

const scoped = (t: SuperposeTarget): Scoped | undefined =>
  t && typeof t === 'object' && 'target' in t && !('loci' in t) ? (t as Scoped) : undefined;

/** The structures that have a selection, in order. */
function selectedStructures(ctx: ApiContext): HierarchyStructure[] {
  const selection = ctx.plugin.managers.structure.selection;
  return loadedStructures(ctx).filter((s) => {
    const data = s.cell.obj?.data;
    if (!data) return false;
    const loci = selection.getLoci(data);
    return StructureElement.Loci.is(loci) && !StructureElement.Loci.isEmpty(loci);
  }) as HierarchyStructure[];
}

async function members(
  ctx: ApiContext,
  targets: readonly SuperposeTarget[],
  wanted: Wanted,
  language?: QueryLanguage,
): Promise<Member[] | ApiFailure> {
  const all = loadedStructures(ctx);
  let items: readonly SuperposeTarget[] = targets;
  if (!items.length) {
    items = selectedStructures(ctx).map((s) => ({ target: 'selection', structure: s.cell.transform.ref }));
    if (items.length < 2) {
      return ApiResult.fail(
        'invalid',
        `Superposing needs ${wanted.what} in each of two or more structures; ${items.length ? `only ${labelOf(ctx, selectedStructures(ctx)[0])} has a selection` : 'nothing is selected'}.`,
        wanted.hint,
      );
    }
  }
  if (items.length < 2) {
    return ApiResult.fail(
      'invalid',
      'Superposing needs a target in each of two or more structures.',
      `Give one per structure: ${wanted.hint}`,
    );
  }
  if (items.length > all.length) {
    return ApiResult.fail(
      'invalid',
      `${items.length} targets for ${all.length} loaded structure${all.length === 1 ? '' : 's'}.`,
      'Give one target per structure.',
    );
  }
  const out: Member[] = [];
  for (let i = 0; i < items.length; i++) {
    const s = scoped(items[i]);
    const target = s ? s.target : (items[i] as TargetValue);
    const structure = s?.structure || all[i].cell.transform.ref;
    const resolved = await ctx.targets.resolve(target, { structure, ...(language ? { language } : {}) });
    if (isFailure(resolved)) return resolved;
    const part = resolved.parts.find((p) => p.selected.atoms > 0);
    if (!part) {
      return ApiResult.fail(
        'invalid',
        `Target ${i + 1} selects nothing in ${resolved.parts[0]?.structure.label ?? structure}.`,
        wanted.hint,
      );
    }
    const found = all.find((x) => x.cell.transform.ref === part.structure.ref);
    const data = found?.cell.obj?.data;
    if (!found || !data) return ApiResult.fail('not-found', `There is no structure '${part.structure.label}'.`);
    if (out.some((m) => m.structure === found)) {
      return ApiResult.fail(
        'invalid',
        `Targets ${out.findIndex((m) => m.structure === found) + 1} and ${i + 1} are both in ${part.structure.label}.`,
        'Superposing takes one target per structure.',
      );
    }
    out.push({
      structure: found,
      data,
      label: labelOf(ctx, found),
      loci: StructureElement.Loci.remap(part.loci, data),
      part,
    });
  }
  return out;
}

// --- the transforms ---------------------------------------------------------------------------------------------

/** The transform a structure has now: its conformation transform's matrix, else none (identity). */
function currentMatrix(s: HierarchyStructure): Mat4 {
  const t = s.transform?.cell;
  const params = t?.transform.params?.transform as
    | { name: string; params: { data?: Mat4; transpose?: boolean } }
    | undefined;
  if (params?.name === 'matrix' && params.params.data) {
    const m = Mat4.clone(params.params.data);
    return params.params.transpose ? Mat4.transpose(m, m) : m;
  }
  const cs = t?.obj?.data.coordinateSystem;
  return cs ? Mat4.clone(cs.matrix) : Mat4.identity();
}

/** What the superposed structures were told to be, before they are written. */
interface Move {
  readonly member: Member;
  readonly matrix: Mat4;
  readonly rmsd: number;
  readonly paired: number;
  readonly extra: Partial<SuperposedStructure>;
}

/**
 * What a superposition records on the transform it writes (tags, so it follows the state: Undo, sessions): what the
 * structure was superposed onto, how, and the RMSD. Mol*'s own panel writes the first tag only.
 */
const Recorded = { onto: 'superposed-onto:', method: 'superposed-by:', rmsd: 'superposed-rmsd:' } as const;

function recordTags(reference: string, method: SuperposeMethod, rmsd: number, kept: readonly string[] = []): string[] {
  return [
    ...kept.filter((t) => t !== SuperpositionTag && !Object.values(Recorded).some((p) => t.startsWith(p))),
    SuperpositionTag,
    `${Recorded.onto}${reference}`,
    `${Recorded.method}${method}`,
    `${Recorded.rmsd}${Number.isFinite(rmsd) ? rmsd.toFixed(3) : 'NaN'}`,
  ];
}

/** Writes the moves as one tree update (an existing conformation transform is updated, else one is inserted). */
async function apply(
  ctx: ApiContext,
  method: SuperposeMethod,
  reference: Member,
  moves: readonly Move[],
): Promise<SuperposedStructure[]> {
  // The reference may have been moved itself: its transform leads, as the classic panel's `coordinateSystem`.
  const cs = reference.structure.transform?.cell.obj?.data.coordinateSystem;
  const lead = cs && !Mat4.isIdentity(cs.matrix) ? cs.matrix : undefined;
  const referenceRef = reference.structure.cell.transform.ref;
  const state = ctx.plugin.state.data;
  const build = state.build();
  let any = false;
  const out: SuperposedStructure[] = [];
  for (const m of moves) {
    const matrix = lead ? Mat4.mul(Mat4(), lead, m.matrix) : m.matrix;
    const moved = !Mat4.areEqual(currentMatrix(m.member.structure), matrix, Epsilon);
    out.push({
      structure: { ref: m.member.structure.cell.transform.ref, label: m.member.label },
      rmsd: m.rmsd,
      paired: m.paired,
      ...m.extra,
      moved,
    });
    if (!moved) continue;
    any = true;
    const params = { transform: { name: 'matrix' as const, params: { data: Mat4.clone(matrix), transpose: false } } };
    const existing = m.member.structure.transform?.cell;
    if (existing) {
      build.to(existing).update(params);
      // The record replaces an earlier one (the builder has no tag removal: its tree does).
      build.currentTree.setTags(
        existing.transform.ref,
        recordTags(referenceRef, method, m.rmsd, existing.transform.tags),
      );
    } else {
      build
        .to(m.member.structure.cell)
        .insert(TransformStructureConformation, params, { tags: recordTags(referenceRef, method, m.rmsd) });
    }
  }
  if (any) await ctx.plugin.runTask(state.updateTree(build));
  return out;
}

// --- the superpositions made ------------------------------------------------------------------------------------

/** A structure that a superposition moved, as its transform records it (`superpose.list`, the Superpose section). */
export interface SuperposedInfo {
  readonly structure: { readonly ref: string; readonly label: string };
  /** What it was superposed onto; undefined when not recorded (Mol*'s own panel) or no longer loaded. */
  readonly reference?: { readonly ref: string; readonly label: string };
  readonly method?: SuperposeMethod;
  readonly rmsd?: number;
}

/** Each method as a person says it. */
export const SuperposeMethodNames: Readonly<Record<SuperposeMethod, string>> = {
  chains: 'chains',
  'tm-align': 'TM-align',
  ligands: 'ligands',
  atoms: 'atoms',
  db: 'UniProt',
};

/** The superposition in words: "1W0E onto 1TQN by chains, RMSD 0.78 Å". */
export function superposedLine(s: SuperposedInfo): string {
  const how = [s.reference ? `onto ${s.reference.label}` : '', s.method ? `by ${SuperposeMethodNames[s.method]}` : '']
    .filter(Boolean)
    .join(' ');
  const rmsd = s.rmsd !== undefined && Number.isFinite(s.rmsd) ? `, RMSD ${fixed(s.rmsd)} Å` : '';
  return `${s.structure.label} ${how || 'superposed'}${rmsd}`;
}

/**
 * The structures a superposition moved and that are still where it put them: their conformation transform carries the
 * superposition's tag and is not the identity (a reset puts them back where they were loaded).
 */
export function readSuperposed(plugin: PluginContext): SuperposedInfo[] {
  const structures = plugin.managers.structure.hierarchy.current.structures;
  const names = structureNames(structures);
  const out: SuperposedInfo[] = [];
  for (const s of structures) {
    const t = s.transform?.cell.transform;
    if (!t || !StateTransform.hasTag(t, SuperpositionTag)) continue;
    if (Mat4.isIdentity(currentMatrix(s), Epsilon)) continue;
    const tag = (prefix: string) => t.tags?.find((x) => x.startsWith(prefix))?.slice(prefix.length);
    const ref = s.cell.transform.ref;
    const onto = tag(Recorded.onto);
    const method = tag(Recorded.method) as SuperposeMethod | undefined;
    const rmsd = tag(Recorded.rmsd);
    out.push({
      structure: { ref, label: names.get(ref) ?? ref },
      ...(onto && names.has(onto) ? { reference: { ref: onto, label: names.get(onto)! } } : {}),
      ...(method && method in SuperposeMethodNames ? { method } : {}),
      ...(rmsd !== undefined && rmsd !== 'NaN' ? { rmsd: Number(rmsd) } : {}),
    });
  }
  return out;
}

const sameSuperposed = (a: SuperposedInfo, b: SuperposedInfo) =>
  a.structure.ref === b.structure.ref &&
  a.structure.label === b.structure.label &&
  a.reference?.ref === b.reference?.ref &&
  a.reference?.label === b.reference?.label &&
  a.method === b.method &&
  a.rmsd === b.rmsd;

const fixed = (n: number) => (Number.isFinite(n) ? n.toFixed(2) : String(n));

/** The result of a method that moved structures onto a reference. */
function finish(
  ctx: ApiContext,
  method: SuperposeMethod,
  reference: Member,
  results: readonly SuperposedStructure[],
  how: string,
  warnings: string[] = [],
) {
  const value: SuperposeValue = {
    method,
    reference: { ref: reference.structure.cell.transform.ref, label: reference.label },
    results,
  };
  const line = results
    .map((r) => {
      const score = r.tmScore !== undefined ? `TM-score ${r.tmScore.toFixed(3)}, ` : '';
      return `${r.structure.label} (${score}RMSD ${fixed(r.rmsd)} Å over ${r.paired} ${method === 'atoms' || method === 'ligands' ? 'atoms' : 'pairs'})`;
    })
    .join('; ');
  const moved = results.filter((r) => r.moved);
  if (!moved.length) {
    return ApiResult.already(`Already superposed onto ${reference.label} ${how}: ${line}.`, { value });
  }
  // A canvas shows the new arrangement whole.
  if (ctx.modules.has('canvas')) ctx.plugin.managers.camera.reset();
  return ApiResult.ok(
    `Superposed ${results.map((r) => r.structure.label).join(', ')} onto ${reference.label} ${how}: ${line}.`,
    { nodes: moved.map((r) => r.structure.ref), view: true },
    { value, ...(warnings.length ? { warnings } : {}) },
  );
}

// --- the methods ------------------------------------------------------------------------------------------------

/**
 * The polymer chain of a target, as the atoms the method uses (CA and equivalent, or all polymer atoms). The classic
 * panel's rule is one polymer chain per structure (or residues within it); a target that holds more than the chain
 * ("chain A" with its waters and ligands) gives its polymer part, which must be one chain.
 */
function chainLoci(ms: readonly Member[], traceOnly: boolean): StructureElement.Loci[] | ApiFailure {
  const { query } = traceOnly ? trace : polymer;
  const out: StructureElement.Loci[] = [];
  for (const m of ms) {
    const part = StructureElement.Loci.remap(polymerPart(m.loci), m.data);
    if (StructureElement.Loci.isEmpty(part)) {
      return ApiResult.fail(
        'invalid',
        `The target in ${m.label} has no polymer (${m.part?.selected.line}).`,
        'Select one polymer chain, or residues within it, in each structure.',
      );
    }
    if (part.elements.length > 1) {
      return ApiResult.fail(
        'invalid',
        `The target in ${m.label} is not one polymer chain: its polymer is in ${part.elements.length} chains or copies (${m.part?.selected.line}).`,
        'Select one polymer chain, or residues within it, in each structure.',
      );
    }
    const sub = StructureElement.Loci.toStructure(part);
    const loci = StructureElement.Loci.remap(
      StructureSelection.toLociWithSourceUnits(query(new QueryContext(sub))),
      m.data,
    );
    const size = StructureElement.Loci.size(loci);
    if (size < MinAtoms) {
      return ApiResult.fail(
        'invalid',
        `The target in ${m.label} has ${size} ${traceOnly ? 'trace ' : ''}atom${size === 1 ? '' : 's'}: too few to superpose (${MinAtoms} are needed).`,
        'Select more of the chain.',
      );
    }
    out.push(loci);
  }
  return out;
}

const noRmsd = (label: string) =>
  ApiResult.fail(
    'failed',
    `The superposition onto ${label} did not converge.`,
    'Check that the chains have atoms in common.',
  );

const needsTwo = (ctx: ApiAreaContext) =>
  ctx.sources.structures.value.length < 2 ? 'Superposing needs two or more structures.' : (true as const);

/** The structures that cannot be placed by UniProt, as a failure. */
function noMapping(ctx: ApiContext, lacking: readonly HierarchyStructure[]) {
  return ApiResult.fail(
    'invalid',
    `${lacking.map((s) => labelOf(ctx, s)).join(', ')} ${lacking.length === 1 ? 'has' : 'have'} no SIFTS UniProt mapping.`,
    'The mapping comes with the updated mmCIF files of the PDB (pdbx_sifts_xref_db_* in atom_site); use chains or TM-align instead.',
  );
}

export class SuperposeApi extends ApiArea {
  readonly id = 'superpose';

  #made: ReadonlyState<readonly SuperposedInfo[]> | undefined;

  /**
   * The structures superposed now (`readSuperposed`), following the state while watched (the Superpose section's list;
   * prism AGENTS.md, hidden UI: nothing listens while nobody looks), read once per animation frame.
   */
  get made(): ReadonlyState<readonly SuperposedInfo[]> {
    const plugin = this.ctx.plugin;
    this.#made ??= watchedState(
      () => readSuperposed(plugin),
      () => [plugin.state.data.events.changed, plugin.managers.structure.hierarchy.behaviors.selection],
      { compare: (a, b) => sameItems(a, b, sameSuperposed) },
    );
    return this.#made;
  }

  readonly list = this.operation({
    id: 'superpose.list',
    title: 'List the superpositions',
    params: {},
    read: true,
    async run(_, ctx) {
      const list = readSuperposed(ctx.plugin);
      if (!list.length) return ApiResult.ok('No structure is superposed.', {}, { value: list });
      return ApiResult.ok(list.map(superposedLine).join('\n'), {}, { value: list });
    },
  });

  // Undoes a superposition where it stands (not by history): the structure goes back where it was loaded (a
  // superposition replaces any earlier transform, so there is nothing else to restore) and its record goes.
  readonly remove = this.operation({
    id: 'superpose.remove',
    title: 'Remove a superposition',
    params: {
      structure: PD.Text('', {
        description: 'A superposed structure by label or state ref; empty: the only one (several: name it, or all)',
      }),
      all: OptionalParam(PD.Boolean(false, { description: 'Every superposed structure' })),
    },
    when: (ctx) => (readSuperposed(ctx.plugin).length ? true : 'No structure is superposed.'),
    async run({ structure, all: every }, ctx) {
      const made = readSuperposed(ctx.plugin);
      if (!made.length) return ApiResult.already('No structure is superposed.');
      let doomed = made;
      // An empty remove is the only one (api.md §1 "Words"); `all` clears.
      if (!structure && !every && made.length > 1) {
        return ApiResult.fail(
          'ambiguous',
          `Which superposition? ${made.length} structures are superposed.`,
          `Name one (${made.map((m) => m.structure.label).join(', ')}), or give all: true.`,
        );
      }
      if (structure) {
        const all = loadedStructures(ctx);
        const s = pickNamed(all, structure, {
          what: 'structure',
          ids: (x) => [x.cell.transform.ref],
          names: (x) => [labelOf(ctx, x), x.cell.obj?.label],
          describe: (x) => `${labelOf(ctx, x)} (${x.cell.transform.ref})`,
          none: { message: 'Nothing is loaded.', hint: 'Load structures first.' },
        });
        if (isFailure(s)) return s;
        doomed = made.filter((m) => m.structure.ref === s.cell.transform.ref);
        if (!doomed.length) {
          return ApiResult.already(`${labelOf(ctx, s)} is not superposed.`);
        }
      }
      const identity = { transform: { name: 'matrix' as const, params: { data: Mat4.identity(), transpose: false } } };
      const state = ctx.plugin.state.data;
      const build = state.build();
      for (const m of doomed) {
        const cell = loadedStructures(ctx).find((x) => x.cell.transform.ref === m.structure.ref)?.transform?.cell;
        if (!cell) continue;
        build.to(cell).update(identity);
        build.currentTree.setTags(
          cell.transform.ref,
          (cell.transform.tags ?? []).filter(
            (t) => t !== SuperpositionTag && !Object.values(Recorded).some((p) => t.startsWith(p)),
          ),
        );
      }
      await ctx.plugin.runTask(state.updateTree(build));
      if (ctx.modules.has('canvas')) ctx.plugin.managers.camera.reset();
      const names = doomed.map((m) => m.structure.label).join(', ');
      return ApiResult.ok(
        `Put ${names} back where ${doomed.length === 1 ? 'it was' : 'they were'} loaded.`,
        { nodes: doomed.map((m) => m.structure.ref), view: true },
        { value: doomed },
      );
    },
  });

  readonly chains = this.operation({
    id: 'superpose.chains',
    title: 'Superpose by chains',
    params: {
      language: LanguageParam(this.ctx.languages),
      targets: TargetsParam('The polymer chain (or residues within it) to superpose'),
      ...SuperpositionOptions,
    },
    when: needsTwo,
    async run({ targets, alignSequences, traceOnly, language }, ctx) {
      const ms = await members(
        ctx,
        targets,
        {
          what: 'a polymer chain',
          hint: 'Select one polymer chain (or residues within it) in each of two or more structures.',
        },
        language,
      );
      if (isFailure(ms)) return ms;
      const locis = chainLoci(ms, traceOnly);
      if (isFailure(locis)) return locis;
      const transforms = alignSequences ? alignAndSuperpose(locis) : superpose(locis);
      const moves: Move[] = [];
      for (let i = 1; i < ms.length; i++) {
        const { bTransform, rmsd, nAlignedElements } = transforms[i - 1];
        if (Number.isNaN(rmsd)) return noRmsd(ms[0].label);
        moves.push({ member: ms[i], matrix: bTransform, rmsd, paired: nAlignedElements, extra: {} });
      }
      const results = await apply(ctx, 'chains', ms[0], moves);
      return finish(ctx, 'chains', ms[0], results, `by chains${alignSequences ? ' (sequence-aligned)' : ''}`);
    },
  });

  readonly tmAlign = this.operation({
    id: 'superpose.tm-align',
    title: 'Superpose by TM-align',
    params: {
      language: LanguageParam(this.ctx.languages),
      targets: TargetsParam(
        'The polymer chain to superpose',
        ' Sequence-independent: for chains of low sequence identity.',
      ),
      traceOnly: SuperpositionOptions.traceOnly,
    },
    when: needsTwo,
    async run({ targets, traceOnly, language }, ctx) {
      const ms = await members(
        ctx,
        targets,
        {
          what: 'a polymer chain',
          hint: 'Select one polymer chain in each of two or more structures.',
        },
        language,
      );
      if (isFailure(ms)) return ms;
      const locis = chainLoci(ms, traceOnly);
      if (isFailure(locis)) return locis;
      const moves = await ctx.plugin.runTask(
        Task.create('TM-align superposition', async (rctx) => {
          const out: Move[] = [];
          for (let i = 1; i < ms.length; i++) {
            if (rctx.shouldUpdate) await rctx.update(`Superposing pair ${i} of ${ms.length - 1}...`);
            const r = tmAlign(locis[0], locis[i]);
            out.push({
              member: ms[i],
              matrix: r.bTransform,
              rmsd: r.rmsd,
              paired: r.alignedLength,
              extra: { tmScore: r.tmScoreA, tmScoreMobile: r.tmScoreB },
            });
          }
          return out;
        }),
      );
      const bad = moves.find((m) => m.paired === 0 || Number.isNaN(m.rmsd));
      if (bad) {
        return ApiResult.fail(
          'failed',
          `TM-align found nothing in common between ${ms[0].label} and ${bad.member.label}.`,
          'Check that the chains are of the same kind of molecule.',
        );
      }
      const results = await apply(ctx, 'tm-align', ms[0], moves);
      return finish(ctx, 'tm-align', ms[0], results, 'by TM-align');
    },
  });

  readonly ligands = this.operation({
    id: 'superpose.ligands',
    title: 'Superpose by ligands',
    params: {
      language: LanguageParam(this.ctx.languages),
      targets: TargetsParam(
        'The ligand (one residue) to superpose',
        ' Matched by atom names when the compound is the same, else by the maximum common connected subgraph.',
      ),
    },
    when: needsTwo,
    async run({ targets, language }, ctx) {
      const ms = await members(
        ctx,
        targets,
        {
          what: 'a ligand',
          hint: 'Select one ligand residue in each of two or more structures.',
        },
        language,
      );
      if (isFailure(ms)) return ms;
      for (const m of ms) {
        if (!isSingleLigandLoci(m.loci)) {
          return ApiResult.fail(
            'invalid',
            `The target in ${m.label} is not one ligand residue (${m.part?.selected.line}).`,
            'Select exactly one ligand residue per structure.',
          );
        }
      }
      const warnings: string[] = [];
      const moves = await ctx.plugin.runTask(
        Task.create('Ligand superposition', async (rctx) => {
          const out: Move[] = [];
          for (let i = 1; i < ms.length; i++) {
            if (rctx.shouldUpdate) await rctx.update(`Superposing ligand ${i} of ${ms.length - 1}...`);
            const r = superposeLigandsByMccs(ms[0].loci, ms[i].loci, DefaultLigandMccsOptions);
            if (!r) {
              warnings.push(`The ligands of ${ms[i].label} and ${ms[0].label} have too little in common to superpose.`);
              continue;
            }
            if (r.truncated) {
              warnings.push(`The search for ${ms[i].label} hit its time budget; the match may not be the best.`);
            }
            out.push({
              member: ms[i],
              matrix: r.bTransform,
              rmsd: r.rmsd,
              paired: r.atomCount,
              extra: { matched: r.method },
            });
          }
          return out;
        }),
      );
      if (!moves.length) {
        return ApiResult.fail(
          'failed',
          warnings.join(' '),
          'Ligands superpose when they share atoms or a common scaffold.',
        );
      }
      const results = await apply(ctx, 'ligands', ms[0], moves);
      return finish(ctx, 'ligands', ms[0], results, 'by ligands', warnings);
    },
  });

  readonly atoms = this.operation({
    id: 'superpose.atoms',
    title: 'Superpose by atoms',
    params: {
      language: LanguageParam(this.ctx.languages),
      targets: TargetsParam(
        'The atoms to superpose',
        ` The same number in each, paired in the order the structure holds them; at least ${MinAtoms}.`,
      ),
    },
    when: needsTwo,
    async run({ targets, language }, ctx) {
      const ms = await members(
        ctx,
        targets,
        {
          what: 'atoms',
          hint: `Select the same atoms, in the same order, in each of two or more structures (at least ${MinAtoms}).`,
        },
        language,
      );
      if (isFailure(ms)) return ms;
      const counts = ms.map((m) => StructureElement.Loci.size(m.loci));
      if (counts.some((n) => n !== counts[0])) {
        return ApiResult.fail(
          'invalid',
          `The targets have different numbers of atoms: ${ms.map((m, i) => `${m.label} ${counts[i]}`).join(', ')}.`,
          'Atoms are paired one to one: select the same number in each structure.',
        );
      }
      if (counts[0] < MinAtoms) {
        return ApiResult.fail(
          'invalid',
          `${counts[0]} atom${counts[0] === 1 ? '' : 's'} per structure: too few to superpose (${MinAtoms} are needed).`,
          'Select more atoms.',
        );
      }
      const transforms = superpose(ms.map((m) => m.loci));
      const moves: Move[] = [];
      for (let i = 1; i < ms.length; i++) {
        const { bTransform, rmsd, nAlignedElements } = transforms[i - 1];
        if (Number.isNaN(rmsd)) return noRmsd(ms[0].label);
        moves.push({ member: ms[i], matrix: bTransform, rmsd, paired: nAlignedElements, extra: {} });
      }
      const results = await apply(ctx, 'atoms', ms[0], moves);
      return finish(ctx, 'atoms', ms[0], results, 'by atoms');
    },
  });

  readonly db = this.operation({
    id: 'superpose.db',
    title: 'Superpose by UniProt (SIFTS)',
    params: {
      structures: PD.Value<readonly string[]>([], {
        description:
          'The structures to superpose by label or state ref, the first is the reference; empty: every loaded one.',
      }),
      traceOnly: SuperpositionOptions.traceOnly,
    },
    when: needsTwo,
    async run({ structures, traceOnly }, ctx) {
      const all = loadedStructures(ctx);
      const hasMapping = (s: HierarchyStructure) =>
        !!s.cell.obj?.data.models.some((m) => SIFTSMapping.Provider.isApplicable(m));
      let picked: HierarchyStructure[] = [...all];
      // Every loaded one means every one the mapping can place: the others are left out, and said so.
      const skipped = structures.length ? [] : all.filter((s) => s.cell.obj && !hasMapping(s));
      if (!structures.length) picked = all.filter((s) => !skipped.includes(s));
      if (structures.length) {
        picked = [];
        for (const name of structures) {
          const s = pickNamed(all, name, {
            what: 'structure',
            ids: (x) => [x.cell.transform.ref],
            names: (x) => [labelOf(ctx, x), x.cell.obj?.label],
            describe: (x) => `${labelOf(ctx, x)} (${x.cell.transform.ref})`,
            none: { message: 'Nothing is loaded.', hint: 'Load structures first.' },
          });
          if (isFailure(s)) return s;
          if (picked.includes(s)) return ApiResult.fail('invalid', `${labelOf(ctx, s)} is named twice.`);
          picked.push(s);
        }
      }
      if (picked.length < 2) {
        if (skipped.length) return noMapping(ctx, skipped);
        return ApiResult.fail(
          'invalid',
          'Superposing needs two or more structures.',
          'Name two or more, or load more.',
        );
      }
      const datas = picked.map((s) => s.cell.obj?.data);
      if (datas.some((d) => !d)) return ApiResult.fail('failed', 'A structure is not ready yet.');
      const lacking = picked.filter((s) => !hasMapping(s));
      if (lacking.length) return noMapping(ctx, lacking);
      const label = (i: number) => labelOf(ctx, picked[i]);
      const { entries, failedPairs, zeroOverlapPairs } = alignAndSuperposeWithSIFTSMapping(datas as Structure[], {
        traceOnly,
      });
      const warnings = [
        ...skipped.map((s) => `${labelOf(ctx, s)} was left out: it has no SIFTS UniProt mapping.`),
        ...zeroOverlapPairs.map(([i, j]) => `${label(i)} and ${label(j)} share no UniProt residues.`),
        ...failedPairs.map(([i, j]) => `The superposition of ${label(j)} onto ${label(i)} failed.`),
      ];
      if (!entries.length) {
        return ApiResult.fail(
          'failed',
          warnings.join(' ') || 'Nothing could be superposed.',
          'The structures share no UniProt entry.',
        );
      }
      const reference = picked[0];
      const members0: Member = {
        structure: reference,
        data: reference.cell.obj!.data,
        label: label(0),
        loci: StructureElement.Loci.all(reference.cell.obj!.data),
      };
      const moves: Move[] = entries.map((e) => ({
        member: {
          structure: picked[e.other],
          data: picked[e.other].cell.obj!.data,
          label: label(e.other),
          loci: StructureElement.Loci.all(picked[e.other].cell.obj!.data),
        },
        matrix: e.transform.bTransform,
        rmsd: e.transform.rmsd,
        paired: e.transform.nAlignedElements,
        extra: {},
      }));
      const results = await apply(ctx, 'db', members0, moves);
      return finish(ctx, 'db', members0, results, 'by UniProt mapping', warnings);
    },
  });
}
