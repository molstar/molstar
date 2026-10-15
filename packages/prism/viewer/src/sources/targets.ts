/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Structure, StructureElement } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import type { ApiErrorCode, ApiFailure } from '../api/area.js';
import type { ResolvedPart, ResolvedTarget, ResolveOptions, TargetResolver } from '../api/targets.js';
import { findStructure, structureNames } from './naming.js';
import { isLociTarget } from '../api/target-param.js';
import {
  evaluateTarget,
  summarize,
  formatTarget,
  type PersonKeyword,
  type Target,
  toMolqlJson,
} from '../targets/index.js';

// The target resolver over the plugin (targets.md §1–2, plan 2.1): a target, in any of its forms, compiled to MolQL
// (the `targets/` port) and evaluated on the structures it applies to. Named structure (a state ref or a label), else
// the focused one, else the one with a selection, else every loaded one; the result says which it was. `selection` and
// `focus` are the person's: Mol*'s selection and focus on that structure.

export function createTargetResolver(plugin: PluginContext): TargetResolver {
  const structures = () => plugin.managers.structure.hierarchy.current.structures;

  const byName = (name: string): StructureRef | undefined => findStructure(structures(), name);
  const nameOf = (s: StructureRef) => structureNames(structures()).get(s.cell.transform.ref) ?? s.cell.transform.ref;

  const focused = (): StructureRef | undefined => {
    const loci = plugin.managers.structure.focus.current?.loci;
    if (!loci) return undefined;
    const root = loci.structure.root;
    return structures().find((s) => s.cell.obj?.data === loci.structure || s.cell.obj?.data.root === root);
  };

  const withSelection = (): StructureRef[] =>
    structures().filter((s) => {
      const data = s.cell.obj?.data;
      if (!data) return false;
      const loci = plugin.managers.structure.selection.getLoci(data);
      return StructureElement.Loci.is(loci) && !StructureElement.Loci.isEmpty(loci);
    });

  const personLoci = (which: PersonKeyword, structure: Structure): StructureElement.Loci | undefined => {
    if (which === 'selection') {
      const loci = plugin.managers.structure.selection.getLoci(structure);
      return StructureElement.Loci.is(loci) && !StructureElement.Loci.isEmpty(loci) ? loci : undefined;
    }
    const focus = plugin.managers.structure.focus.current?.loci;
    if (!focus || focus.structure.root !== structure.root) return undefined;
    const remapped = StructureElement.Loci.remap(focus, structure);
    return StructureElement.Loci.is(remapped) ? remapped : undefined;
  };

  return {
    async resolve(target: unknown, options: ResolveOptions = {}): Promise<ResolvedTarget | ApiFailure> {
      if (target === undefined || target === null || target === '') {
        return fail('invalid', 'Which part?', 'Give a target: a keyword (ligand), an address (A:120-130) or a query.');
      }
      const all = structures();
      if (!all.length) return fail('not-found', 'Nothing is loaded.', 'Load a structure first.');
      // Loci from code in the page (the sequence link): on the structure they belong to, never stored.
      if (isLociTarget(target)) {
        const loci = target.loci;
        const s = all.find((x) => x.cell.obj?.data.root === loci.structure.root);
        if (!s) return fail('not-found', 'Those loci are on no loaded structure.');
        const selected = summarize(loci);
        return {
          ok: true,
          parts: [{ structure: { ref: s.cell.transform.ref, label: nameOf(s) }, loci, selected }],
          scope: 'named',
          summary: selected.line,
          empty: selected.atoms === 0,
          stored: undefined,
        };
      }

      let scope: StructureRef[];
      let scoping: ResolvedTarget['scope'];
      if (options.structure) {
        const s = byName(options.structure);
        if (!s) {
          return fail(
            'not-found',
            `There is no structure '${options.structure}'.`,
            `Loaded: ${all.map(nameOf).join(', ')}.`,
          );
        }
        scope = [s];
        scoping = 'named';
      } else {
        const f = focused();
        const sel = withSelection();
        if (f) {
          scope = [f];
          scoping = 'focused';
        } else if (target === 'selection' && sel.length) {
          scope = sel;
          scoping = 'selected';
        } else if (sel.length === 1) {
          scope = sel;
          scoping = 'selected';
        } else {
          scope = all;
          scoping = all.length === 1 ? 'only' : 'all';
        }
      }

      const parts: ResolvedPart[] = [];
      let failure: ApiFailure | undefined;
      let stored: unknown;
      for (const s of scope) {
        const structure = s.cell.obj?.data;
        if (!structure) continue;
        const r = evaluateTarget(structure, target as Target, {
          ...(options.language ? { language: options.language } : {}),
          resolvePerson: (which) => {
            const loci = personLoci(which, structure);
            return loci ? StructureElement.Bundle.toExpression(StructureElement.Bundle.fromLoci(loci)) : undefined;
          },
        });
        if (!r.ok) {
          failure ??= fail('invalid', r.message, r.hint);
          continue;
        }
        stored ??= toMolqlJson(r.target.molql);
        parts.push({
          structure: { ref: s.cell.transform.ref, label: nameOf(s) },
          loci: r.loci,
          selected: r.selected,
        });
      }
      if (!parts.length) return failure ?? fail('not-found', 'Nothing to resolve the target on.');
      const nonEmpty = parts.filter((p) => p.selected.atoms > 0);
      const summary =
        nonEmpty.length === 0
          ? `${formatTarget(target as Target)} selects nothing${parts.length > 1 ? ` in ${parts.length} structures` : ` in ${parts[0].structure.label}`}.`
          : nonEmpty
              .map(
                (p) => `${p.selected.line}${parts.length > 1 || scoping !== 'only' ? ` in ${p.structure.label}` : ''}`,
              )
              .join('; ');
      return {
        ok: true,
        parts,
        scope: scoping,
        summary,
        empty: nonEmpty.length === 0,
        stored,
      };
    },
  };
}

function fail(code: ApiErrorCode, message: string, hint?: string): ApiFailure {
  return hint ? { ok: false, code, message, hint } : { ok: false, code, message };
}
