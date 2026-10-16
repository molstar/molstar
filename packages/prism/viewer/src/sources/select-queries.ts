/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import type { Loci } from '@molstar/model/model/loci';
import { type Structure, StructureElement, StructureProperties, Unit } from '@molstar/model/model/structure';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { PluginContext } from '@molstar/plugin/context';
import {
  ElementSymbolQuery,
  EntityDescriptionQuery,
  ResidueQuery,
} from '@molstar/plugin/state/queries/structure/dynamic';
import type { StructureSelectionQuery } from '@molstar/plugin/state/queries/structure/query';
import type { CommandGroup, IconName, PreviewSink } from '@molstar/prism-ui/core';
import type { Expression } from '@molstar/query-language/language/expression';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';

// What the Select group's commands are made of (commands.md §5; the commands: `commands/select.ts`): the plugin's
// selection-query registry and what the loaded structures hold, as queries; and the preview sink that lights a
// highlighted row's atoms. The plugin is read here (a source, state.md §4), not in the commands.

/** A highlight to draw while a row is browsed (`PreviewSink`): the loci, resolved when it is drawn. */
export interface LociPreview {
  readonly kind: 'loci';
  resolve(): Promise<readonly Loci[]>;
}

const isLociPreview = (p: unknown): p is LociPreview =>
  !!p &&
  typeof p === 'object' &&
  (p as LociPreview).kind === 'loci' &&
  typeof (p as LociPreview).resolve === 'function';

/**
 * The viewer's preview sink (commands.md §4–5): a `LociPreview` lights its loci with the hover highlight; a newer preview
 * or a clear drops an older one still resolving. Other previews draw nothing.
 */
export function lociPreviewSink(plugin: PluginContext): PreviewSink {
  let token = 0;
  let shown = false;
  return {
    set(preview) {
      const mine = ++token;
      if (!isLociPreview(preview)) return;
      void preview.resolve().then(
        (loci) => {
          if (mine !== token) return;
          const highlights = plugin.managers.interactivity.lociHighlights;
          loci.forEach((l, i) =>
            i === 0 ? highlights.highlightOnly({ loci: l }, false) : highlights.highlight({ loci: l }, false),
          );
          if (!loci.length) highlights.clearHighlights();
          shown = true;
        },
        () => {},
      );
    },
    clear() {
      ++token;
      if (!shown) return;
      shown = false;
      plugin.managers.interactivity.lociHighlights.clearHighlights();
    },
  };
}

/** One selection query as a command: its title, where it sits, its expression. */
export interface QueryEntry {
  readonly id: string;
  readonly title: string;
  readonly group: string;
  readonly icon: IconName;
  readonly expression: Expression;
  readonly ensure?: StructureSelectionQuery['ensureCustomProperties'];
  readonly alias?: readonly string[];
}

const slug = (text: string) =>
  text
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '') || 'x';

/** Ids made unique within a list (a second "x" is "x-2"). */
function uniqueIds<T extends { id: string }>(list: readonly T[]): T[] {
  const seen = new Map<string, number>();
  return list.map((e) => {
    const n = (seen.get(e.id) ?? 0) + 1;
    seen.set(e.id, n);
    return n === 1 ? e : { ...e, id: `${e.id}-${n}` };
  });
}

/**
 * The registry's queries that select by themselves (not the current selection's: those are `select.around`,
 * `select.expand`), by category; each category a subgroup of By query.
 */
export function registryQueries(plugin: PluginContext): { groups: CommandGroup[]; entries: QueryEntry[] } {
  const groups: CommandGroup[] = [];
  const entries: QueryEntry[] = [];
  const list = plugin.query.structure.registry.list.filter(
    (q) => !q.isHidden && !q.referencesCurrent && q.category && q.category !== 'Internal',
  );
  for (const q of list) {
    const group = `select.query.${slug(q.category)}`;
    if (!groups.some((g) => g.id === group))
      groups.push({ id: group, title: q.category, parent: 'select.query', order: groups.length });
    entries.push({
      id: `select.query.${slug(q.category)}.${slug(q.label)}`,
      title: q.label,
      group,
      icon: 'mouse-pointer-click',
      expression: q.expression,
      ...(q.ensureCustomProperties ? { ensure: q.ensureCustomProperties } : {}),
    });
  }
  return { groups, entries: uniqueIds(entries) };
}

/** Standard residues (amino acids, nucleotides, water) are not listed as residues of their own. */
const Standard = new Set([
  ...'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL SEC PYL UNK'.split(' '),
  ...'A C G I U DA DC DG DI DT DU N DN'.split(' '),
  ...'HOH DOD WAT H2O'.split(' '),
]);

/** At most this many of each kind (a ribosome's chains): the rest are found by `Select…` with text. */
const MaxEach = 120;

/** What the loaded structures hold, as queries: entities, chains, non-standard residues, elements. */
export function structureQueries(structures: readonly Structure[]): QueryEntry[] {
  const entities = new Map<string, string[]>();
  const chains = new Map<string, string>();
  const residues = new Map<string, string>();
  const elements = new Set<string>();
  const l = StructureElement.Location.create();
  for (const structure of structures) {
    l.structure = structure;
    for (const ug of structure.unitSymmetryGroups) {
      l.unit = ug.units[0];
      l.element = ug.elements[0];
      const type = StructureProperties.entity.type(l);
      if (type === 'polymer' || type === 'branched') {
        const description = StructureProperties.entity.pdbx_description(l);
        entities.set(description.join(', '), description);
      }
      // Author chain ids are atomic chains' (a coarse model's are its label chains: `Select…` reaches them).
      const chain = Unit.isAtomic(l.unit) ? StructureProperties.chain.auth_asym_id(l) : '';
      if (chain && !chains.has(chain)) chains.set(chain, StructureProperties.entity.pdbx_description(l).join(', '));
    }
    for (const r of structure.uniqueResidueNames) {
      if (Standard.has(r) || residues.has(r)) continue;
      const name = structure.models.map((m) => m.properties.chemicalComponentMap.get(r)?.name).find((n) => !!n);
      residues.set(r, name ?? r);
    }
    for (const e of structure.uniqueElementSymbols) elements.add(e);
  }
  const first = <T>(xs: Iterable<T>) => [...xs].slice(0, MaxEach);
  const out: QueryEntry[] = [];
  for (const [key, names] of first(entities)) {
    const q = EntityDescriptionQuery([names, key], '', 0);
    out.push({
      id: `select.entity.${slug(key)}`,
      title: key || 'Unnamed entity',
      group: 'select.structure.entities',
      icon: 'pick-entity',
      expression: q.expression,
    });
  }
  for (const [id, description] of first(chains)) {
    out.push({
      id: `select.chain.${id}`,
      title: `Chain ${id}`,
      group: 'select.structure.chains',
      icon: 'pick-chain',
      alias: description ? [description] : [],
      expression: MS.struct.generator.atomGroups({ 'chain-test': MS.core.rel.eq([MS.ammp('auth_asym_id'), id]) }),
    });
  }
  for (const [code, name] of first(residues)) {
    const q = ResidueQuery([[code], name], '', 0);
    out.push({
      id: `select.residue.${code}`,
      title: name === code ? code : `${code}: ${name}`,
      group: 'select.structure.residues',
      icon: 'ligand',
      alias: [code],
      expression: q.expression,
    });
  }
  for (const e of first(elements)) {
    const q = ElementSymbolQuery([[e], e], '', 0);
    out.push({
      id: `select.element.${e}`,
      title: q.label,
      group: 'select.structure.elements',
      icon: 'pick-atom',
      expression: q.expression,
    });
  }
  return uniqueIds(out);
}

/** The structures loaded now (their data). */
export function loadedStructureData(plugin: PluginContext): Structure[] {
  return plugin.managers.structure.hierarchy.current.structures
    .map((s) => s.cell.obj?.data)
    .filter((s): s is Structure => !!s);
}

/** The Select group's queries from the plugin: the registry's and the loaded structures' (read when asked). */
export interface SelectQuerySource {
  /** The registry's categories (By query's subgroups) and every query now. */
  list(): { groups: CommandGroup[]; entries: QueryEntry[] };
  /** Computes what a query needs (secondary structure for Helix) on the loaded structures. */
  ensure(entry: QueryEntry): Promise<void>;
}

export function selectQuerySource(plugin: PluginContext): SelectQuerySource {
  return {
    list() {
      const registry = registryQueries(plugin);
      const structures = loadedStructureData(plugin);
      return {
        groups: registry.groups,
        entries: [...registry.entries, ...(structures.length ? structureQueries(structures) : [])],
      };
    },
    async ensure(entry) {
      if (!entry.ensure) return;
      const structures = loadedStructureData(plugin);
      await plugin.runTask(
        Task.create('Selection query', async (runtime) => {
          const ctx: CustomProperty.Context = { runtime, assetManager: plugin.managers.asset };
          for (const s of structures) await entry.ensure!(ctx, s);
        }),
      );
    },
  };
}
