/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Unit } from '../../unit.js';
import { BondProvider } from './provider.js';
import { EntityType } from '../../../model/types.js';
import { stringToWords } from '@molstar/core/util/string';

/**
 * The generic bond providers (bond-providers plan §2.4, §5): `'none'`, the built-in bonds without the distance-based
 * ones, and `'per-unit'`, other providers chosen per unit by rules.
 */
export { NoComputedBondsProvider, PerUnitBondProvider };

const NoComputedBondsName = 'none';

const NoComputedBondsInstance: BondProvider.Instance = {
  provider: undefined as any as BondProvider,
  cacheKey: NoComputedBondsName,
  intra: (_, ctx) => ctx.builtIn({ noDistance: true }),
  inter: (_, ctx) => ctx.builtIn({ noDistance: true }),
};

/**
 * The bonds the format states (bond tables, CONECT, struct_conn, the component dictionary) and none computed from
 * distances: never no bonds at all. A structure without them (a PDB entry without `chem_comp_bond`) has no peptide
 * bonds within its chains.
 */
const NoComputedBondsProvider: BondProvider<{}> = {
  name: NoComputedBondsName,
  label: 'None (format bonds only)',
  description: 'Only the bonds the file states; none computed from distances',
  getParams: () => ({}),
  create: () => NoComputedBondsInstance,
};
(NoComputedBondsInstance as { provider: BondProvider }).provider = NoComputedBondsProvider;

const PerUnitName = 'per-unit';

/** The mmCIF entity types (`EntityType`, the dictionary's), as options. */
const EntityTypes = Object.keys(EntityType)
  .filter((k) => Number.isNaN(Number(k)) && k !== 'unknown')
  .map((k) => [k, stringToWords(k)] as [string, string]);

/** What a rule matches: properties the same for every element of a chain and every operator (plan §5). */
const MatchParams = PD.MappedStatic('entity-type', {
  'entity-type': PD.MultiSelect(['non-polymer'], EntityTypes, { label: 'Entity types' }),
  'label-asym-id': PD.Text('', { label: 'Chains', description: 'label_asym_id, comma-separated' }),
  'entity-id': PD.Text('', { label: 'Entities', description: 'Entity ids, comma-separated' }),
});
type Match = PD.Values<{ match: typeof MatchParams }>['match'];

function getPerUnitParams(model: Parameters<BondProvider['getParams']>[0], ctx?: BondProvider.Context) {
  const exclude = [PerUnitName];
  // without a registry only the built-in bonds can be chosen
  const choice = () =>
    ctx?.choiceParam(model, { exclude }) ??
    PD.Mapped<any>(BondProvider.DefaultName, [[BondProvider.DefaultName, BondProvider.Default.label]], () =>
      PD.Group({}),
    );
  return {
    rules: PD.ObjectList(
      { match: MatchParams, provider: choice() },
      (r) =>
        `${r.match.name}: ${Array.isArray(r.match.params) ? r.match.params.join(', ') : r.match.params} → ${r.provider.name}`,
      { defaultValue: [] as { match: Match; provider: BondProvider.Choice }[] },
    ),
    inter: PD.Group({ provider: choice() }, { label: 'Between units', isFlat: true }),
  };
}
type PerUnitParams = ReturnType<typeof getPerUnitParams>;

const list = (text: string) =>
  new Set(
    text
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s.length > 0),
  );

/**
 * Whether a unit matches, read on its first element: a function of the model and the elements alone, as the intra-unit
 * bond cache requires (equal elements, equal bonds). A symmetry copy and a frame have the elements and so the match of
 * their unit; a child of a unit with several chains (merged waters) matches by its own first chain.
 */
function matcher(match: Match): (unit: Unit.Atomic) => boolean {
  const entityOf = (u: Unit.Atomic) => u.model.atomicHierarchy.index.getEntityFromChain(u.chainIndex[u.elements[0]]);
  switch (match.name) {
    case 'entity-type': {
      const types = new Set(match.params);
      return (u) => types.has(u.model.entities.data.type.value(entityOf(u)));
    }
    case 'label-asym-id': {
      const ids = list(match.params);
      return (u) => ids.has(u.model.atomicHierarchy.chains.label_asym_id.value(u.chainIndex[u.elements[0]]));
    }
    case 'entity-id': {
      const ids = list(match.params);
      return (u) => ids.has(u.model.entities.data.id.value(entityOf(u)));
    }
  }
}

/**
 * Other providers chosen per unit: the first rule matching a unit decides its intra-unit bonds, units no rule
 * matches get the built-in ones, and one provider decides the bonds between units. Serializable as plain params, so
 * per-unit choice needs nothing of the model layer beyond units with their own instances (plan §5).
 */
const PerUnitBondProvider: BondProvider<PerUnitParams> = {
  name: PerUnitName,
  label: 'By unit',
  description: 'Choose the provider of each unit by entity type, chain or entity',
  getParams: getPerUnitParams,
  create(props, ctx) {
    const resolve = (choice: BondProvider.Choice | undefined) => ctx?.resolve(choice) ?? BondProvider.DefaultInstance;
    const rules = props.rules.map((r) => ({ matches: matcher(r.match), instance: resolve(r.provider) }));
    const inter = resolve(props.inter.provider);
    const instances = [...rules.map((r) => r.instance), inter];
    const keys = instances.map((i) => i.cacheKey);
    const instanceFor = (unit: Unit.Atomic) => {
      for (const r of rules) if (r.matches(unit)) return r.instance;
      return BondProvider.DefaultInstance;
    };
    return {
      provider: PerUnitBondProvider,
      // equal rules and equal children give equal bonds
      cacheKey: keys.some((k) => k === undefined)
        ? undefined
        : `${PerUnitName}:${JSON.stringify(props.rules.map((r) => r.match))}:${keys.join('|')}`,
      // one setting for every unit: a child computes its root's bonds when any rule's provider asks for that
      childBonds: instances.some((i) => i.childBonds === 'root') ? 'root' : 'self',
      intra(unit, c) {
        const instance = instanceFor(unit);
        return instance.intra ? instance.intra(unit, c) : c.builtIn();
      },
      inter: inter.inter ? (structure, c) => inter.inter!(structure, c) : undefined,
    };
  },
};
