/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The extension's one registry entry (seqspec:format#extensions; seqspec:libraries#interactions): the
// `interactions:fingerprint` provider and layer type, the `interactions:kind` color and glyph themes, and its issue
// code. Pure data.

import {
  list,
  literal,
  nullable,
  OptionalField,
  type ProviderEntry,
  RequiredField,
  ref,
  SimpleParamsSchema,
  type SvsRegistryEntry,
  type ThemeRef,
} from '@molstar/seq-core/format';
import {
  type AnnotationInput,
  type CoordinateSystem,
  Issue,
  KeyKind,
  type KeyTarget,
  Mappings,
} from '@molstar/seq-core/model';
import { interactionColors } from './builder.js';
import { type Fingerprint, Fingerprints } from './fingerprint.js';
import { KindMask, type KindName, Kinds } from './vocabulary.js';

/** The runtime property (seqspec:model#custom-data) the provider attaches to its matrix annotation. */
export const FingerprintProp = 'interactions:fingerprint';
export const fingerprintOf = (a: { custom: { props: Map<string, unknown> } }) =>
  a.custom.props.get(FingerprintProp) as Fingerprint | undefined;

interface Options {
  pairs: string[];
  ligands: string[];
  system: string;
  kinds: KindName[] | null;
}

export const FingerprintProvider: ProviderEntry = {
  name: 'interactions:fingerprint',
  geometry: 'matrix',
  options: SimpleParamsSchema({
    pairs: RequiredField(list(ref('annotation')), 'pair annotations'),
    ligands: RequiredField(list(ref('unit')), 'ligands or poses, a row each'),
    system: RequiredField(
      ref('unit', 'sequence', 'alignment'),
      "the columns' coordinate system: a unit, a sequence or an alignment's columns",
    ),
    kinds: OptionalField(nullable(list(literal(...Kinds))), null, 'only these kinds; null: all'),
  }),
  /**
   * A matrix over the frame: a row per contacted position, the ligands as states, the number of kinds as the value;
   * a ligand with no contact at a position is a missing value (null), drawn empty, as KLIFS shows it.
   */
  async compute(ctx, options) {
    const o = options as Options;
    const empty: AnnotationInput = { key: '', geometry: 'matrix', states: ['none'], columns: {} };
    const set = ctx.set;
    if (!set) return empty;
    const t: KeyTarget = { kind: -1, id: -1, row: -1 };
    const idOf = (key: string, kind: number) => (set.byName(key, t) && t.kind === kind ? t.id : -1);
    const bad = (path: string, what: string) =>
      ctx.report(Issue.error('model.provider.input', `/params/options/${path}`, `interactions:fingerprint: ${what}.`));
    const pairs = o.pairs.flatMap((key, j) => {
      const a = set.annotations[idOf(key, KeyKind.annotation)];
      if (!a?.pairs) bad(`pairs/${j}`, `"${key}" is not a pairs annotation`);
      return a?.pairs ? [a] : [];
    });
    const states = o.ligands.filter((key, j) => {
      const known = idOf(key, KeyKind.unit) >= 0;
      if (!known) bad(`ligands/${j}`, `no unit "${key}"`);
      return known;
    });
    const system: CoordinateSystem = { kind: -1, id: -1, row: -1 };
    if (!Mappings.systemOf(set, o.system, system) || system.row >= 0) {
      bad('system', `no unit, sequence or alignment "${o.system}"`);
      return empty;
    }
    if (!states.length) return empty;
    const issues: Issue[] = [];
    const kinds = o.kinds ? KindMask.of(...o.kinds) : KindMask.all;
    const ligands = Int32Array.from(states, (key) => idOf(key, KeyKind.unit));
    const fp = Fingerprints.build(set, pairs, ligands, system, { kinds, issues });
    for (const i of issues) ctx.report({ ...i, path: '/params/options/system' });
    const C = fp.columns.length;
    const value = Array.from(fp.columns, (_, c) =>
      Array.from(ligands, (_, r) => (fp.cells[r * C + c] ? KindMask.count(fp.cells[r * C + c]) : null)),
    );
    return {
      key: '',
      geometry: 'matrix',
      semanticType: 'interactions:fingerprint',
      label: 'Interaction fingerprint',
      target: o.system,
      states,
      columns: { index: Array.from(fp.columns), value },
      props: { [FingerprintProp]: fp },
    };
  },
};

/** The fingerprint's `counts` channel: the kind counts of its matrix fallback, `numericDefaultColor`'s Blues over [0, hi]. */
const countsColor: ThemeRef = {
  theme: 'field_palette',
  params: { palette: { kind: 'continuous', colors: 'Blues', value_domain: [0, null] } },
};

export const InteractionsRegistry: SvsRegistryEntry = {
  providers: [FingerprintProvider],
  layerTypes: [
    {
      name: 'interactions:fingerprint',
      description: 'the kinds of interaction per position and ligand, as KLIFS shows them',
      options: SimpleParamsSchema({}),
      accepts: ['matrix'],
      layouts: ['tracks'],
      axes: ['blocks', 'columns'],
      sizing: 'free',
      fallback: 'matrix',
      // `color`: the kinds. `counts`: the number of kinds, which a block without the fingerprint property draws (the
      // `matrix` fallback's cells); no glyph, the fingerprint draws none.
      channels: [
        { name: 'color', kind: 'color', unit: null, default: interactionColors() },
        { name: 'counts', kind: 'color', unit: null, default: countsColor },
      ],
    },
  ],
  // One entry for color and glyph; the references the builder writes carry their kernel fallbacks.
  themes: [
    {
      name: 'interactions:kind',
      description: "an interaction's kind (its category) through this vocabulary's colors and glyphs",
      params: { color: SimpleParamsSchema({}), glyph: SimpleParamsSchema({}) },
      fallback: null,
    },
  ],
  issueCodes: [{ code: 'interactions.fingerprint.unmapped', severity: 'warning' }],
};
