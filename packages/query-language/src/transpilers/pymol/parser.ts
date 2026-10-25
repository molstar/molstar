/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Panagiotis Tourlas <panagiot_tourlov@hotmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * Adapted from MolQL project
 */

// https://pymol.org/dokuwiki/doku.php?id=selection:alpha
// https://pymolwiki.org/index.php/Selection_Algebra
// https://github.com/evonove/pymol/blob/master/pymol/layer3/Selector.cpp

import * as P from '@molstar/core/util/monadic-parser';
import * as h from '../helper.js';
import { MolScriptBuilder } from '../../language/builder.js';
const B = MolScriptBuilder;
import { properties } from './properties.js';
import { operators } from './operators.js';
import { keywords } from './keywords.js';
import type { AtomGroupArgs } from '../types.js';
import type { TranspileOptions, Transpiler } from '../transpiler.js';
import type { Expression } from '../../language/expression.js';

const propertiesDict = h.getPropertyRules(properties);

const slash = P.MonadicParser.string('/');

/** The options of the call being parsed: the grammar is made once, so a call's names are read from here. */
let current: TranspileOptions | undefined;

/** An object or a named selection, through the caller's `resolveName`. */
function resolveName(name: string): Expression {
  const resolve = current?.resolveName;
  if (!resolve) throw new Error(`property 'object' not supported, value '${name}'`);
  const e = resolve(name);
  if (!e) throw new Error(`no object or selection named '${name}'`);
  return e;
}

function orNull(rule: P.MonadicParser<any>) {
  return rule.or(P.MonadicParser.of(null));
}

function atomSelectionQuery(x: any) {
  const tests: AtomGroupArgs = {};
  const props: { [k: string]: any[] } = {};

  // The macro's object field (`/1abc//A`) names an object: its atoms, intersected with the rest of the macro.
  const object = x.object ? resolveName(x.object) : undefined;
  let rest = false;
  for (const k in x) {
    if (k === 'object') continue;
    const ps = properties[k];
    if (!ps) {
      throw new Error(`property '${k}' not supported, value '${x[k]}'`);
    }
    if (x[k] === null) continue;
    rest = true;
    if (!props[ps.level]) props[ps.level] = [];
    props[ps.level].push(x[k]);
  }

  for (const p in props) {
    tests[p] = h.andExpr(props[p]);
  }

  const groups = B.struct.generator.atomGroups(tests);
  if (!object) return groups;
  return rest ? B.struct.modifier.intersectBy({ 0: groups, by: object }) : object;
}

/** The prefix operators (`not`, `byres`, `bound_to`, …): their words and maps, as `combineOperators` applies them. */
const Prefixes = operators
  .filter((o) => o.type === h.prefix)
  .map((o) => ({ rule: o.rule, map: o.isUnsupported ? h.makeError(`operator '${o.name}' not supported`) : o.map }));

const lang = P.MonadicParser.createLanguage({
  Parens: function (r: any) {
    // An operator chain holds a lone group or expression too (`((a) within 5 of b)`, `(chain A)`): one alternative, so
    // an unclosed group is read once per level, not once per alternative (three a level was exponential in the depth).
    return r.Operator.wrap(P.MonadicParser.string('('), P.MonadicParser.string(')'));
  },

  Expression: function (r: any) {
    return P.MonadicParser.alt(
      r.Keywords,
      r.AtomSelectionMacro.map(atomSelectionQuery),
      r.NamedAtomProperties,
      r.Pepseq,
      r.Rep,
      r.Object,
    );
  },

  AtomSelectionMacro: function (r: any) {
    return P.MonadicParser.alt(
      slash.then(
        P.MonadicParser.alt(
          P.MonadicParser.seq(
            orNull(r.ObjectProperty).skip(slash),
            orNull(propertiesDict.segi).skip(slash),
            orNull(propertiesDict.chain).skip(slash),
            orNull(propertiesDict.resi).skip(slash),
            orNull(propertiesDict.name),
          ).map((x) => {
            return { object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4] };
          }),
          P.MonadicParser.seq(
            orNull(r.ObjectProperty).skip(slash),
            orNull(propertiesDict.segi).skip(slash),
            orNull(propertiesDict.chain).skip(slash),
            orNull(propertiesDict.resi),
          ).map((x) => {
            return { object: x[0], segi: x[1], chain: x[2], resi: x[3] };
          }),
          P.MonadicParser.seq(
            orNull(r.ObjectProperty).skip(slash),
            orNull(propertiesDict.segi).skip(slash),
            orNull(propertiesDict.chain),
          ).map((x) => {
            return { object: x[0], segi: x[1], chain: x[2] };
          }),
          P.MonadicParser.seq(orNull(r.ObjectProperty).skip(slash), orNull(propertiesDict.segi)).map((x) => {
            return { object: x[0], segi: x[1] };
          }),
          P.MonadicParser.seq(orNull(r.ObjectProperty)).map((x) => {
            return { object: x[0] };
          }),
        ),
      ),
      P.MonadicParser.alt(
        P.MonadicParser.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name),
        ).map((x) => {
          return { object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4] };
        }),
        P.MonadicParser.seq(
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name),
        ).map((x) => {
          return { segi: x[0], chain: x[1], resi: x[2], name: x[3] };
        }),
        P.MonadicParser.seq(
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name),
        ).map((x) => {
          return { chain: x[0], resi: x[1], name: x[2] };
        }),
        P.MonadicParser.seq(orNull(propertiesDict.resi).skip(slash), orNull(propertiesDict.name)).map((x) => {
          return { resi: x[0], name: x[1] };
        }),
      ),
    );
  },

  NamedAtomProperties: function () {
    return P.MonadicParser.alt(...h.getNamedPropertyRules(properties));
  },

  Keywords: () => P.MonadicParser.alt(...h.getKeywordRules(keywords)),

  ObjectProperty: () => {
    const w = h.getReservedWords(properties, keywords, operators).sort(h.strLenSortFn).map(h.escapeRegExp).join('|');
    // Not a reserved word (a keyword, property or operator) as a whole word: a name may start with one (`b2`, `both`).
    return P.MonadicParser.regexp(new RegExp(`(?!(${w})(?![A-Z0-9_]))[A-Z0-9_]+`, 'i'));
  },
  Object: (r: any) => {
    return r.ObjectProperty.notFollowedBy(slash).map((x: any) => resolveName(x));
  },

  // Selects peptide sequence matching upper-case one-letter
  // sequence SEQ (see also FindSeq).
  // PEPSEQ seq
  Pepseq: () => {
    return P.MonadicParser.regexp(/(PEPSEQ|ps\.)\s+([a-z]+)/i, 2).map(h.makeError(`operator 'pepseq' not supported`));
  },

  // Selects atoms which show representation rep.
  // REP rep
  Rep: () => {
    return P.MonadicParser.regexp(
      /REP\s+(lines|spheres|mesh|ribbon|cartoon|sticks|dots|surface|labels|extent|nonbonded|nb_spheres|slice|extent|slice|dashes|angles|dihedrals|cgo|cell|callback|everything)/i,
      1,
    ).map(h.makeError(`operator 'rep' not supported`));
  },

  Operator: function (r: any) {
    // An operand may itself be a prefix operator's (`byres bound_to resn CYS`): the operator read here, then a whole
    // `Operator` after it. Each step reads the operator's word first, so the recursion always moves on; the rule used
    // to list `Operator` itself, which re-entered where nothing read and recursed without end (a stack overflow on
    // `chain`, `b` or an empty text).
    const prefixed = P.MonadicParser.alt(
      ...Prefixes.map((p) => P.MonadicParser.seq(p.rule, r.Operator).map((x) => p.map(...x))),
    );
    return h.combineOperators(operators, P.MonadicParser.alt(r.Parens, r.Expression, prefixed));
  },

  Query: function (r: any) {
    return P.MonadicParser.alt(r.Operator, r.Parens, r.Expression).trim(P.MonadicParser.optWhitespace);
  },
});

export const transpiler: Transpiler = (str, options) => {
  const previous = current;
  current = options;
  try {
    return lang.Query.tryParse(str);
  } finally {
    current = previous;
  }
};
