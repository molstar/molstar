/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// https://pymol.org/dokuwiki/doku.php?id=selection:alpha
// https://github.com/evonove/pymol/blob/master/pymol/layer3/Selector.cpp

import * as P from 'parsimmon'

import * as h from '../helper'
import { AtomGroupArgs } from '../types'

import properties from './properties'
import operators from './operators'
import keywords from './keywords'

import Transpiler from '../transpiler'
import B from '../../molql/builder'

const propertiesDict = h.getPropertyRules(properties)

const slash = P.string('/')

function orNull(rule: P.Parser<any>) {
  return rule.or(P.of(null))
}

function atomSelectionQuery(x: any) {
  const tests: AtomGroupArgs = {}
  const props: {[k: string]: any[]} = {}

  for (let k in x) {
    const ps = properties[k]
    if (!ps) {
      throw new Error(`property '${k}' not supported, value '${x[k]}'`)
    }
    if (x[k] === null) continue
    if (!props[ps.level]) props[ps.level] = []
    props[ps.level].push(x[k])
  }

  for (let p in props) {
    tests[p] = h.andExpr(props[p])
  }

  return B.struct.generator.atomGroups(tests)
}

const lang = P.createLanguage({
  Parens: function (r) {
    return P.alt(
      r.Parens,
      r.Operator,
      r.Expression
    ).wrap(P.string('('), P.string(')'))
  },

  Expression: function(r) {
    return P.alt(
      r.AtomSelectionMacro.map(atomSelectionQuery),
      r.NamedAtomProperties,
      r.Pepseq,
      r.Rep,
      r.Keywords,
      r.Object
    )
  },

  AtomSelectionMacro: function(r) {
    return P.alt(
      slash.then(P.alt(
        P.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name)
        ).map(x => {return {object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4]}}),
        P.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi)
        ).map(x => {return {object: x[0], segi: x[1], chain: x[2], resi: x[3]}}),
        P.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain)
        ).map(x => {return {object: x[0], segi: x[1], chain: x[2]}}),
        P.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi)
        ).map(x => {return {object: x[0], segi: x[1]}}),
        P.seq(
          orNull(r.ObjectProperty)
        ).map(x => {return {object: x[0]}}),
      )),
      P.alt(
        P.seq(
          orNull(r.ObjectProperty).skip(slash),
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name)
        ).map(x => {return {object: x[0], segi: x[1], chain: x[2], resi: x[3], name: x[4]}}),
        P.seq(
          orNull(propertiesDict.segi).skip(slash),
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name)
        ).map(x => {return {segi: x[0], chain: x[1], resi: x[2], name: x[3]}}),
        P.seq(
          orNull(propertiesDict.chain).skip(slash),
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name)
        ).map(x => {return {chain: x[0], resi: x[1], name: x[2]}}),
        P.seq(
          orNull(propertiesDict.resi).skip(slash),
          orNull(propertiesDict.name)
        ).map(x => {return {resi: x[0], name: x[1]}}),
      )
    )
  },

  NamedAtomProperties: function() {
    return P.alt(...h.getNamedPropertyRules(properties))
  },

  Keywords: () => P.alt(...h.getKeywordRules(keywords)),

  ObjectProperty: () => {
    const w = h.getReservedWords(properties, keywords, operators)
      .sort(h.strLenSortFn).map(h.escapeRegExp).join('|')
    return P.regex(new RegExp(`(?!(${w}))[A-Z0-9_]+`, 'i'))
  },

  Object: (r) => {
    return r.ObjectProperty.notFollowedBy(slash)
      .map(x => { throw new Error(`property 'object' not supported, value '${x}'`) })
  },

  // Selects peptide sequence matching upper-case one-letter
  // sequence SEQ (see also FindSeq).
  // PEPSEQ seq
  Pepseq: () => {
    return P.regex(/(PEPSEQ|ps\.)\s+([a-z]+)/i, 2)
      .map(h.makeError(`operator 'pepseq' not supported`))
  },

  // Selects atoms which show representation rep.
  // REP rep
  Rep: () => {
    return P.regex(/REP\s+(lines|spheres|mesh|ribbon|cartoon|sticks|dots|surface|labels|extent|nonbonded|nb_spheres|slice|extent|slice|dashes|angles|dihedrals|cgo|cell|callback|everything)/i, 1)
      .map(h.makeError(`operator 'rep' not supported`))
  },

  Operator: function(r) {
    return h.combineOperators(operators, P.alt(r.Parens, r.Expression, r.Operator))
  },

  Query: function(r) {
    return P.alt(
      r.Operator,
      r.Parens,
      r.Expression
    ).trim(P.optWhitespace)
  }
})

const transpiler: Transpiler = str => lang.Query.tryParse(str)
export default transpiler
