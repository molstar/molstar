/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as P from 'parsimmon'

import * as h from '../helper'
import { OperatorList } from '../types'

import properties from './properties'
import { sstrucMap, sstrucDict } from './properties'
import operators from './operators'
import keywords from './keywords'
import functions from './functions'

import Transpiler from '../transpiler'
import B from '../../molql/builder'

// const propertiesDict = h.getPropertyRules(properties)

// <, <=, = or ==, >=, >, and !=
// lt, le, eq, ge, gt, and ne, =~
const valueOperators: OperatorList = [
  {
    '@desc': 'multiplication, division',
    '@examples': [],
    name: 'mul-div',
    type: h.binaryLeft,
    rule: P.regex(/\s*(\*|\/)\s*/, 1),
    map: (op, e1, e2) => {
      switch (op) {
        case '*': return B.core.math.mult([e1, e2])
        case '/': return B.core.math.div([e1, e2])
        default: throw new Error(`value operator '${op}' not supported`);
      }
    }
  },
  {
    '@desc': 'addition, substraction',
    '@examples': [],
    name: 'add-sub',
    type: h.binaryLeft,
    rule: P.regex(/\s*(-|\+)\s*/, 1),
    map: (op, e1, e2) => {
      switch (op) {
        case '-': return B.core.math.sub([e1, e2])
        case '+': return B.core.math.add([e1, e2])
        default: throw new Error(`value operator '${op}' not supported`);
      }
    }
  },
  {
    '@desc': 'value comparisons',
    '@examples': [],
    name: 'comparison',
    type: h.binaryLeft,
    rule: P.alt(P.regex(/\s*(=~|==|>=|<=|=|!=|>|<)\s*/, 1), P.whitespace.result('=')),
    map: (op, e1, e2) => {
      // console.log(op, e1, e2)
      let expr
      if (e1.head === 'structure.atom-property.macromolecular.secondary-structure-flags') {
        expr = B.core.flags.hasAny([e1, sstrucMap(e2)])
      } else if (e2.head === 'structure.atom-property.macromolecular.secondary-structure-flags') {
        expr = B.core.flags.hasAny([e2, sstrucMap(e1)])
      } else if (e1.head === 'core.type.regex') {
        expr = B.core.str.match([ e1, B.core.type.str([e2]) ])
      } else if (e2.head === 'core.type.regex') {
        expr = B.core.str.match([ e2, B.core.type.str([e1]) ])
      } else if (op === '=~') {
        if (e1.head) {
          expr = B.core.str.match([
            B.core.type.regex([`^${e2}$`, 'i']),
            B.core.type.str([e1])
          ])
        } else {
          expr = B.core.str.match([
            B.core.type.regex([`^${e1}$`, 'i']),
            B.core.type.str([e2])
          ])
        }
      }
      if (!expr) {
        if (e1.head) e2 = h.wrapValue(e1, e2)
        if (e2.head) e1 = h.wrapValue(e2, e1)
        switch (op) {
          case '=':
          case '==':
            expr = B.core.rel.eq([e1, e2])
            break
          case '!=':
            expr = B.core.rel.neq([e1, e2])
            break
          case '>':
            expr = B.core.rel.gr([e1, e2])
            break
          case '<':
            expr = B.core.rel.lt([e1, e2])
            break
          case '>=':
            expr = B.core.rel.gre([e1, e2])
            break
          case '<=':
            expr = B.core.rel.lte([e1, e2])
            break
          default: throw new Error(`value operator '${op}' not supported`);
        }
      }
      return B.struct.generator.atomGroups({ 'atom-test': expr })
    }
  }
]

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
      r.RangeListProperty,
      r.ValueQuery,
      r.Keywords,
    )
  },

  Keywords: () => P.alt(...h.getKeywordRules(keywords)),

  ValueRange: function(r) {
    return P.seq(
      r.Value
        .skip(P.regex(/\s+TO\s+/i)),
      r.Value
    ).map(x => ({range: x}))
  },

  RangeListProperty: function(r) {
    return P.seq(
      P.alt(...h.getPropertyNameRules(properties, /\s/))
        .skip(P.whitespace),
      P.alt(
        r.ValueRange,
        r.Value
      ).sepBy1(P.whitespace)
    ).map(x => {
      const [property, values] = x
      const listValues: (string|number)[] = []
      const rangeValues: any[] = []

      values.forEach(v => {
        if (v.range) {
          rangeValues.push(
            B.core.rel.inRange([property, v.range[0], v.range[1]])
          )
        } else {
          listValues.push(h.wrapValue(property, v, sstrucDict))
        }
      })

      const rangeTest = h.orExpr(rangeValues)
      const listTest = h.valuesTest(property, listValues)

      let test
      if (rangeTest && listTest) {
        test = B.core.logic.or([rangeTest, listTest])
      } else {
        test = rangeTest ? rangeTest : listTest
      }

      return B.struct.generator.atomGroups({ [h.testLevel(property)]: test })
    })
  },

  Operator: function(r) {
    return h.combineOperators(operators, P.alt(r.Parens, r.Expression, r.ValueQuery))
  },

  Query: function(r) {
    return P.alt(
      r.Operator,
      r.Parens,
      r.Expression
    ).trim(P.optWhitespace)
  },

  Number: function () {
    return P.regexp(/-?(0|[1-9][0-9]*)([.][0-9]+)?([eE][+-]?[0-9]+)?/)
      .map(Number)
      .desc('number')
  },

  String: function () {
    const w = h.getReservedWords(properties, keywords, operators)
      .sort(h.strLenSortFn).map(h.escapeRegExp).join('|')
    return P.alt(
      P.regex(new RegExp(`(?!(${w}))[A-Z0-9_]+`, 'i')),
      P.regex(/'((?:[^"\\]|\\.)*)'/, 1),
      P.regex(/"((?:[^"\\]|\\.)*)"/, 1).map(x => B.core.type.regex([`^${x}$`, 'i']))
    ).desc('string')
  },

  Value: function (r) {
    return P.alt(r.Number, r.String)
  },

  ValueParens: function (r) {
    return P.alt(
      r.ValueParens,
      r.ValueOperator,
      r.ValueExpressions
    ).wrap(P.string('('), P.string(')'))
  },

  ValuePropertyNames: function() {
    return P.alt(...h.getPropertyNameRules(properties, /=~|==|>=|<=|=|!=|>|<|\)|\s|\+|-|\*|\//i))
  },

  ValueOperator: function(r) {
    return h.combineOperators(valueOperators, P.alt(r.ValueParens, r.ValueExpressions))
  },

  ValueExpressions: function(r) {
    return P.alt(
      r.ValueFunctions,
      r.Value,
      r.ValuePropertyNames
    )
  },

  ValueFunctions: function(r) {
    return P.alt(...h.getFunctionRules(functions, r.ValueOperator))
  },

  ValueQuery: function(r) {
    return P.alt(
      r.ValueOperator.map(x => {
        // if (!x.head || x.head.startsWith('core.math') || x.head.startsWith('structure.atom-property')) {
        if (!x.head || !x.head.startsWith('structure.generator')) {
          throw new Error(`values must be part of an comparison, value '${x}'`);
        } else {
          return x
        }
      })
    )
  }
})

const transpiler: Transpiler = str => lang.Query.tryParse(str)
export default transpiler
