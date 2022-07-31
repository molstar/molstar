
import Transpiler from './transpilers/transpiler'
import _transpiler from './transpilers/all'
import compile from './runtime/query/compiler'
//import compile from './reference-implementation/molql/compiler'

const transpiler: {[index: string]: Transpiler} = _transpiler

const util = require('util')

function parse(lang: string, str: string) {
  try {
    const query = transpiler[lang](str)
    console.log(str)
    console.log(util.inspect(query, {depth: 20, color: true}))
    console.log('\n')
    compile(query)
  } catch (e) {
    console.log(str)
    console.log(e.message)
    console.log('\n')
  }
}

const [,,lang, str] = process.argv

if (lang && str) {
  parse(lang, str)
}
