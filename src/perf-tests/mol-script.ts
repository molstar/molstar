import Examples from 'mol-script/parsers/mol-script/examples'
import transpile from 'mol-script/parsers/mol-script/parser'

const expr = transpile(Examples[Examples.length - 1].value);
console.log(expr);