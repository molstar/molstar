import Examples from 'mol-script/script/mol-script/examples'
import parse from 'mol-script/script/mol-script/parser'

const expr = parse(Examples[Examples.length - 1].value);
console.log(expr);