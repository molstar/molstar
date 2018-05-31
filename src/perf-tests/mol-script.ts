import Examples from 'mol-script/script/mol-script/examples'
import { parseMolScript } from 'mol-script/script/mol-script/parser'
import * as util from 'util'
//import { compileAST } from 'mol-script/script/mol-script/compile';

for (const e of Examples) {
    const expr = parseMolScript(e.value)[0];
    console.log(e.name, util.inspect(expr, true, 10, true));
}
// const exprs = parseMolScript(`(sel.atom.atom-groups
//     :residue-test (= atom.auth_comp_id ALA)
//     ;; ho ho ho
//     :atom-test (set.has { _C _N } atom.el)) ; comm
//     ;; this is a comment
//     ((hi) (ho))`);

// console.log(util.inspect(exprs, true, 10, true));
// //console.log(expr);