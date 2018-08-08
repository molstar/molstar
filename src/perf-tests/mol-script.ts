import { MolScriptBuilder } from 'mol-script/language/builder';
import { compile, QuerySymbolRuntime, DefaultQueryRuntimeTable } from 'mol-script/runtime/query/compiler';
import { QueryContext, Structure, StructureQuery, ModelPropertyDescriptor } from 'mol-model/structure';
import { readCifFile, getModelsAndStructure } from '../apps/structure-info/model';
import { CustomPropSymbol } from 'mol-script/language/symbol';
import Type from 'mol-script/language/type';

// import Examples from 'mol-script/script/mol-script/examples'
// import { parseMolScript } from 'mol-script/script/mol-script/parser'
// import * as util from 'util'
// //import { compileAST } from 'mol-script/script/mol-script/compile';

// for (const e of Examples) {
//     const expr = parseMolScript(e.value)[0];
//     console.log(e.name, util.inspect(expr, true, 10, true));
// }
// const exprs = parseMolScript(`(sel.atom.atom-groups
//     :residue-test (= atom.auth_comp_id ALA)
//     ;; ho ho ho
//     :atom-test (set.has { _C _N } atom.el)) ; comm
//     ;; this is a comment
//     ((hi) (ho))`);

// console.log(util.inspect(exprs, true, 10, true));
// //console.log(expr);

const expr = MolScriptBuilder.core.math.add([1, 2, 3]);
const compiled = compile<number>(expr);
const result = compiled(new QueryContext(Structure.Empty));
console.log(result);

const CustomProp = ModelPropertyDescriptor({
    name: 'test_prop',
    isStatic: true,
    cifExport: { categories: [ ]},
    symbols: {
        residueIndex: QuerySymbolRuntime.Dynamic(CustomPropSymbol('custom.test-prop', 'residue-index', Type.Num), ctx => {
            const e = ctx.element;
            //console.log(e.element, e.unit.model.atomicHierarchy.residueAtomSegments.index[e.element])
            return e.unit.model.atomicHierarchy.residueAtomSegments.index[e.element];
        })
    }
});

DefaultQueryRuntimeTable.addCustomProp(CustomProp);

async function testQ() {
    const frame = await readCifFile('e:/test/quick/1cbs_updated.cif');
    const { structure } = await getModelsAndStructure(frame);

    const expr = MolScriptBuilder.struct.generator.atomGroups({
        'atom-test': MolScriptBuilder.core.rel.eq([
            MolScriptBuilder.struct.atomProperty.core.elementSymbol(),
            MolScriptBuilder.es('C')
        ]),
        // 'residue-test': MolScriptBuilder.core.rel.eq([
        //     MolScriptBuilder.struct.atomProperty.macromolecular.label_comp_id(),
        //     'REA'
        // ])
        'residue-test': MolScriptBuilder.core.rel.inRange([CustomProp.symbols.residueIndex.symbol(), 1, 5])
    });

    const compiled = compile<StructureQuery>(expr);
    const result = compiled(new QueryContext(structure));

    console.log(result);
}

testQ();