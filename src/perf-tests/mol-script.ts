import { MolScriptBuilder } from '../mol-script/language/builder';
import { compile, QuerySymbolRuntime, DefaultQueryRuntimeTable } from '../mol-script/runtime/query/compiler';
import { QueryContext, Structure, StructureQuery } from '../mol-model/structure';
import { readCifFile, getModelsAndStructure } from '../cli/structure-info/model';
import { CustomPropSymbol } from '../mol-script/language/symbol';
import Type from '../mol-script/language/type';
import { parseMolScript } from '../mol-script/language/parser';
import * as util from 'util';
import { transpileMolScript } from '../mol-script/script/mol-script/symbols';
import { formatMolScript } from '../mol-script/language/expression-formatter';
import { StructureQualityReport, StructureQualityReportProvider } from '../extensions/pdbe/structure-quality-report/prop';
import fetch from 'node-fetch';
import { CustomPropertyDescriptor } from '../mol-model/custom-property';

// import Examples from 'mol-script/script/mol-script/examples'
// import { parseMolScript } from 'mol-script/script/mol-script/parser'
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

// ;; :residue-test (= atom.label_comp_id REA)
const exprs = parseMolScript(`(sel.atom.atom-groups
    :residue-test (> pdbe.structure-quality.issue-count 0)
    :atom-test (= atom.el _C))`);

const tsp = transpileMolScript(exprs[0]);

// console.log(util.inspect(exprs, true, 10, true));

console.log(util.inspect(tsp, true, 10, true));

console.log(formatMolScript);
console.log(formatMolScript(tsp));
// //console.log(expr);

const expr = MolScriptBuilder.core.math.add([1, 2, 3]);
const compiled = compile<number>(expr);
const result = compiled(new QueryContext(Structure.Empty));
console.log(result);

const CustomProp = CustomPropertyDescriptor({
    name: 'test_prop',
    cifExport: { prefix: '', categories: [ ]},
    symbols: {
        residueIndex: QuerySymbolRuntime.Dynamic(CustomPropSymbol('custom.test-prop', 'residue-index', Type.Num), ctx => {
            const e = ctx.element;
            // console.log(e.element, e.unit.model.atomicHierarchy.residueAtomSegments.index[e.element])
            return e.unit.model.atomicHierarchy.residueAtomSegments.index[e.element];
        })
    }
});

DefaultQueryRuntimeTable.addCustomProp(CustomProp);

DefaultQueryRuntimeTable.addCustomProp(StructureQualityReportProvider.descriptor);

export async function testQ() {
    const frame = await readCifFile('e:/test/quick/1cbs_updated.cif');
    const { structure } = await getModelsAndStructure(frame);
    const model = structure.models[0];

    const rawData = await fetch(`https://www.ebi.ac.uk/pdbe/api/validation/residuewise_outlier_summary/entry/${model.entryId.toLowerCase()}`, { timeout: 1500 });
    const data = StructureQualityReport.fromJson(model, await rawData.json());

    StructureQualityReportProvider.set(model, { serverUrl: '' }, data);

    let expr = MolScriptBuilder.struct.generator.atomGroups({
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

    expr = tsp;

    const compiled = compile<StructureQuery>(expr);
    const result = compiled(new QueryContext(structure));

    console.log(result);
}

testQ();