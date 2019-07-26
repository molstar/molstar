/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement } from '../../../mol-model/structure';
import { EmptyLoci, isEmptyLoci } from '../../../mol-model/loci';
import { MolScriptBuilder } from '../../../mol-script/language/builder';
import { formatMolScript } from '../../../mol-script/language/expression-formatter';

export function getExpression(loci: StructureElement.Loci | EmptyLoci) {
    const scriptExpression = isEmptyLoci(loci)
        ? MolScriptBuilder.struct.generator.empty()
        : StructureElement.Loci.toScriptExpression(loci)
    return formatMolScript(scriptExpression)
}