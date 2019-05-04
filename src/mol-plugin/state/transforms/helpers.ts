/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure, StructureSelection, QueryContext } from 'mol-model/structure';
import { Color } from 'mol-util/color';
import { Overpaint } from 'mol-theme/overpaint';
import { parseMolScript } from 'mol-script/language/parser';
import { transpileMolScript } from 'mol-script/script/mol-script/symbols';
import { compile } from 'mol-script/runtime/query/compiler';
import { Transparency } from 'mol-theme/transparency';
import { ComputedSecondaryStructure } from 'mol-model-props/computed/secondary-structure';

type Script = { language: string, expression: string }

function scriptToLoci(structure: Structure, script: Script) {
    const parsed = parseMolScript(script.expression)
    if (parsed.length === 0) throw new Error('No query')
    const query = transpileMolScript(parsed[0])

    const compiled = compile<StructureSelection>(query)
    const result = compiled(new QueryContext(structure))
    return StructureSelection.toLoci2(result)
}

export function getStructureOverpaint(structure: Structure, scriptLayers: { script: Script, color: Color }[], alpha: number): Overpaint {
    const layers: Overpaint.Layer[] = []
    for (let i = 0, il = scriptLayers.length; i < il; ++i) {
        const { script, color } = scriptLayers[i]
        layers.push({ loci: scriptToLoci(structure, script), color })
    }
    return { layers, alpha }
}

export function getStructureTransparency(structure: Structure, script: Script, value: number, variant: Transparency.Variant): Transparency {
    return { loci: scriptToLoci(structure, script), value, variant }
}

/**
 * Attaches ComputedSecondaryStructure property when unavailable in sourceData
 */
export async function ensureSecondaryStructure(s: Structure) {
    if (s.model && s.model.sourceData.kind === 'mmCIF') {
        if (!s.model.sourceData.data.struct_conf.id.isDefined && !s.model.sourceData.data.struct_sheet_range.id.isDefined) {
            await ComputedSecondaryStructure.attachFromCifOrCompute(s)
        }
    }
}