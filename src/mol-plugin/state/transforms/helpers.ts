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

type Script = { language: string, expression: string }

export function getStructureOverpaint(structure: Structure, scriptLayers: { script: Script, color: Color }[], alpha: number): Overpaint {
    const layers: Overpaint.Layer[] = []
    for (let i = 0, il = scriptLayers.length; i < il; ++i) {
        const { script, color } = scriptLayers[i]
        const parsed = parseMolScript(script.expression)
        if (parsed.length === 0) throw new Error('No query')
        const query = transpileMolScript(parsed[0])

        const compiled = compile<StructureSelection>(query)
        const result = compiled(new QueryContext(structure))
        const loci = StructureSelection.toLoci2(result)

        layers.push({ loci, color })
    }
    return { layers, alpha }
}

export function getStructureTransparency(structure: Structure, scriptLayers: { script: Script, value: number }[]): Transparency {
    const layers: Transparency.Layer[] = []
    for (let i = 0, il = scriptLayers.length; i < il; ++i) {
        const { script, value } = scriptLayers[i]
        const parsed = parseMolScript(script.expression)
        if (parsed.length === 0) throw new Error('No query')
        const query = transpileMolScript(parsed[0])

        const compiled = compile<StructureSelection>(query)
        const result = compiled(new QueryContext(structure))
        const loci = StructureSelection.toLoci2(result)

        layers.push({ loci, value })
    }
    return { layers }
}