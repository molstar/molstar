/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Tokens } from '../../../mol-io/reader/common/text/tokenizer';
import { CifCategory, CifField } from '../../../mol-io/reader/cif';
import { WaterNames } from '../../../mol-model/structure/model/types';

const Spec = {
    'MOL_ID': '',
    'MOLECULE': '',
    'CHAIN': '',
    'FRAGMENT': '',
    'SYNONYM': '',
    'EC': '',
    'ENGINEERED': '',
    'MUTATION': '',
    'OTHER_DETAILS': ''
}
type Spec = keyof typeof Spec

type Compound = { chains: string[], name: string }

export function parseCmpnd(lines: Tokens, lineStart: number, lineEnd: number) {
    const getLine = (n: number) => lines.data.substring(lines.indices[2 * n], lines.indices[2 * n + 1])

    let currentSpec: Spec | undefined
    let currentCompound: Compound = { chains: [], name: '' }
    const Compounds: Compound[] = []

    for (let i = lineStart; i < lineEnd; i++) {
        let line = getLine(i)
        // COLUMNS       DATA TYPE       FIELD         DEFINITION
        // ----------------------------------------------------------------------------------
        //  1 -  6       Record name     "COMPND"
        //  8 - 10       Continuation    continuation  Allows concatenation of multiple records.
        // 11 - 80       Specification   compound      Description of the molecular components.
        //               list

        const cmpnd = line.substr(10, 70).trim()
        const cmpndSpecEnd = cmpnd.indexOf(':')
        const cmpndSpec = cmpnd.substring(0, cmpndSpecEnd)

        let value: string

        if (cmpndSpec in Spec) {
            currentSpec = cmpndSpec as Spec
            value = cmpnd.substring(cmpndSpecEnd + 2)
        } else {
            value = cmpnd
        }
        value = value.replace(/;$/, '')

        if (currentSpec === 'MOL_ID') {
            currentCompound = {
                chains: [],
                name: ''
            }
            Compounds.push(currentCompound)
        } else if (currentSpec === 'MOLECULE') {
            if (currentCompound.name) currentCompound.name += ' '
            currentCompound.name += value
        } else if (currentSpec === 'CHAIN') {
            Array.prototype.push.apply(currentCompound.chains, value.split(/\s*,\s*/))
        }
    }

    return Compounds
}

export class EntityBuilder {
    private count = 0
    private ids: string[] = []
    private types: string[] = []
    private descriptions: string[] = []

    private compoundsMap = new Map<string, string>()
    private heteroMap = new Map<string, string>()
    private chainMap = new Map<string, string>()
    private waterId?: string

    private set(type: string, description: string) {
        this.count += 1
        this.ids.push(`${this.count}`)
        this.types.push(type)
        this.descriptions.push(description)
    }

    getEntityId(residueName: string, chainId: string, isHet: boolean): string {
        if (isHet) {
            if (WaterNames.has(residueName)) {
                if (this.waterId === undefined) {
                    this.set('water', 'Water')
                    this.waterId = `${this.count}`
                }
                return this.waterId;
            } else {
                if (!this.heteroMap.has(residueName)) {
                    this.set('non-polymer', residueName)
                    this.heteroMap.set(residueName, `${this.count}`)
                }
                return this.heteroMap.get(residueName)!
            }
        } else if (this.compoundsMap.has(chainId)) {
            return this.compoundsMap.get(chainId)!
        } else {
            if (!this.chainMap.has(chainId)) {
                this.set('polymer', chainId)
                this.chainMap.set(chainId, `${this.count}`)
            }
            return this.chainMap.get(chainId)!
        }
    }

    getEntityCategory() {
        const entity = {
            id: CifField.ofStrings(this.ids),
            type: CifField.ofStrings(this.types),
            pdbx_description: CifField.ofStrings(this.descriptions)
        }
        return CifCategory.ofFields('entity', entity)
    }

    setCompounds(compounds: Compound[]) {
        for (let i = 0, il = compounds.length; i < il; ++i) {
            const { chains, name } = compounds[i]
            this.set('polymer', name)
            for (let j = 0, jl = chains.length; j < jl; ++j) {
                this.compoundsMap.set(chains[j], `${this.count}`)
            }
        }
    }
}