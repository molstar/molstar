/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Column from '../../../mol-base/collections/column'

export type Table<Data> = { [E in keyof Data]: Column<Data[E]> }

export interface ElementSymbol extends String { '@type': 'element-symbol' }
export function ElementSymbol(s: string): ElementSymbol {
    // TODO: optimize?
    return s.toUpperCase() as any;
}


export interface Atoms extends Table<{
    name: string,
    elementSymbol: ElementSymbol,
    
}> { }
 
interface Common {

}

export default Common