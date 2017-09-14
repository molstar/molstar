/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { File } from './file'
import { TextBlock } from './text-block'

export class TextFile implements File<TextBlock> {
    data: string;
    blocks: TextBlock[] = [];

    constructor(data: string) {
        this.data = data;
    }
}
