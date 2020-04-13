/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface Writer {
    writeString(data: string): boolean,
    writeBinary(data: Uint8Array): boolean
}

namespace Writer {

}

export default Writer;