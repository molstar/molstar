/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Category } from "../../../../mol-io/writer/cif/encoder";
import { LigandExplorer } from "../ligand-encoder";
import { StringBuilder } from "../../../../mol-util";

// specification: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
export class Mol2Encoder extends LigandExplorer {
    private meta: StringBuilder;
    private encoded = false;
    private error = false;

    writeCategory<Ctx>(category: Category<Ctx>, context?: Ctx): void {
        if (this.encoded) {
            throw new Error('The writer contents have already been encoded, no more writing.');
        }

    }

    encode(): void {
        // write meta-information, do so after ctab
        if (this.error || this.metaInformation) {
            StringBuilder.writeSafe(this.builder, StringBuilder.getString(this.meta));
        }

        this.encoded = true;
    }

    constructor(readonly encoder: string, readonly metaInformation: boolean, readonly hydrogens: boolean) {
        super(encoder, hydrogens);
        this.meta = StringBuilder.create();
    }
}