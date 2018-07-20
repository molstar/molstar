import { resolveJob } from './server/query';
import * as fs from 'fs'
import { StructureCache } from './server/structure-wrapper';
import { createJob } from './server/jobs';

function wrapFile(fn: string) {
    const w = {
        open(this: any) {
            if (this.opened) return;
            this.file = fs.openSync(fn, 'w');
            this.opened = true;
        },
        writeBinary(this: any, data: Uint8Array) {
            this.open();
            fs.writeSync(this.file, new Buffer(data));
            return true;
        },
        writeString(this: any, data: string) {
            this.open();
            fs.writeSync(this.file, data);
            return true;
        },
        end(this: any) {
            if (!this.opened || this.ended) return;
            fs.close(this.file, function () { });
            this.ended = true;
        },
        file: 0,
        ended: false,
        opened: false
    };

    return w;
}

async function run() {
    try {
        const request = createJob('_local_', 'e:/test/quick/1cbs_updated.cif', 'residueInteraction', { label_comp_id: 'REA' });
        const encoder = await resolveJob(request);
        const writer = wrapFile('e:/test/mol-star/1cbs_full.cif');
        encoder.writeTo(writer);
        writer.end();
    } finally {
        StructureCache.expireAll();
    }
}

run();