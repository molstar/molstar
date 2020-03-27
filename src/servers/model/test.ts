// import { resolveJob } from './server/query';
// import * as fs from 'fs'
// import * as path from 'path'
// import { StructureCache } from './server/structure-wrapper';
// import { createJob } from './server/jobs';

// function wrapFile(fn: string) {
//     const w = {
//         open(this: any) {
//             if (this.opened) return;
//             this.file = fs.openSync(fn, 'w');
//             this.opened = true;
//         },
//         writeBinary(this: any, data: Uint8Array) {
//             this.open();
//             fs.writeSync(this.file, Buffer.from(data));
//             return true;
//         },
//         writeString(this: any, data: string) {
//             this.open();
//             fs.writeSync(this.file, data);
//             return true;
//         },
//         end(this: any) {
//             if (!this.opened || this.ended) return;
//             fs.close(this.file, function () { });
//             this.ended = true;
//         },
//         file: 0,
//         ended: false,
//         opened: false
//     };

//     return w;
// }

// export const basePath = path.join(__dirname, '..', '..', '..', '..')
// export const examplesPath = path.join(basePath, 'examples')
// export const outPath = path.join(basePath, 'build', 'test')
// if (!fs.existsSync(outPath)) fs.mkdirSync(outPath);

// async function run() {
//     try {
//         // const testFile = '1crn.cif'
//         // const testFile = '1grm_updated.cif'
//         // const testFile = 'C:/Projects/mol-star/molstar/build/test/in/1grm_updated.cif'
//         // const request = createJob({
//         //     entryId: testFile,
//         //     queryName: 'full',
//         //     queryParams: { },
//         // });
//         const testFile = '1cbs_updated.cif'
//         const request = createJob({
//             entryId: path.join(examplesPath, testFile),
//             queryName: 'full',
//             queryParams: { }
//         });

//         // const request = createJob({
//         //     entryId: path.join(examplesPath, testFile),
//         //     queryName: 'atoms',
//         //     queryParams: {
//         //         atom_site: { label_comp_id: 'ALA' }
//         //     }
//         // });
//         // const request = createJob({
//         //     entryId: path.join(examplesPath, testFile),
//         //     queryName: 'residueInteraction',
//         //     queryParams: {
//         //         atom_site: { label_comp_id: 'REA' },
//         //         radius: 5
//         //     }
//         // });
//         const encoder = await resolveJob(request);
//         const writer = wrapFile(path.join(outPath, testFile));
//         // const writer = wrapFile(path.join(outPath, '1grm_test.cif'));
//         encoder.writeTo(writer);
//         writer.end();
//     } catch (e) {
//         console.error(e)
//     } finally {
//         StructureCache.expireAll();
//     }
// }

// run();