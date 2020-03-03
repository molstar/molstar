/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// import { DataSourceProvider } from './provider';
// import { MmcifFormatProvider } from './formats';
// import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
// import { Download } from './actions';

// TODO: basic types string, binary (both download and from file)
// TODO: decompress functionality (string|binary --> string|binary)

// export const PDBeUpdatedMmcifDataSource = DataSourceProvider({
//     id: 'pdbe-updated-mmcif',
//     display: { name: 'PDBe Updated mmCIF', group: 'Molecule' },
//     format: MmcifFormatProvider,
//     params() {
//         return { id: PD.Text('1cbs', { label: 'PDB Id(s)', description: 'One or more comma separated PDB ids.' }) }
//     },
//     async apply(ctx, params) {
//         const download = Download(ctx.state.build().toRoot(), { url: `https://www.ebi.ac.uk/pdbe/static/entry/${params.id.toLowerCase()}_updated.cif` }, ctx);
//         await ctx.state.updateTree(download).runInContext(ctx.ctx);
//         return 0;
//     }
// })