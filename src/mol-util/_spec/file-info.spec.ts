/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Russell Parker <russell@benchling.com>
 */

import { getFileNameInfo } from '../file-info';

describe('getFileNameInfo', () => {
    it('handles empty string', () => {
        expect(getFileNameInfo('')).toEqual({ path: '', name: '', ext: '', base: '', dir: '', protocol: '', query: '' });
    });

    it('handles url', () => {
        expect(getFileNameInfo('https://models.rcsb.org/4KTC.bcif')).toEqual({ path: 'models.rcsb.org/4KTC.bcif', name: '4KTC.bcif', ext: 'bcif', base: '4KTC', dir: 'models.rcsb.org/', protocol: 'https', query: '' });
    });

    it('handles compressed url', () => {
        expect(getFileNameInfo('https://files.rcsb.org/download/7QPD.cif.gz?foo=bar')).toEqual({ path: 'files.rcsb.org/download/7QPD.cif.gz', name: '7QPD.cif.gz', ext: 'cif', base: '7QPD', dir: 'files.rcsb.org/download/', protocol: 'https', query: '?foo=bar' });
    });

    it('handles local path', () => {
        expect(getFileNameInfo('/usr/local/data/structure.pdb')).toEqual({ path: '/usr/local/data/structure.pdb', name: 'structure.pdb', ext: 'pdb', base: 'structure', dir: '/usr/local/data/', protocol: '', query: '' });
    });

    it('handles local path with protocol', () => {
        expect(getFileNameInfo('file:///usr/local/data/structure.pdb')).toEqual({ path: '/usr/local/data/structure.pdb', name: 'structure.pdb', ext: 'pdb', base: 'structure', dir: '/usr/local/data/', protocol: 'file', query: '' });
    });
});
