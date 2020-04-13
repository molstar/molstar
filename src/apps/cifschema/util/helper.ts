/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export type Import = { save?: string, file?: string }

export function parseImportGet(s: string): Import[] {
    // [{'save':hi_ang_Fox_coeffs  'file':templ_attr.cif}   {'save':hi_ang_Fox_c0  'file':templ_enum.cif}]
    // [{"file":'templ_enum.cif' "save":'H_M_ref'}]
    return s.trim().substring(2, s.length - 2).split(/}[ \n\t]*{/g).map(s => {
        const save = s.match(/('save'|"save"):([^ \t\n]+)/);
        const file = s.match(/('file'|"file"):([^ \t\n]+)/);
        return {
            save: save ? save[0].substr(7).replace(/['"]/g, '') : undefined,
            file: file ? file[0].substr(7).replace(/['"]/g, '') : undefined
        };
    });
}