export default /* GraphQL */ `
query AssemblySymmetry($assembly_id: String!, $entry_id: String!) {
    assembly(assembly_id: $assembly_id, entry_id: $entry_id) {
        rcsb_struct_symmetry {
            clusters {
                avg_rmsd
                members {
                    asym_id
                    pdbx_struct_oper_list_ids
                }
            }
            kind
            oligomeric_state
            rotation_axes {
                order
                start
                end
            }
            stoichiometry
            symbol
            type
        }
    }
}
`;