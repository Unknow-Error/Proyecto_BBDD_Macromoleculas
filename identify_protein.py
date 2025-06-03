def identify_protein(protein):
    if protein[1] == "P" and protein[2] == "_":
        uniprot_ids = map_ncbi_to_uni(protein)
        if not uniprot_ids:
            return "ID NCBI sin correspondencia en UniProt", {}

        for id in uniprot_ids:
          features = fetch_features(id)
          pdbs = fetch_structures(id)
          return "Identificador de NCBI", {
              "ncbi_id": protein,
              "uniprot_id": id,
              "features": features,
              "pdb_ids": pdbs
          }

    elif protein[1].isdigit():

        features = fetch_features(protein)
        pdbs = fetch_structures(protein)
        return "Identificador de UniProt", {
            "uniprot_id": protein,
            "features": features,
            "pdb_ids": pdbs
        }

    else:
        # Búsqueda por texto libre
        #text = fetch_text(protein)
        #return "Búsqueda por texto", text
        return "Búsqueda por texto", {}
