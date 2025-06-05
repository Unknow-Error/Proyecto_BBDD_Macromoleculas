def main():
    print("=== BUSCADOR DE PROTEÍNAS ===\n")
    user_input = input("Ingrese un identificador (UniProt o NCBI RefSeq): ").strip()

    tipo, datos = identify_protein(user_input)
    print(f"\nTipo detectado: {tipo}")

    if not datos:
        print("No se pudieron obtener datos.")
        return

    print(f"\nIdentificador UniProt: {datos.get('uniprot_id', '-')}")
    print(f"Identificador NCBI: {datos.get('ncbi_id', '-')}")
    print(f"\nFeatures encontradas: {len(datos['features'])}")
    for f in datos['features'][:10]:  # Mostrar solo las primeras 10
        print(f"- {f['type']}: {f['description']} ({f['start']} - {f['end']})")

    print(f"\nPDBs asociados ({len(datos['pdb_ids'])}): {', '.join(datos['pdb_ids']) if datos['pdb_ids'] else 'Ninguno'}")


if __name__ == "__main__":
    main()
