#!/usr/bin/env python3
"""
Script de prueba para la funcionalidad de búsqueda de PDB
"""

from data.fetch_pdb import buscar_pdb_por_accession


def test_buscar_pdb():
    """Prueba la búsqueda de PDB con diferentes accessions"""

    # Test 1: Accession válido de UniProt
    print("=== Test 1: Accession válido ===")
    accession1 = "P01308"  # Insulina humana
    resultado1 = buscar_pdb_por_accession(accession1)
    print(f"Resultado: {resultado1}")
    print()

    # Test 2: Accession inválido
    print("=== Test 2: Accession inválido ===")
    accession2 = "INVALID123"
    resultado2 = buscar_pdb_por_accession(accession2)
    print(f"Resultado: {resultado2}")
    print()

    # Test 3: Accession que no existe
    print("=== Test 3: Accession inexistente ===")
    accession3 = "P99999"
    resultado3 = buscar_pdb_por_accession(accession3)
    print(f"Resultado: {resultado3}")


if __name__ == "__main__":
    test_buscar_pdb()
