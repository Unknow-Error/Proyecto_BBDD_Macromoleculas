import click

from utils import features_search as fs
from utils import pdb_search as pdb
from utils import pdb_viewer as pdbv
from utils import prote_search as ps
from utils import rmsd_analysis as rmsd


# CLI para buscar proteínas en bases de datos biológicas
@click.group()
def cli():
    pass


# Busca información de una proteína por su ID
@cli.command()
@click.argument("prompt")
def buscar(prompt):
    print(ps.buscar(prompt))


# Busca estructuras PDB asociadas a un accession de UniProt
@cli.command()
@click.argument("accession")
def buscar_pdb(accession):
    print(pdb.lista_pdb(accession))


# Analiza RMSD local entre dos estructuras PDB
@cli.command()
@click.argument("pdb1")
@click.argument("pdb2")
@click.option("--cadena1", "-c1", help="ID de la cadena del PDB1 a analizar (opcional)")
@click.option("--cadena2", "-c2", help="ID de la cadena del PDB2 a analizar (opcional)")
@click.option(
    "--ventana", "-w", default=5, help="Tamaño de la ventana deslizante (default: 5)"
)
def rmsd_pdb(pdb1, pdb2, cadena1, cadena2, ventana):
    resultado = rmsd.analizar_rmsd_local(pdb1, pdb2, cadena1, cadena2, ventana)
    if resultado[0]:
        print(f"\nAnálisis completado exitosamente!")
        print(f"Archivo generado: {resultado[0]}")


# Descarga features de una proteína por accession de UniProt
@cli.command()
@click.argument("accession")
@click.option(
    "--formato", "-f", default="json", help="Formato de salida (json, txt, xml, gff)"
)
def features(accession, formato):
    print(fs.descargar_features(accession, formato))


# Mostrar la estructura de la proteína según su pdb.
@cli.command()
@click.argument("codigopdb")
@click.option("--cadena", "-c", help="ID de la cadena a analizar (opcional)")
def mostrar_PDB_simple(codigopdb, cadena):
    pdbMostrador = pdbv.PDB_Viewer(codigopdb)
    pdbMostrador.mostrar_pdb_desde_id(cadena)


@cli.command()
@click.argument("codigopdb")
@click.option(
    "--feature",
    "-f",
    nargs=4,
    type=(str, int, int, str),
    help="Color, comienzo y fin, y nombre de la feature a añadir. Color hex, inicio, fin, nombre. Ej: -f '#FF0000' 10 50 dominio",
)
def mostrar_PDB_features(codigopdb, feature):
    pdbMostrador = pdbv.PDB_Viewer(codigopdb)
    pdbMostrador.mostrar_pdb_domains_regiones(feature)


# Mostrar el alineamiento entre dos proteínas según sus códigos pdb.
@cli.command()
@click.argument("codigopdb1")
@click.argument("codigopdb2")
@click.option("--cadena", "-c", help="ID de la cadena a analizar (opcional)")
@click.option(
    "--ventana", "-w", default=5, help="Tamaño de la ventana deslizante (default: 5)"
)
@click.option(
    "--colores",
    "-col",
    nargs=2,
    type=(str, str),
    help="Colores de las cadenas alineadas. Color_1 Color_2. Ej: blue green (Opcional)",
)
def mostrar_alineamiento(codigopdb1, codigopdb2, cadena, colores, ventana):
    pdbMostrador = pdbv.PDB_Viewer(codigopdb1)
    pdbMostrador.mostrar_alineamiento_pdb(codigopdb2, cadena, colores, ventana)


if __name__ == "__main__":
    cli()
