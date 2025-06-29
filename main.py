import click

from utils import features_search as fs
from utils import pdb_search as pdb
from utils import prote_search as ps
from utils import rmsd_analysis as rmsd
from utils import pdb_viewer as pdbv

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
@click.option("--cadena", "-c", help="ID de la cadena a analizar (opcional)")
@click.option(
    "--ventana", "-w", default=5, help="Tamaño de la ventana deslizante (default: 5)"
)
def rmsd_pdb(pdb1, pdb2, cadena, ventana):
    resultado = rmsd.analizar_rmsd_local(pdb1, pdb2, cadena, ventana)
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
def mostrar_PDB(codigopdb):
    pdbMostrador = pdbv.PDB_Viewer(codigopdb)
    pdbMostrador.mostrar_pdb_domains_regiones()

if __name__ == "__main__":
    cli()