import click

from utils import prote_search as ps


@click.group()
def cli():
    pass


@cli.command()
@click.argument("prompt")
def buscar(prompt):
    print(ps.buscar(prompt))


@cli.command()
@click.argument("prompt")
def descargar(id):
    print(ps.buscar(id))


if __name__ == "__main__":
    cli()
