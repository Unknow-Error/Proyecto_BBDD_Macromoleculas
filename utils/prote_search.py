import re  # re es una libreria para evaluar expreciones regulares

from data.fetch_ncbi import buscar_acn_ncbi
from data.fetch_uniprot import buscar_id_uniprot, buscar_uniprot


# identifica si el texto es un id de uniprot utilizando re y las expresiones de los id de uniprot
# Entrada = texto
# Salida = booleano
def es_id_uniprot(texto):
    return bool(
        re.fullmatch(
            r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[A-Z0-9]{7}$",
            texto,
        )
    )


# identifica si el texto es un id de ncbi utilizando re y las expresiones de los id de proteinas de ncbi (NP_, XP_, YP_, etc...)
# Entrada = texto
# Salida = booleano
def es_id_ncbi(texto):
    return bool(re.fullmatch(r"^[A-Z]{1,3}_?\d{5,9}(\.\d+)?$", texto))


def buscar(prompt):

    if es_id_uniprot(prompt):
        return buscar_id_uniprot(prompt)
    elif es_id_ncbi(prompt):
        return buscar_acn_ncbi(prompt)
    else:
        return buscar_uniprot(prompt)
