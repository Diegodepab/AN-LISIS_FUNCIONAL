import pandas as pd
import requests
import warnings
import gzip
import shutil
import networkx as nx
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from collections import defaultdict

# Leer el archivo TSV (dataset original que contiene un listado de los genes diferencialmente expresados que has obtenido previamente)
try:
    df = pd.read_csv('datos.tsv', sep='\t')
except pd.errors.ParserError as e:
    warnings.warn(f"Error al analizar el archivo: {e}")

# Filtrar genes sobreexpresados (logFC > 0 y adj.P.Val < 0.05)
overexpressed_genes = df[(df['logFC'] > 0) & (df['adj.P.Val'] < 0.05)]['ID'].tolist()

# Función para obtener identificador de STRING a partir del ID del gen
def get_string_id(gene_id):
    url = f"https://string-db.org/api/json/get_string_ids?identifiers={gene_id}&species=3702"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if data:
            # Devolver el primer ID de STRING encontrado
            return data[0]['stringId'] if 'stringId' in data[0] else None
    else:
        warnings.warn(f"No se encontró identificador de STRING para {gene_id}.")
    return None

# Obtener los identificadores de STRING para los genes sobreexpresados
gene_string_ids = {}
for gene in overexpressed_genes:
    string_id = get_string_id(gene)
    if string_id:
        gene_string_ids[gene] = string_id


# Definir el organismo (Arabidopsis thaliana) para STRING, 
organism_taxon_id = 3702

# URL base de STRING para redes de interacción proteica
base_url = "https://stringdb-downloads.org/download/protein.links.v12.0/"

# Nombre del archivo de la red completa para el organismo
file_name = f"{organism_taxon_id}.protein.links.v12.0.txt.gz"

# URL completa 
url = base_url + file_name

# Descargar de la red de arabiidopsis thaliana
print(f"Descargando datos de la red de interacción proteica para el organismo con taxon ID: {organism_taxon_id} (Arabidopsis Thaliana)")
response = requests.get(url, stream=True)

# Guardar el archivo comprimido localmente
with open(file_name, 'wb') as f:
    f.write(response.content)

# Descomprimir el archivo descargado
output_file = file_name.replace('.gz', '')
print("Descomprimiendo archivo...")
with gzip.open(file_name, 'rb') as f_in:
    with open(output_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

print(f"Archivo descomprimido: {output_file}")

# Dataframe que contiene las relaciones en el organismo buscado 
network_df = pd.read_csv(output_file, sep=' ')

# Filtrar por combined score escogí 400, pero 700 o 900 pueden ser valores interesantes para estudiar relaciones más directas
filtered_network_df = network_df[network_df['combined_score'] > 400]

# Obtener los STRING IDs del DataFrame filtrado
filtered_proteins = set(filtered_network_df['protein1']).union(set(filtered_network_df['protein2']))

# Comprobar si los STRING IDs están en la red filtrada y emitir advertencias si no se encuentran
for gene, string_id in gene_string_ids.items():
    if string_id not in filtered_proteins:
        warnings.warn(f"El STRING ID {string_id} para el gen {gene} no se encontró en la red filtrada.")
    
    

#----------------------DIAMOnD-------------------------
alpha= 1 

#funcion que crea el grafo a partir del dataframe estudiado
def read_input_from_df(df):
    G = nx.Graph()

    # Recorre el DataFrame fila por fila, agregando cada arista al grafo
    for index, row in df.iterrows():
        node1 = str(row.iloc[0])  # Acceder a la primera columna por posición
        node2 = str(row.iloc[1])  # Acceder a la segunda columna por posición
        G.add_edge(node1, node2)

    return G

################################### main##############################

# Leer el grafo desde el DataFrame
G = read_input_from_df(df)


#G, seed_genes, n, alpha
# G, seed_genes, diamond_genes


# ------------------- FUNCIONES DE DIAMOND -------------------

def compute_all_gamma_ln(N):
    """
    Calcula el logaritmo del valor gamma para cada número de 1 a N.
    El valor gamma se usa en las distribuciones hipergeométricas para evitar cálculos repetidos.
    """
    gamma_ln = {i: scipy.special.gammaln(i) for i in range(1, N + 1)}
    return gamma_ln

def gauss_hypergeom(x, r, b, n, gamma_ln):
    """
    Calcula el valor de la distribución hipergeométrica (logaritmo de la probabilidad)
    usando los valores precomputados de gamma_ln.
    """
    max_index = len(gamma_ln) - 1  # El último índice permitido
    if r + b > max_index or x > max_index or r < x or b < (n - x):
        # Si los índices están fuera de rango, retornar un valor predeterminado (0)
        return 0  # Este valor puede cambiar a 1 si deseas un valor de p-valor de máximo riesgo

    # Fórmula de la distribución hipergeométrica
    return np.exp(gamma_ln[r] - (gamma_ln[x] + gamma_ln[r - x]) +
                  gamma_ln[b] - (gamma_ln[n - x] + gamma_ln[b - n]) -
                  gamma_ln[r + b])

def pvalue(kb, k, N, s, gamma_ln):
    """
    Calcula el p-valor sumando los resultados de la distribución hipergeométrica
    para cada valor en el rango [kb, k].
    """
    return sum(gauss_hypergeom(n, s, N - s, k, gamma_ln) for n in range(kb, k + 1))

def diamond_iteration_of_first_X_nodes(G, S, X, alpha=1):
    """
    Realiza la iteración DIAMOnD sobre los primeros X nodos seleccionados, basándose en un conjunto de nodos semilla S.
    Utiliza el p-valor para seleccionar los nodos de forma iterativa.
    """
    added_nodes = []  # Lista para almacenar los nodos seleccionados
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}  # Diccionario de vecinos de cada nodo
    degrees = dict(G.degree())  # Obtener grados de nodos (número de conexiones para cada nodo)
    cluster_nodes = set(S)  # Conjunto de nodos semilla
    gamma_ln = compute_all_gamma_ln(len(G.nodes))  # Precomputar valores de logaritmo de gamma

    # Iterar hasta seleccionar X nodos
    while len(added_nodes) < X:
        min_p = float('inf')  # Inicializar el p-valor más bajo
        next_node = None  # Nodo a añadir en la siguiente iteración
        for node in set(G.nodes) - cluster_nodes:  # Iterar sobre nodos que no están en el cluster
            k = degrees[node]  # Grado del nodo (número de conexiones)
            kb = sum((1 for neighbor in neighbors[node] if neighbor in cluster_nodes))  # Grado en el cluster de semillas

            # Evitar divisiones por cero o nodos sin conexiones
            if k == 0:
                continue

            try:
                # Calcular el p-valor para este nodo
                p = pvalue(kb, k, len(G.nodes), len(cluster_nodes), gamma_ln)
            except ValueError as e:
                print(f"Error calculando pvalue para nodo {node}: {e}")
                continue

            # Si el p-valor es más bajo, seleccionar este nodo
            if p < min_p:
                min_p = p
                next_node = node

        # Añadir el nodo con el p-valor mínimo
        if next_node:
            added_nodes.append(next_node)
            cluster_nodes.add(next_node)  # Añadir a los nodos seleccionados

    return added_nodes  # Retornar la lista de nodos seleccionados

# ------------------- VISUALIZACIÓN -------------------

def graficar_red_enriquecida(G, seed_genes, diamond_genes):
    """
    Grafica la red enriquecida con los nodos semilla y los nodos seleccionados por DIAMOnD.
    Los nodos semilla se representan en azul y los nodos DIAMOnD en naranja.
    """
    pos = nx.spring_layout(G)  # Disposición de los nodos en el espacio (layout de la red)
    plt.figure(figsize=(10, 10))  # Tamaño de la figura

    # Nodos de semillas (color azul claro)
    nx.draw_networkx_nodes(G, pos, nodelist=seed_genes, node_color='lightblue', node_size=500, label="Seed Genes")

    # Nodos DIAMOnD (color naranja)
    nx.draw_networkx_nodes(G, pos, nodelist=diamond_genes, node_color='orange', node_size=300, label="DIAMOnD Genes")

    # Enlaces entre nodos (conjunto de nodos semilla y DIAMOnD)
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v in G.edges() if u in (seed_genes | set(diamond_genes)) and v in (seed_genes | set(diamond_genes))], alpha=0.5)

    # Etiquetas de los nodos
    nx.draw_networkx_labels(G, pos, font_size=8)

    # Mostrar leyenda, título y la figura
    plt.legend(loc="best")
    plt.title("Red de Genes Enriquecida usando DIAMOnD")
    plt.show()

# ------------------- EJECUCIÓN -------------------

# Ejemplo de los datos cargados y creación de la red
# seed_genes = seed_genes_grado  # Asignar el conjunto de genes semilla

# Realizar la iteración DIAMOnD para obtener los genes seleccionados
diamond_genes = diamond_iteration_of_first_X_nodes(G, seed_genes, n, alpha)

# Graficar la red enriquecida
graficar_red_enriquecida(G, seed_genes, diamond_genes)



