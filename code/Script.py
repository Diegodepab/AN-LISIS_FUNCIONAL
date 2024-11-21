############################################################
########## Descarga o llamar las librerías #################
############################################################

# Función para instalar una librería si no está disponible
def install_package(package):
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    except Exception as e:
        print(f"Error installing package {package}: {e}")

# Importar librerías con manejo de errores
try:
    import pandas as pd
except ImportError:
    install_package('pandas')
    import pandas as pd

try:
    import requests
except ImportError:
    install_package('requests')
    import requests

try:
    import warnings
except ImportError:
    install_package('warnings')  
    import warnings

try:
    import gzip
except ImportError:
    install_package('gzip')  
    import gzip

try:
    import shutil
except ImportError:
    install_package('shutil')  
    import shutil

try:
    import networkx as nx
except ImportError:
    install_package('networkx')
    import networkx as nx

try:
    import matplotlib.pyplot as plt
except ImportError:
    install_package('matplotlib')
    import matplotlib.pyplot as plt

try:
    import scipy.stats
except ImportError:
    install_package('scipy')
    import scipy.stats

try:
    import numpy as np
except ImportError:
    install_package('numpy')
    import numpy as np

try:
    from collections import defaultdict
except ImportError:
    install_package('collections')  # Parte de la librería estándar
    from collections import defaultdict

try:
    import os
except ImportError:
    install_package('os')  # Parte de la librería estándar
    import os
    
try:
    import subprocess
except ImportError:
    print("subprocess no está instalado. Intentando instalar...")
    # subprocess es parte de la librería estándar, así que esto no debería ser necesario
    install_package('subprocess')
    import subprocess

try:
    import sys
except ImportError:
    print("sys no está instalado. Intentando instalar...")
    # sys es parte de la librería estándar, así que esto no debería ser necesario
    install_package('sys')
    import sys

try:
    import csv
except ImportError:
    print("sys no está instalado. Intentando instalar...")
    # sys es parte de la librería estándar, así que esto no debería ser necesario
    install_package('csv')
    import csv

############################################################
### Llamar los datos a usar (datos sobreexpresados #########
############################################################
    



# Datos de ejemplo si no se puede leer el archivo CSV
try:
    df = pd.read_csv('datos.tsv', sep='\t')
    # Filtrar genes sobreexpresados (logFC > 0 y adj.P.Val < 0.05)
    overexpressed_genes = df[(df['logFC'] > 0) & (df['adj.P.Val'] < 0.05)]['ID'].tolist()
except pd.errors.ParserError as e:
    warnings.warn(f"Error al analizar el archivo: {e}, se usará datos guardados")
    overexpressed_genes = [
        'AT5G25980', 'AT2G25510', 'AT5G26000', 'AT3G22231', 'AT1G56510', 
        'AT5G02500', 'AT3G58780', 'AT2G14610', 'AT3G52700', 'AT5G23010', 'AT3G17390'
    ]

# Definir el organismo (Arabidopsis thaliana) para STRING
organism_taxon_id = 3702

############################################################
### Aplicar análisis funcional a los genes dados ##########
############################################################

print("Se realizará un análisis funcional de los genes sobrexpresados, para contrastar contra los resultados finales")

# STRINGdb API URL and method details
STRING_API_URL = "https://version-11-5.string-db.org/api"
OUTPUT_FORMAT = "json"
METHOD = "enrichment"

# Función para construir la URL de la solicitud
def build_request_url():
    return f"{STRING_API_URL}/{OUTPUT_FORMAT}/{METHOD}"

# Función para construir los parámetros de la solicitud
def build_params(gene_ids, organism_taxon_id):
    return {
        "identifiers": "%0d".join(gene_ids),  # Formatear la lista de identificadores de genes
        "species": organism_taxon_id,  # Usar el ID taxonómico de Arabidopsis thaliana
        "caller_identity": "test_HAB"  # Reemplazar con un identificador propio si es necesario
    }

# Función para hacer la solicitud POST a la API de STRINGdb
def get_functional_enrichment_results(gene_ids, organism_taxon_id):
    request_url = build_request_url()
    params = build_params(gene_ids, organism_taxon_id)
    
    response = requests.post(request_url, data=params)
    
    if response.status_code != 200:
        raise Exception(f"Error en la solicitud: {response.status_code}")
    
    return response.json()

# Función para filtrar los resultados con FDR < 0.01 y categoría "Process"
def filter_significant_processes(data):
    categorias = ["Process", "Component", "Function"]  # Categorías funcionales
    return [
        row for row in data
        if row["category"] in categorias and float(row["fdr"]) < 0.01
    ]

# Función para guardar los resultados en un archivo CSV
def save_results_to_csv(filtered_data, filename="datos_iniciales_string_db_results.csv"):
    # Definir las cabeceras del archivo CSV
    headers = ["Term", "Preferred Names", "FDR", "Category", "Description"]
    
    # Abrir el archivo CSV en modo escritura
    with open(filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(headers)  # Escribir las cabeceras
        
        # Escribir cada fila de datos
        for row in filtered_data:
            term = row["term"]
            preferred_names = ",".join(row["preferredNames"])
            fdr = float(row["fdr"])
            description = row["description"]
            category = row["category"]
            
            writer.writerow([term, preferred_names, fdr, category, description])
    
    print(f"Resultados guardados en {filename}")

# Función principal para ejecutar el análisis y guardar los resultados
def analyze_and_save_functional_enrichment(gene_ids, organism_taxon_id, filename="string_db_results.csv"):
    try:
        # Obtener los resultados de la API
        data = get_functional_enrichment_results(gene_ids, organism_taxon_id)
        
        # Filtrar los resultados para incluir solo GO Biological Processes con FDR < 0.01
        filtered_data = filter_significant_processes(data)
        
        # Guardar los resultados en el archivo CSV
        save_results_to_csv(filtered_data, filename)
    
    except Exception as e:
        print(f"Error durante el análisis: {e}")
        
# Ejecutar el análisis y guardar los resultados en CSV
analyze_and_save_functional_enrichment(overexpressed_genes, organism_taxon_id, "datos_iniciales_string_db_results.csv")


############################################################
####### Descarga de red para la propagación ################
############################################################

# Función para obtener identificador de STRING a partir del ID del gen
def get_string_id(gene_id, organism_taxon_id):
    url = f"https://string-db.org/api/json/get_string_ids?identifiers={gene_id}&species={organism_taxon_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if data:
            # Devolver el primer ID de STRING encontrado, quitando el prefijo '3702.'
            return data[0]['stringId'].split('.')[1] if 'stringId' in data[0] else None
    else:
        warnings.warn(f"No se encontró identificador de STRING para {gene_id}.")
    return None

# Obtener los identificadores de STRING para los genes sobreexpresados
gene_string_ids = {}
for gene in overexpressed_genes:
    string_id = get_string_id(gene, organism_taxon_id)
    if string_id:
        gene_string_ids[gene] = string_id


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
filtered_network_df = network_df[network_df['combined_score'] > 400].copy()  # Hacer una copia explícita

# Eliminar el prefijo '3702.' de las columnas 'protein1' y 'protein2' usando .loc
filtered_network_df['protein1'] = filtered_network_df['protein1'].str.replace('3702.', '', regex=False)
filtered_network_df['protein2'] = filtered_network_df['protein2'].str.replace('3702.', '', regex=False)

# Obtener los STRING IDs del DataFrame filtrado
filtered_proteins = set(filtered_network_df['protein1']).union(set(filtered_network_df['protein2']))

# Comprobar si los STRING IDs están en la red filtrada y emitir advertencias si no se encuentran
for gene, string_id in gene_string_ids.items():
    if string_id not in filtered_proteins:
        warnings.warn(f"El STRING ID {string_id} para el gen {gene} no se encontró en la red filtrada.")


############################################################
################ Propagación de red ########################
############################################################

print("comienzo del DIAMOnD, puede durar unos minutos")
alpha= 1 
n=40
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
G = read_input_from_df(filtered_network_df)

seed_genes = set(gene_string_ids.values())

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
    added_nodes = list(S)  # Incluir los nodos semilla directamente en la lista de nodos seleccionados
    neighbors = {node: set(G.neighbors(node)) for node in G.nodes}  # Diccionario de vecinos de cada nodo
    degrees = dict(G.degree())  # Obtener grados de nodos (número de conexiones para cada nodo)
    cluster_nodes = set(S)  # Conjunto de nodos semilla
    gamma_ln = compute_all_gamma_ln(len(G.nodes))  # Precomputar valores de logaritmo de gamma

    # Iterar hasta seleccionar X nodos (incluyendo los semilla)
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

def graficar_y_guardar_resultados(G, seed_genes, diamond_genes, carpeta_resultados='resultados_propagacion'):
    """
    Grafica la red enriquecida con los nodos semilla y los nodos seleccionados por DIAMOnD.
    Guarda tanto la gráfica como los genes seleccionados en una carpeta especificada.
    """
    # Crear carpeta si no existe
    if not os.path.exists(carpeta_resultados):
        os.makedirs(carpeta_resultados)
    
    # Guardar genes DIAMOnD en un archivo
    with open(os.path.join(carpeta_resultados, 'diamond_genes.txt'), 'w') as f:
        for gene in diamond_genes:
            f.write(f"{gene}\n")
    
    # Configurar la disposición de la red para la gráfica
    pos = nx.spring_layout(G)  # Disposición de los nodos en el espacio (layout de la red)
    plt.figure(figsize=(500, 500))  # Tamaño de la figura

    # Nodos de semillas (color azul claro)
    nx.draw_networkx_nodes(G, pos, nodelist=seed_genes, node_color='lightblue', node_size=500, label="Seed Genes")

    # Nodos DIAMOnD (color naranja)
    nx.draw_networkx_nodes(G, pos, nodelist=diamond_genes, node_color='orange', node_size=100, label="DIAMOnD Genes")

    # Enlaces entre nodos (semillas y DIAMOnD)
    nx.draw_networkx_edges(G, pos, edgelist=[(u, v) for u, v in G.edges() if u in (seed_genes | set(diamond_genes)) and v in (seed_genes | set(diamond_genes))], alpha=0.5)

    # Etiquetas de los nodos
    nx.draw_networkx_labels(G, pos, font_size=8)

    # Mostrar leyenda y título
    plt.legend(loc="best")
    plt.title("Red de Genes Enriquecida usando DIAMOnD")

    # Guardar la gráfica en la carpeta
    plt.savefig(os.path.join(carpeta_resultados, 'red_enriquecida.png'))
    plt.close()  # Cerrar la figura para liberar memoria

# ------------------- EJECUCIÓN -------------------

# Realizar la iteración DIAMOnD para obtener los genes seleccionados
diamond_genes = diamond_iteration_of_first_X_nodes(G, seed_genes, n, alpha)


# Ruta donde deseas guardar los genes y la grafica
carpeta_resultados = 'resultados_propagacion'

# Graficar la red enriquecida
graficar_y_guardar_resultados(G, seed_genes, diamond_genes, carpeta_resultados)

"""
def random_walk_with_restart(G, seed_nodes, restart_prob=0.85, max_iter=100):
    #Algoritmo de Random Walk with Restart
    nodes = list(G.nodes)
    N = len(nodes)
    r = {n: 1 / N for n in nodes}  # Distribución de probabilidad inicial uniforme
    p = {n: 1 if n in seed_nodes else 0 for n in nodes}  # Reiniciar en nodos semilla
    
    for _ in range(max_iter):
        r_next = {n: (1 - restart_prob) * sum(r[neighbor] / G.degree(neighbor) for neighbor in G.neighbors(n)) + restart_prob * p[n] for n in nodes}
        r = r_next

    # Ordenar nodos por importancia
    return sorted(r.items(), key=lambda x: x[1], reverse=True)

def heat_diffusion(G, seed_nodes, diffusion_time=10):
    #Algoritmo de difusión de calor
    heat = {n: 1 if n in seed_nodes else 0 for n in G.nodes}  # Calor inicial en nodos semilla
    for _ in range(diffusion_time):
        new_heat = {n: sum(heat[neighbor] / G.degree(neighbor) for neighbor in G.neighbors(n)) for n in G.nodes}
        heat = new_heat
    return sorted(heat.items(), key=lambda x: x[1], reverse=True)

def personalized_pagerank(G, seed_nodes, alpha=0.85):
    personalization = {n: 1 if n in seed_nodes else 0 for n in G.nodes}
    return nx.pagerank(G, alpha=alpha, personalization=personalization)


"""
