import os, re
from Bio import Entrez, SeqIO

# Função para ler o arquivo Input
def ler_arquivo(nome_arquivo):
    try:
        with open(nome_arquivo, 'r') as arquivo:
            linhas = arquivo.readlines()
            valores = [linha.strip() for linha in linhas]
            return valores
    except FileNotFoundError:
        print("Arquivo não encontrado.")
        return []

# Organização dos dados
def organizar_dados(genoma_id):

    arquivo = open('resultados_localizacao_genes/Genes_'+str(genoma_id)+'.txt','w')
    arquivo.write('#START\tSTOP\tLABEL\tCOLOUR\n')
    
    for rec in SeqIO.parse('../Banco_de_dados/banco_de_dados_genomas/genoma_'+str(genoma_id)+'.gb', 'genbank'):
        
        for feature in rec.features:

            if feature.type == 'gene':

                for record in SeqIO.parse('../Banco_de_dados/Banco_de_dados_AgroBioTech/genes_pgp_referencias/genes_pgp_referencias_'+str(genoma_id)+'.fasta', 'fasta'):

                    desc = record.description
                    indice = desc.find('locus_tag:')

                    if not indice == -1:

                        locus = desc[(indice+11):(indice+22)]

                        if locus == feature.qualifiers['locus_tag'][0]:

                            local = str(feature.location)
                            regex = re.search('(\d+)\D(\d+)', local).groups()
                            start = regex[0]
                            end = regex[1]

                            arquivo.write(start+'\t'+end+'\t'+str(record.name)+'\tblack\n')

    arquivo.close()

# Leitura do Input do Script
if len(sys.argv) != 2:

    print("Uso: python blastp.py input.txt")

else:

    nome_arquivo = sys.argv[1]
    valores = ler_arquivo(nome_arquivo)

    if valores:

        for genoma_id in valores:

            organizar_dados(genoma_id)

#FIM
